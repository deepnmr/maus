"""MAUS — Methyl Assignments Using Satisfiability (clean-room reimplementation).

Based on Nerli, De Paula, McShan & Sgourakis, "Backbone-independent NMR
resonance assignments of methyl probes in large proteins", Nat. Commun. 12:691
(2021).  https://doi.org/10.1038/s41467-021-20984-0

Inputs are *peak lists*, as in a real experiment — nothing about the answer is
baked into the indexing:

  * HMQC peak list  (data-graph nodes) : each 2D methyl peak carries its
    (1H, 13C) shift and residue *type*.  Its structural identity is unknown and
    is what MAUS recovers.
  * NOESY peak list (data-graph edges) : 3D (H)CCH methyl-methyl cross peaks
    `label C1 C2 H2`.  The observed methyl is matched by (H2, C2); the partner
    by carbon C1 only.  Shift degeneracy — especially the carbon-only partner —
    yields genuine, irreducible ambiguity.
  * Structure graph G : one node per methyl carbon from the PDB; edges
    classified geminal / short (< short_cut) / long (< long_cut).

Hard constraints (no scoring): every HMQC peak maps to exactly one methyl of the
matching type, the map is injective, and every NOESY cross peak that can be
assigned to a definite pair of HMQC peaks must map onto a G edge within long_cut.
For each peak the set of valid methyls (those in >=1 satisfying assignment) is
enumerated with the solver's assumption interface — the paper's iterative ansatz
— so the ground truth is provably never excluded.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from pysat.solvers import Glucose3


THREE_TO_ONE = {'ALA': 'A', 'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'THR': 'T', 'VAL': 'V'}
# residue type -> [(methyl carbon atom, geminal-partner atom or None)]
METHYL_ATOMS = {
  'A': [('CB', None)],
  'I': [('CD1', None)],
  'L': [('CD1', 'CD2'), ('CD2', 'CD1')],
  'M': [('CE', None)],
  'T': [('CG2', None)],
  'V': [('CG1', 'CG2'), ('CG2', 'CG1')],
}


@dataclass(frozen=True)
class Methyl:
  index: int
  label: str          # e.g. L135CD1
  res_type: str       # one-letter
  res_num: int
  atom: str
  coord: Tuple[float, float, float]
  geminal_atom: str   # partner methyl atom name, or '' if none


def parse_structure(pdb_lines, labeling: Dict[str, List[Tuple[str, Optional[str]]]]) -> List[Methyl]:
  coords: Dict[Tuple[str, int], Dict[str, Tuple[float, float, float]]] = {}
  for line in pdb_lines:
    if not line.startswith('ATOM'):
      continue
    resn = line[17:20].strip()
    if resn not in THREE_TO_ONE:
      continue
    chain = line[21]
    if chain not in (' ', 'A'):
      continue
    try:
      resi = int(line[22:26])
    except ValueError:
      continue
    atom = line[12:16].strip()
    coords.setdefault((THREE_TO_ONE[resn], resi), {})[atom] = (
      float(line[30:38]), float(line[38:46]), float(line[46:54])
    )

  methyls: List[Methyl] = []
  for (one, resi) in sorted(coords, key=lambda k: (k[1], k[0])):
    for atom, gem in labeling.get(one, []):
      if atom not in coords[(one, resi)]:
        continue
      methyls.append(Methyl(
        index=len(methyls),
        label=f'{one}{resi}{atom}',
        res_type=one,
        res_num=resi,
        atom=atom,
        coord=coords[(one, resi)][atom],
        geminal_atom=gem or '',
      ))
  return methyls


def build_structure_graph(methyls: List[Methyl], short_cut: float, long_cut: float):
  """Return sets of edge classes: geminal / short / long (all symmetric)."""
  gem, short, long = set(), set(), set()
  for a, b in combinations(methyls, 2):
    same_res = a.res_num == b.res_num and a.res_type == b.res_type
    if same_res and (a.atom == b.geminal_atom):
      gem.add((a.index, b.index)); gem.add((b.index, a.index))
      continue
    d = math.dist(a.coord, b.coord)
    if d < short_cut:
      short.add((a.index, b.index)); short.add((b.index, a.index))
    elif d <= long_cut:
      long.add((a.index, b.index)); long.add((b.index, a.index))
  return gem, short, long


@dataclass(frozen=True)
class Peak:
  index: int
  peak_id: str
  res_type: str
  h_ppm: float
  c_ppm: float


def load_hmqc(path: str) -> List[Peak]:
  """Input HMQC peak list: label  H_ppm  C_ppm  res_type.  `label` is an
  anonymous peak id and carries nothing about the answer."""
  peaks: List[Peak] = []
  for line in Path(path).read_text().splitlines():
    if not line.strip() or line.startswith('#') or line.startswith('label'):
      continue
    label, h, c, rtype = line.split('\t')
    peaks.append(Peak(index=len(peaks), peak_id=label, res_type=rtype.strip(),
                      h_ppm=float(h), c_ppm=float(c)))
  return peaks


def load_truth(path: str) -> Dict[str, str]:
  """Truth key: label  H_ppm  C_ppm  res_type  True.  Returns {label: True}."""
  truth: Dict[str, str] = {}
  for line in Path(path).read_text().splitlines():
    if not line.strip() or line.startswith('#') or line.startswith('label'):
      continue
    cols = line.split('\t')
    truth[cols[0]] = cols[4].strip()
  return truth


def load_ccch(path: str, id_prefix: str) -> List[Tuple[float, float, float]]:
  """Load a 3D (H)CCH-style peak list: label C1 C2 H2.  The observed methyl is
  (H2, C2); the partner methyl contributes carbon C1 only.  Returns [(c1,c2,h2)].
  Used for both NOESY and HMBC-HMQC (same layout)."""
  rows = []
  for line in Path(path).read_text().splitlines():
    if not line.strip() or line.startswith('#') or line.startswith('label'):
      continue
    _label, c1, c2, h2 = line.split('\t')[:4]
    rows.append((float(c1), float(c2), float(h2)))
  return rows


def load_noesy(path: str):
  return load_ccch(path, 'X')


def load_hmbc(path: str):
  return load_ccch(path, 'B')


def _match_rows(peaks: List[Peak], rows, tol_h: float, tol_c: float):
  """Resolve each row to a pair of HMQC peaks: the observed methyl by both
  dimensions (H2, C2), the partner by carbon C1 only (its proton is not in the
  row).  A firm edge is kept only when BOTH endpoints resolve to a unique HMQC
  peak; otherwise it is ambiguous and dropped (conservative — never excludes a
  valid assignment).  Returns (edges, stats)."""
  def cand_hc(h, c):
    return [p.index for p in peaks
            if abs(p.h_ppm - h) <= tol_h and abs(p.c_ppm - c) <= tol_c]

  def cand_c(c):
    return [p.index for p in peaks if abs(p.c_ppm - c) <= tol_c]

  edges = set()
  firm = ambiguous = unmatched = 0
  for (c1, c2, h2) in rows:
    a = cand_hc(h2, c2)     # observed methyl (proton + carbon)
    b = cand_c(c1)          # partner methyl (carbon only)
    if not a or not b:
      unmatched += 1
      continue
    if len(a) == 1 and len(b) == 1 and a[0] != b[0]:
      i, j = a[0], b[0]
      edges.add((min(i, j), max(i, j)))
      firm += 1
    else:
      ambiguous += 1
  return edges, {'firm': firm, 'ambiguous': ambiguous, 'unmatched': unmatched}


def match_noe(peaks, rows, tol_h, tol_c):
  return _match_rows(peaks, rows, tol_h, tol_c)


def match_hmbc(peaks, rows, tol_h, tol_c):
  return _match_rows(peaks, rows, tol_h, tol_c)


class MAUS:
  def __init__(self, methyls, peaks, gem, short_g, long_g, noe_edges,
               gem_links=None):
    self.methyls = methyls
    self.peaks = peaks
    self.gem, self.short_g, self.long_g = gem, short_g, long_g
    self.noe_edges = noe_edges            # NOESY firm edges (single class)
    self.gem_links = gem_links or set()   # HMBC geminal same-residue links
    self.sites_by_type: Dict[str, List[int]] = {}
    for m in methyls:
      self.sites_by_type.setdefault(m.res_type, []).append(m.index)
    # SAT variable ids: x[(peak_i, methyl_g)]
    self.var: Dict[Tuple[int, int], int] = {}
    self.domain: Dict[int, List[int]] = {}
    nid = 1
    for p in peaks:
      dom = self.sites_by_type.get(p.res_type, [])
      self.domain[p.index] = dom
      for g in dom:
        self.var[(p.index, g)] = nid
        nid += 1
    self.n_vars = nid - 1

  def _base_clauses(self) -> List[List[int]]:
    clauses: List[List[int]] = []
    # (1) each peak assigned to exactly one methyl of its type
    for p in self.peaks:
      lits = [self.var[(p.index, g)] for g in self.domain[p.index]]
      if not lits:
        continue
      clauses.append(lits)                                  # at least one
      for a, b in combinations(lits, 2):
        clauses.append([-a, -b])                            # at most one
    # (2) each methyl used by at most one peak
    peaks_of: Dict[int, List[int]] = {}
    for p in self.peaks:
      for g in self.domain[p.index]:
        peaks_of.setdefault(g, []).append(self.var[(p.index, g)])
    for g, lits in peaks_of.items():
      for a, b in combinations(lits, 2):
        clauses.append([-a, -b])
    # (3) every firm NOE edge must map onto a structure edge within long_cut
    #     (single mixing class: geminal / short / long all acceptable).
    allowed = self.gem | self.short_g | self.long_g
    for (i, j) in self.noe_edges:
      di, dj = self.domain[i], self.domain[j]
      for gi in di:
        for gj in dj:
          if gi == gj:
            continue
          if (gi, gj) not in allowed:
            clauses.append([-self.var[(i, gi)], -self.var[(j, gj)]])
    # (4) optional HMBC geminal links: the two linked peaks must map to the two
    #     geminal methyls of one residue -> only geminal G-edges allowed.
    for (i, j) in self.gem_links:
      di, dj = self.domain[i], self.domain[j]
      for gi in di:
        for gj in dj:
          if gi == gj:
            continue
          if (gi, gj) not in self.gem:
            clauses.append([-self.var[(i, gi)], -self.var[(j, gj)]])
    return clauses

  def solve_options(self) -> Dict[int, List[int]]:
    # Build the CNF once; enumerate per-peak valid methyls incrementally with
    # the solver's assumption interface (the paper's iterative ansatz). A methyl
    # g is a valid option for peak i iff the CNF is still satisfiable when
    # x(i,g) is asserted true.
    base = self._base_clauses()
    solver = Glucose3(bootstrap_with=base)
    options: Dict[int, List[int]] = {}
    for p in self.peaks:
      valid = [g for g in self.domain[p.index]
               if solver.solve(assumptions=[self.var[(p.index, g)]])]
      options[p.index] = valid
    solver.delete()
    return options


def parse_labeling(spec: str) -> Dict[str, List[Tuple[str, Optional[str]]]]:
  labeling: Dict[str, List[Tuple[str, Optional[str]]]] = {}
  for chunk in spec.split(';'):
    chunk = chunk.strip()
    if not chunk:
      continue
    one = chunk.split(',')[0].strip()
    labeling[one] = METHYL_ATOMS.get(one, [])
  return labeling


def main(argv=None):
  ap = argparse.ArgumentParser(description='MAUS: SAT-based methyl assignment (clean-room).')
  ap.add_argument('pdb')
  ap.add_argument('hmqc', help='HMQC peak list TSV: label H_ppm C_ppm res_type')
  ap.add_argument('noesy', help='NOESY peak list TSV: label C1 C2 H2')
  ap.add_argument('--hmbc', default=None,
                  help='optional HMBC-HMQC peak list TSV (label C1 C2 H2): '
                       'geminal same-residue links for Leu/Val')
  ap.add_argument('--truth', default=None,
                  help='truth key TSV (label ... True) for scoring only')
  ap.add_argument('--short-cut', type=float, default=6.0, help='structure short-range cutoff (A)')
  ap.add_argument('--long-cut', type=float, default=10.0, help='structure long-range cutoff (A)')
  ap.add_argument('--tol-h', type=float, default=0.02, help='1H match tolerance (ppm)')
  ap.add_argument('--tol-c', type=float, default=0.20, help='13C match tolerance (ppm)')
  ap.add_argument('--labeling', default='A;I;L;M;T;V')
  ap.add_argument('--out', default=None, help='write per-peak options TSV here')
  args = ap.parse_args(argv)

  labeling = parse_labeling(args.labeling)
  methyls = parse_structure(Path(args.pdb).read_text().splitlines(), labeling)
  peaks = load_hmqc(args.hmqc)
  crosses = load_noesy(args.noesy)

  gem, short_g, long_g = build_structure_graph(methyls, args.short_cut, args.long_cut)
  noe_edges, nstat = match_noe(peaks, crosses, args.tol_h, args.tol_c)
  gem_links, hstat = set(), None
  if args.hmbc:
    gem_links, hstat = match_hmbc(peaks, load_hmbc(args.hmbc), args.tol_h, args.tol_c)
  maus = MAUS(methyls, peaks, gem, short_g, long_g, noe_edges, gem_links)
  options = maus.solve_options()

  label_by_index = {m.index: m.label for m in methyls}
  truth = load_truth(args.truth) if args.truth else {}

  def truth_of(p):
    return truth.get(p.peak_id, '')

  if args.out:
    with open(args.out, 'w') as f:
      f.write('label\tres_type\tn_options\toptions\ttruth\ttruth_in_set\n')
      for p in peaks:
        labels = [label_by_index[g] for g in options[p.index]]
        t = truth_of(p)
        in_set = '' if not truth else int(t in labels)
        f.write(f'{p.peak_id}\t{p.res_type}\t{len(labels)}\t'
                f'{",".join(labels)}\t{t}\t{in_set}\n')

  n = len(peaks)
  unique = amb23 = amb_more = unassigned = 0
  for p in peaks:
    k = len(options[p.index])
    if k == 0:
      unassigned += 1
    elif k == 1:
      unique += 1
    elif k <= 3:
      amb23 += 1
    else:
      amb_more += 1

  print(f'methyls(G nodes)={len(methyls)}  HMQC peaks={n}  NOESY cross peaks={len(crosses)}')
  print(f'G edges: geminal={len(gem)//2} short={len(short_g)//2} long={len(long_g)//2}')
  print(f'NOE match (tol H+-{args.tol_h}/C+-{args.tol_c}): '
        f'firm={nstat["firm"]} ambiguous(dropped)={nstat["ambiguous"]} unmatched={nstat["unmatched"]}')
  print(f'firm NOE data edges: {len(noe_edges)}')
  if hstat is not None:
    print(f'HMBC geminal links: firm={hstat["firm"]} '
          f'ambiguous={hstat["ambiguous"]} unmatched={hstat["unmatched"]}')
  print(f'unique(1 option)      = {unique}/{n}')
  print(f'ambiguous(2-3 options)= {amb23}/{n}')
  print(f'ambiguous(>3 options) = {amb_more}/{n}')
  print(f'unassigned            = {unassigned}/{n}')

  if truth:
    correct_in_set = sum(1 for p in peaks
                         if truth_of(p) in [label_by_index[g] for g in options[p.index]])
    unique_correct = sum(1 for p in peaks
                         if len(options[p.index]) == 1
                         and label_by_index[options[p.index][0]] == truth_of(p))
    denom = unique or 1
    print(f'truth in option set   = {correct_in_set}/{n} = {100*correct_in_set/n:.1f}%  (error rate {100*(n-correct_in_set)/n:.1f}%)')
    print(f'unique calls correct  = {unique_correct}/{unique} = {100*unique_correct/denom:.1f}%')
  else:
    print('(no --truth given; scoring skipped)')
  return 0


if __name__ == '__main__':
  raise SystemExit(main())
