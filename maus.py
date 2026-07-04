"""MAUS — Methyl Assignments Using Satisfiability (clean-room reimplementation).

Based on Nerli, De Paula, McShan & Sgourakis, "Backbone-independent NMR
resonance assignments of methyl probes in large proteins", Nat. Commun. 12:691
(2021).  https://doi.org/10.1038/s41467-021-20984-0

Unlike MAGIC (a scoring/exhaustive-search method), MAUS casts methyl assignment
as a *subgraph isomorphism* solved with a SAT solver:

  * Structure graph G: one node per methyl carbon; edges classified as
    geminal / short-range (< short_cut) / long-range (< long_cut) from the PDB.
  * Data graph H: 2D reference peaks (residue type known) with symmetric NOE
    edges, split into short- and long-mixing-time classes.
  * Hard constraints (no scoring): every peak maps to exactly one methyl of the
    matching residue type, the map is injective, geminal edges are preserved,
    and every NOE edge must map onto a G edge of a compatible distance class.

For each peak the set of *valid* methyls (those appearing in at least one
satisfying assignment) is enumerated with the solver's assumption interface —
the paper's iterative ansatz — giving 1 / 2-3 / >3 option counts and, because
the enumeration is exact, a guarantee that the ground truth is never excluded
when the inputs are consistent.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Tuple

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


def parse_structure(pdb_lines, labeling: Dict[str, List[Tuple[str, str]]]) -> List[Methyl]:
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
  """Return dicts of edge classes: geminal / short / long (all symmetric)."""
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


@dataclass
class Peak:
  index: int
  peak_id: str
  res_type: str
  truth_label: str


def simulate_noe(methyls: List[Methyl], short_cut: float, long_cut: float, keep_k: int):
  """Simulate symmetric short/long NOE edges from the structure (sparse: the
  KEEP_K nearest partners per methyl within long_cut), mimicking real data."""
  short_e, long_e = set(), set()
  for a in methyls:
    dists = sorted((math.dist(a.coord, b.coord), b.index) for b in methyls if b.index != a.index)
    for d, j in dists[:keep_k]:
      if d > long_cut:
        break
      pair = (a.index, j)
      if d < short_cut:
        short_e.add(pair)
      else:
        long_e.add(pair)
  # symmetrize (keep only edges seen in both directions — MAUS symmetrization)
  short_sym = {(i, j) for (i, j) in short_e if (j, i) in short_e}
  long_sym = {(i, j) for (i, j) in long_e if (j, i) in long_e}
  return short_sym, long_sym


class MAUS:
  def __init__(self, methyls, peaks, gem, short_g, long_g, short_noe, long_noe):
    self.methyls = methyls
    self.peaks = peaks
    self.gem, self.short_g, self.long_g = gem, short_g, long_g
    self.short_noe, self.long_noe = short_noe, long_noe
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
    # (3) every NOE edge must map onto a compatible structure edge.
    #     short NOE -> geminal or short G-edge; long NOE -> any G-edge (<= long_cut)
    short_ok = self.gem | self.short_g
    long_ok = self.gem | self.short_g | self.long_g
    for noe, allowed in ((self.short_noe, short_ok), (self.long_noe, long_ok)):
      for (i, j) in noe:
        if i >= j:            # each unordered edge once
          continue
        di, dj = self.domain[i], self.domain[j]
        for gi in di:
          for gj in dj:
            if gi == gj:
              continue
            if (gi, gj) not in allowed:
              # forbid this simultaneous assignment
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


def parse_labeling(spec: str) -> Dict[str, List[Tuple[str, str]]]:
  labeling: Dict[str, List[Tuple[str, str]]] = {}
  for chunk in spec.split(';'):
    chunk = chunk.strip()
    if not chunk:
      continue
    tokens = [t.strip() for t in chunk.split(',')]
    one = tokens[0]
    labeling[one] = METHYL_ATOMS.get(one, [])
  return labeling


def main(argv=None):
  ap = argparse.ArgumentParser(description='MAUS: SAT-based methyl assignment (clean-room).')
  ap.add_argument('pdb')
  ap.add_argument('peaks', help='TSV: peak_id  res_type  truth_label')
  ap.add_argument('--short-cut', type=float, default=6.0)
  ap.add_argument('--long-cut', type=float, default=10.0)
  ap.add_argument('--noe-short', type=float, default=6.0)
  ap.add_argument('--noe-long', type=float, default=8.0)
  ap.add_argument('--keep-k', type=int, default=8)
  ap.add_argument('--labeling', default='A;I;L;M;T;V')
  ap.add_argument('--out', default=None, help='write per-peak options TSV here')
  args = ap.parse_args(argv)

  labeling = parse_labeling(args.labeling)
  methyls = parse_structure(Path(args.pdb).read_text().splitlines(), labeling)
  peaks = []
  for line in Path(args.peaks).read_text().splitlines():
    if not line.strip() or line.startswith('#'):
      continue
    pid, rtype, truth = line.split('\t')
    peaks.append(Peak(index=len(peaks), peak_id=pid, res_type=rtype, truth_label=truth))

  gem, short_g, long_g = build_structure_graph(methyls, args.short_cut, args.long_cut)
  short_noe, long_noe = simulate_noe(methyls, args.noe_short, args.noe_long, args.keep_k)
  maus = MAUS(methyls, peaks, gem, short_g, long_g, short_noe, long_noe)
  options = maus.solve_options()

  label_by_index = {m.index: m.label for m in methyls}
  if args.out:
    with open(args.out, 'w') as f:
      f.write('peak_id\tres_type\tn_options\toptions\ttruth\ttruth_in_set\n')
      for p in peaks:
        labels = [label_by_index[g] for g in options[p.index]]
        f.write(f'{p.peak_id}\t{p.res_type}\t{len(labels)}\t'
                f'{",".join(labels)}\t{p.truth_label}\t'
                f'{int(p.truth_label in labels)}\n')
  n = len(peaks)
  unique = amb23 = amb_more = unassigned = correct_in_set = 0
  for p in peaks:
    opts = options[p.index]
    labels = [label_by_index[g] for g in opts]
    if len(opts) == 0:
      unassigned += 1
    elif len(opts) == 1:
      unique += 1
    elif len(opts) <= 3:
      amb23 += 1
    else:
      amb_more += 1
    if p.truth_label in labels:
      correct_in_set += 1
  # accuracy of the unique calls
  unique_correct = sum(1 for p in peaks
                       if len(options[p.index]) == 1
                       and label_by_index[options[p.index][0]] == p.truth_label)

  print(f'methyls(G nodes)={len(methyls)}  peaks={n}')
  print(f'G edges: geminal={len(gem)//2} short={len(short_g)//2} long={len(long_g)//2}')
  print(f'NOE edges: short={len(short_noe)//2} long={len(long_noe)//2}')
  print(f'unique(1 option)      = {unique}/{n}')
  print(f'ambiguous(2-3 options)= {amb23}/{n}')
  print(f'ambiguous(>3 options) = {amb_more}/{n}')
  print(f'unassigned            = {unassigned}/{n}')
  print(f'truth in option set   = {correct_in_set}/{n} = {100*correct_in_set/n:.1f}%  (error rate {100*(n-correct_in_set)/n:.1f}%)')
  denom = unique or 1
  print(f'unique calls correct  = {unique_correct}/{unique} = {100*unique_correct/denom:.1f}%')
  return 0


if __name__ == '__main__':
  raise SystemExit(main())
