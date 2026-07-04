"""Build HMQC + NOESY peak lists for MAUS from a BMRB chemical-shift file
and a PDB structure.

This replaces the previous circular setup (NOEs simulated from the same
coordinates *and indexed identically to the methyls*).  Here the two peak
lists are the only things MAUS sees, and the link back to the structure is
made **through chemical-shift frequency matching**, not through identity —
so shift degeneracy produces genuine assignment ambiguity, exactly as in the
paper.

Outputs
-------
hmqc.tsv (input) : label  H_ppm  C_ppm  res_type
    What MAUS actually sees.  `label` is an anonymous peak id (P1, P2, ...) that
    leaks nothing about the answer.  (1H, 13C) = data-graph node coordinates;
    residue *type* is known (from the labeling / spectral region).

hmqc_true.tsv (truth key) : label  H_ppm  C_ppm  res_type  True
    Same rows as hmqc.tsv plus a `True` column giving the real methyl identity
    for each anonymous label.  Used only for scoring, never as input.

noesy.tsv : label  C1  C2  H2  intensity
    Methyl-methyl 3D (H)CCH NOESY cross peaks.  The observed methyl is (H2, C2);
    the NOE partner contributes carbon C1 only (its proton is not in the peak).
    intensity is the cross-peak height (~ r^-6); read but not used for matching.
    Both directions are emitted for each close pair (symmetric NOESY).  MAUS
    matches the observed methyl by (H2,C2) and the partner by carbon C1.

hmbc.tsv : label  C1  C2  H2
    Optional HMBC-HMQC geminal links: one row per Leu/Val residue, same layout
    as noesy.tsv.  The observed methyl is (H2, C2); its geminal partner
    contributes carbon C1.  MAUS forces the pair onto a geminal structure edge.
"""

from __future__ import annotations

import argparse
import math
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Tuple

THREE_TO_ONE = {'ALA': 'A', 'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'THR': 'T', 'VAL': 'V'}
# residue one-letter -> [(methyl carbon, geminal-partner carbon or None)]
METHYL_ATOMS = {
  'A': [('CB', None)],
  'I': [('CD1', None)],
  'L': [('CD1', 'CD2'), ('CD2', 'CD1')],
  'M': [('CE', None)],
  'T': [('CG2', None)],
  'V': [('CG1', 'CG2'), ('CG2', 'CG1')],
}


def methyl_protons(carbon: str) -> List[str]:
  """Methyl proton names for a methyl carbon: CD1 -> HD11/HD12/HD13, CB -> HB1..."""
  stem = 'H' + carbon[1:]
  return [f'{stem}{i}' for i in (1, 2, 3)]


def parse_bmrb_shifts(star_lines) -> Dict[Tuple[int, str, str], float]:
  """Return {(seq_id, comp_3letter, atom): shift}. Parses the
  _Atom_chem_shift loop rows (whitespace-delimited NMR-STAR 3.1)."""
  shifts: Dict[Tuple[int, str, str], float] = {}
  for line in star_lines:
    tok = line.split()
    # loop row layout: ID . EAsm Ent CompIdx Seq Comp Atom Type Iso Val ...
    if len(tok) < 11:
      continue
    comp = tok[6]
    if comp not in THREE_TO_ONE:
      continue
    try:
      seq = int(tok[5])
      val = float(tok[10])
    except ValueError:
      continue
    atom = tok[7]
    shifts[(seq, comp, atom)] = val
  return shifts


def parse_structure_coords(pdb_lines):
  """{(one_letter, resi): {atom: (x,y,z)}} for methyl-bearing residues."""
  coords: Dict[Tuple[str, int], Dict[str, Tuple[float, float, float]]] = {}
  for line in pdb_lines:
    if not line.startswith('ATOM'):
      continue
    resn = line[17:20].strip()
    if resn not in THREE_TO_ONE:
      continue
    if line[21] not in (' ', 'A'):
      continue
    try:
      resi = int(line[22:26])
    except ValueError:
      continue
    atom = line[12:16].strip()
    coords.setdefault((THREE_TO_ONE[resn], resi), {})[atom] = (
      float(line[30:38]), float(line[38:46]), float(line[46:54]))
  return coords


def build_methyls(coords, shifts):
  """One record per methyl that has BOTH structure coords and BMRB shifts.
  Returns list of dicts with label, res_type, geminal carbon, xyz, H, C."""
  three = {v: k for k, v in THREE_TO_ONE.items()}
  out = []
  for (one, resi) in sorted(coords, key=lambda k: (k[1], k[0])):
    comp = three[one]
    for carbon, gem in METHYL_ATOMS[one]:
      if carbon not in coords[(one, resi)]:
        continue
      c_shift = shifts.get((resi, comp, carbon))
      h_vals = [shifts.get((resi, comp, hp)) for hp in methyl_protons(carbon)]
      h_vals = [v for v in h_vals if v is not None]
      if c_shift is None or not h_vals:
        continue  # unassigned in BMRB -> no HMQC peak
      out.append({
        'label': f'{one}{resi}{carbon}',
        'res_type': one,
        'carbon': carbon,
        'geminal': gem or '',
        'xyz': coords[(one, resi)][carbon],
        'H': round(sum(h_vals) / len(h_vals), 3),
        'C': round(c_shift, 3),
      })
  return out


def main(argv=None):
  ap = argparse.ArgumentParser()
  ap.add_argument('pdb')
  ap.add_argument('bmrb', help='NMR-STAR chemical shift file (bmrXXXX_3.str)')
  ap.add_argument('--noe-long', type=float, default=8.0, help='NOE distance cutoff (A)')
  ap.add_argument('--keep-k', type=int, default=12, help='nearest-K NOE partners per methyl')
  ap.add_argument('--hmqc-out', default='examples/mbp/hmqc.tsv')
  ap.add_argument('--truth-out', default='examples/mbp/hmqc_true.tsv')
  ap.add_argument('--noesy-out', default='examples/mbp/noesy.tsv')
  ap.add_argument('--hmbc-out', default='examples/mbp/hmbc.tsv')
  args = ap.parse_args(argv)

  coords = parse_structure_coords(Path(args.pdb).read_text().splitlines())
  shifts = parse_bmrb_shifts(Path(args.bmrb).read_text().splitlines())
  methyls = build_methyls(coords, shifts)

  # --- HMQC peak lists ---
  # input: anonymous label (P1..) so nothing about the answer leaks;
  # truth key: same rows + a `True` column mapping label -> real methyl.
  in_lines = ['label\tH_ppm\tC_ppm\tres_type']
  true_lines = ['label\tH_ppm\tC_ppm\tres_type\tTrue']
  for i, m in enumerate(methyls, 1):
    pid = f'P{i}'
    row = f'{pid}\t{m["H"]:.3f}\t{m["C"]:.3f}\t{m["res_type"]}'
    in_lines.append(row)
    true_lines.append(f'{row}\t{m["label"]}')
  Path(args.hmqc_out).write_text('\n'.join(in_lines) + '\n')
  Path(args.truth_out).write_text('\n'.join(true_lines) + '\n')

  # --- NOESY peak list (3D (H)CCH): label C1 C2 H2 intensity ---
  #   observed methyl = (H2, C2); partner methyl contributes carbon C1 only.
  #   intensity ~ r^-6 (NOE buildup), scaled to readable numbers.
  #   Symmetric NOESY: emit both directions for each close pair.
  noe = {}  # unordered (i,j), i<j -> distance, within long cut
  for a in range(len(methyls)):
    dists = sorted(
      (math.dist(methyls[a]['xyz'], methyls[b]['xyz']), b)
      for b in range(len(methyls)) if b != a)
    for d, b in dists[:args.keep_k]:
      if d > args.noe_long:
        break
      noe[(min(a, b), max(a, b))] = d

  noesy_lines = ['label\tC1\tC2\tH2\tintensity']
  k = 0
  for (a, b), d in sorted(noe.items()):
    intensity = round(d ** -6 * 1e6, 1)      # r^-6 buildup, arbitrary scale
    for obs, par in ((a, b), (b, a)):        # both directions
      k += 1
      mo, mp = methyls[obs], methyls[par]
      noesy_lines.append(
        f'X{k}\t{mp["C"]:.3f}\t{mo["C"]:.3f}\t{mo["H"]:.3f}\t{intensity}')
  Path(args.noesy_out).write_text('\n'.join(noesy_lines) + '\n')

  # --- HMBC-HMQC geminal links (3D (H)CCH): label C1 C2 H2 ---
  #   observed methyl = (H2, C2); geminal partner contributes carbon C1 only.
  by_label = {m['label']: m for m in methyls}
  hmbc_pairs = []
  for m in methyls:
    if not m['geminal']:
      continue
    if m['carbon'] > m['geminal']:      # emit once, from the lower carbon name
      continue
    partner = by_label.get(m['label'][:-len(m['carbon'])] + m['geminal'])
    if partner is not None:
      hmbc_pairs.append((m, partner))
  hmbc_lines = ['label\tC1\tC2\tH2']
  for k, (mo, mp) in enumerate(hmbc_pairs, 1):
    # observed methyl mo: proton H2 + carbon C2; partner mp: carbon C1
    hmbc_lines.append(
      f'B{k}\t{mp["C"]:.3f}\t{mo["C"]:.3f}\t{mo["H"]:.3f}')
  Path(args.hmbc_out).write_text('\n'.join(hmbc_lines) + '\n')

  # --- degeneracy report ---
  seen: Dict[Tuple[float, float], int] = {}
  for m in methyls:
    key = (round(m['H'], 2), round(m['C'], 2))
    seen[key] = seen.get(key, 0) + 1
  degenerate = sum(1 for m in methyls
                   if seen[(round(m['H'], 2), round(m['C'], 2))] > 1)
  print(f'methyls with BMRB shift + structure = {len(methyls)}')
  print(f'HMQC input peaks    = {len(methyls)}  -> {args.hmqc_out}')
  print(f'HMQC truth key      = {len(methyls)}  -> {args.truth_out}')
  print(f'NOESY cross peaks   = {2*len(noe)}  ({len(noe)} pairs x2 directions)  -> {args.noesy_out}')
  print(f'HMBC geminal links  = {len(hmbc_pairs)}  -> {args.hmbc_out}')
  print(f'near-degenerate HMQC peaks (>=2 within 0.01/0.01 ppm bin) = {degenerate}')
  return 0


if __name__ == '__main__':
  raise SystemExit(main())
