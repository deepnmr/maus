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
hmqc.tsv : peak_id  res_type  H_ppm  C_ppm  truth_label
    One 2D methyl peak per methyl group.  (1H, 13C) = data-graph node
    coordinates; residue *type* is known (from the labeling / spectral
    region); the *identity* (truth_label) is what MAUS must recover and is
    carried only for scoring.

noesy.tsv : peak_id  H1  C1  H2  C2  mix
    Methyl-methyl NOESY cross peaks.  Endpoint coordinates are the two
    methyls' (1H,13C) shifts; `mix` in {short,long} tags the mixing-time
    class.  Cross peaks are generated for structurally close methyl pairs;
    MAUS re-matches each endpoint back to hmqc.tsv by frequency.
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
  ap.add_argument('--noe-short', type=float, default=6.0, help='short-mix NOE cutoff (A)')
  ap.add_argument('--noe-long', type=float, default=8.0, help='long-mix NOE cutoff (A)')
  ap.add_argument('--keep-k', type=int, default=12, help='nearest-K NOE partners per methyl')
  ap.add_argument('--hmqc-out', default='examples/mbp/hmqc.tsv')
  ap.add_argument('--noesy-out', default='examples/mbp/noesy.tsv')
  args = ap.parse_args(argv)

  coords = parse_structure_coords(Path(args.pdb).read_text().splitlines())
  shifts = parse_bmrb_shifts(Path(args.bmrb).read_text().splitlines())
  methyls = build_methyls(coords, shifts)

  # --- HMQC peak list ---
  hmqc_lines = ['peak_id\tres_type\tH_ppm\tC_ppm\ttruth_label']
  for i, m in enumerate(methyls, 1):
    hmqc_lines.append(f'P{i}\t{m["res_type"]}\t{m["H"]:.3f}\t{m["C"]:.3f}\t{m["label"]}')
  Path(args.hmqc_out).write_text('\n'.join(hmqc_lines) + '\n')

  # --- NOESY peak list: cross peaks for close methyl pairs (nearest-K, symmetric) ---
  noe = {}  # (i,j) i<j -> 'short'/'long'
  for a in range(len(methyls)):
    dists = sorted(
      (math.dist(methyls[a]['xyz'], methyls[b]['xyz']), b)
      for b in range(len(methyls)) if b != a)
    for d, b in dists[:args.keep_k]:
      if d > args.noe_long:
        break
      key = (min(a, b), max(a, b))
      mix = 'short' if d < args.noe_short else 'long'
      # keep the tighter class if seen twice
      if key not in noe or (noe[key] == 'long' and mix == 'short'):
        noe[key] = mix

  noesy_lines = ['peak_id\tH1\tC1\tH2\tC2\tmix']
  for k, ((a, b), mix) in enumerate(sorted(noe.items()), 1):
    ma, mb = methyls[a], methyls[b]
    noesy_lines.append(
      f'X{k}\t{ma["H"]:.3f}\t{ma["C"]:.3f}\t{mb["H"]:.3f}\t{mb["C"]:.3f}\t{mix}')
  Path(args.noesy_out).write_text('\n'.join(noesy_lines) + '\n')

  # --- degeneracy report ---
  seen: Dict[Tuple[float, float], int] = {}
  for m in methyls:
    key = (round(m['H'], 2), round(m['C'], 2))
    seen[key] = seen.get(key, 0) + 1
  degenerate = sum(1 for m in methyls
                   if seen[(round(m['H'], 2), round(m['C'], 2))] > 1)
  print(f'methyls with BMRB shift + structure = {len(methyls)}')
  print(f'HMQC peaks written  = {len(methyls)}  -> {args.hmqc_out}')
  print(f'NOESY cross peaks   = {len(noe)}  (short={sum(v=="short" for v in noe.values())}'
        f' long={sum(v=="long" for v in noe.values())})  -> {args.noesy_out}')
  print(f'near-degenerate HMQC peaks (>=2 within 0.01/0.01 ppm bin) = {degenerate}')
  return 0


if __name__ == '__main__':
  raise SystemExit(main())
