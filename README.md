# maus

Clean-room implementation of **MAUS** — *Methyl Assignments Using Satisfiability*
(Nerli, De Paula, McShan & Sgourakis, *Nat. Commun.* 12:691, 2021,
[doi:10.1038/s41467-021-20984-0](https://doi.org/10.1038/s41467-021-20984-0)).

MAUS assigns methyl NMR resonances by **subgraph isomorphism** solved with a SAT
solver, rather than by scoring. It returns, for each 2D peak, the *set* of methyls
consistent with every hard constraint — never a single guess unless the data force
one, and it provably never excludes a valid assignment.

Companion to the sibling [`../magic/`](../magic) project (the scoring-based MAGIC
method). See [`COMPARISON.md`](COMPARISON.md) for a head-to-head on maltose-binding
protein.

## How it works

MAUS reads real **peak lists** — the same data an experiment produces. Nothing
about the answer is baked into the indexing.

- **Structure graph G** — one node per methyl carbon (from the PDB); edges
  classified geminal / short-range (`< short_cut`) / long-range (`< long_cut`).
- **Data graph H** — nodes are **HMQC peaks** (`label`, residue type, ¹H/¹³C
  shift); edges are **NOESY cross peaks** (3D `(H)CCH`: `label C1 C2 H2`). The
  observed methyl is matched back to an HMQC peak by `(H2,C2)`; the NOE partner
  by carbon `C1` only (its proton is absent from a 3D peak).
- **SAT encoding** — variable `x(i,j)` = HMQC peak *i* → methyl *j*. Hard
  constraints: exactly-one methyl per peak, injective map, and every NOESY cross
  peak that resolves to a definite pair of HMQC peaks mapped onto a structure
  edge of a compatible distance class.
- **Per-peak options** — a methyl is a valid option iff the CNF stays satisfiable
  when `x(i,j)` is asserted; enumerated incrementally with the solver's assumption
  interface (the paper's iterative ansatz).

Because NOESY endpoints are resolved through chemical shift (not identity), shift
**degeneracy** produces genuine, irreducible ambiguity — as in the real method.

## Install

```bash
python3 -m pip install python-sat
```

## Usage

```bash
python maus.py PDB HMQC.tsv NOESY.tsv [--hmbc HMBC.tsv] [--truth TRUTH.tsv] [options]
```

- `HMQC.tsv` (input) — `label ⇥ H_ppm ⇥ C_ppm ⇥ res_type`. `label` is an
  anonymous peak id (`P1, P2, …`) that leaks nothing about the answer. The
  `res_type` cell may be a bare type (`L`) **or a tentative assignment**
  (`L45D2` = Leu45 Cδ2; the `C` may be omitted). A tentative cell pins that
  peak's domain to the named methyl, and the constraint propagates through the
  NOE network to sharpen the rest — a few anchors go a long way (on MBP 24
  anchors lift unique 51 → 79). An out-of-structure tentative label is ignored
  (falls back to the residue type).
- `NOESY.tsv` — `label ⇥ C1 ⇥ C2 ⇥ H2` (3D `(H)CCH`). Observed methyl = `(H2,C2)`;
  NOE partner = carbon `C1` only. Single distance class (≤ `long_cut`).
- `HMBC.tsv` (optional) — same layout `label ⇥ C1 ⇥ C2 ⇥ H2`. Each row links one
  Leu/Val residue's two prochiral methyls (observed `(H2,C2)`, geminal partner
  carbon `C1`); MAUS forces the pair onto a geminal structure edge. Pass with
  `--hmbc`.
- `TRUTH.tsv` (optional, scoring only) — `label ⇥ H_ppm ⇥ C_ppm ⇥ res_type ⇥ True`,
  where `True` is the real methyl for each label. Pass with `--truth` to score.

Build all three from a BMRB shift file plus a PDB with `make_peaklists.py`
(see below).

| option | default | meaning |
|---|---|---|
| `--short-cut` | 6.0 | structure short-range edge cutoff (Å) |
| `--long-cut` | 10.0 | structure long-range edge cutoff (Å) |
| `--tol-h` | 0.02 | ¹H NOESY→HMQC match tolerance (ppm) |
| `--tol-c` | 0.20 | ¹³C NOESY→HMQC match tolerance (ppm) |
| `--labeling` | `A;I;L;M;T;V` | residue types present |
| `--out` | – | write per-peak options TSV |

## Building peak lists from BMRB

```bash
python make_peaklists.py PDB bmrXXXX_3.str \
    --noe-short 6.0 --noe-long 8.0 --keep-k 12
```

Writes `hmqc.tsv` (input, anonymous labels), `hmqc_true.tsv` (truth key with the
`True` column), `noesy.tsv` (methyl-methyl cross peaks for structurally close
pairs), and `hmbc.tsv` (one geminal link per Leu/Val residue).

> **Note on HMBC.** The geminal-partner endpoint is matched on **carbon only**
> (its proton is not in the row), which is highly degenerate — on MBP only ~1 of
> 50 links resolves uniquely, the rest drop as ambiguous. Even a firm link
> couples a Leu/Val pair (which two peaks share a residue) but not *which*
> residue, so it does not change the per-peak option *counts* on MBP — the
> residual ambiguity is prochiral (CD1↔CD2 / CG1↔CG2, intrinsically
> unresolvable) or shift-degenerate. The constraint is applied and the
> never-exclude guarantee is preserved.

## Example — maltose-binding protein

```bash
python maus.py examples/mbp/1ANF.pdb examples/mbp/hmqc.tsv examples/mbp/noesy.tsv \
    --truth examples/mbp/hmqc_true.tsv --tol-h 0.01 --tol-c 0.05 \
    --out examples/mbp/mbp_options.tsv
```

```
methyls(G nodes)=192  HMQC peaks=192  NOESY cross peaks=1650
NOE match (tol H±0.01/C±0.05): firm=502 ambiguous(dropped)=1148 unmatched=0
unique(1 option)      = 51/192
ambiguous(2-3 options)= 81/192
ambiguous(>3 options) = 60/192
unassigned            = 0/192
truth in option set   = 192/192 = 100.0%  (error rate 0.0%)
```

Real BMRB 7114 shifts + 1ANF geometry, fed as HMQC + 3D `(H)CCH` NOESY peak
lists. `truth in option set = 100%` is the MAUS guarantee (a valid assignment is
never excluded).

The resolution is limited by the **3D data format**: a NOESY peak carries the
observed methyl's proton but only the *carbon* of the NOE partner, and ¹³C alone
is highly degenerate — so most cross peaks (1148 of 1650 even at ±0.05 ppm) match
several candidate partners and drop as ambiguous. Carbon tolerance is the main
lever (unique: 29 → 51 → 70 at ±0.1 / ±0.05 / ±0.02 ppm). A 4D experiment that
also resolves the partner's proton would recover far more; MAUS never guesses in
either case. See [`examples/mbp/`](examples/mbp) and
[`COMPARISON.md`](COMPARISON.md).
