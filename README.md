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
- **Data graph H** — nodes are **HMQC peaks** (`peak_id`, residue type, ¹H/¹³C
  shift); edges are **NOESY cross peaks**, each endpoint matched back to an HMQC
  peak *by frequency* within tolerance and tagged short/long by mixing time.
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
python maus.py PDB HMQC.tsv NOESY.tsv [--truth TRUTH.tsv] [options]
```

- `HMQC.tsv` (input) — `label ⇥ H_ppm ⇥ C_ppm ⇥ res_type`. `label` is an
  anonymous peak id (`P1, P2, …`) that leaks nothing about the answer.
- `NOESY.tsv` — `peak_id ⇥ H1 ⇥ C1 ⇥ H2 ⇥ C2 ⇥ mix` (`mix` ∈ short/long)
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
`True` column), and `noesy.tsv` (methyl-methyl cross peaks for structurally close
pairs, endpoint coordinates = the two methyls' shifts).

## Example — maltose-binding protein

```bash
python maus.py examples/mbp/1ANF.pdb examples/mbp/hmqc.tsv examples/mbp/noesy.tsv \
    --truth examples/mbp/hmqc_true.tsv --tol-h 0.02 --tol-c 0.2 \
    --out examples/mbp/mbp_options.tsv
```

```
methyls(G nodes)=192  HMQC peaks=192  NOESY cross peaks=825
NOE match: firm=508 ambiguous(dropped)=317 unmatched=0
unique(1 option)      = 130/192
ambiguous(2-3 options)= 19/192
ambiguous(>3 options) = 43/192
unassigned            = 0/192
truth in option set   = 192/192 = 100.0%  (error rate 0.0%)
```

Real BMRB 7114 shifts + 1ANF geometry, fed as HMQC + NOESY peak lists. `truth in
option set = 100%` is the MAUS guarantee (a valid assignment is never excluded).
The residual ambiguity is real: 317 of 825 NOESY cross peaks have a degenerate
endpoint that matches more than one HMQC peak, so they cannot be pinned to a
methyl pair and methyls with degenerate shifts keep several options. Tightening
`--tol-h/--tol-c` recovers more unique calls (170/192 at ±0.01/0.1) — the honest
resolution/degeneracy trade-off. See [`examples/mbp/`](examples/mbp) and
[`COMPARISON.md`](COMPARISON.md).
