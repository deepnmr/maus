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

- **Structure graph G** — one node per methyl carbon (from the PDB); edges
  classified geminal / short-range (`< short_cut`) / long-range (`< long_cut`).
- **Data graph H** — 2D reference peaks (residue type known) with symmetric NOE
  edges, split into short- and long-mixing classes.
- **SAT encoding** — variable `x(i,j)` = peak *i* → methyl *j*. Hard constraints:
  exactly-one methyl per peak, injective map, geminal edges preserved, and every
  NOE edge mapped onto a structure edge of a compatible distance class.
- **Per-peak options** — a methyl is a valid option iff the CNF stays satisfiable
  when `x(i,j)` is asserted; enumerated incrementally with the solver's assumption
  interface (the paper's iterative ansatz).

## Install

```bash
python3 -m pip install python-sat
```

## Usage

```bash
python maus.py PDB PEAKS.tsv [options]
```

`PEAKS.tsv` is tab-separated: `peak_id ⇥ residue_type ⇥ truth_label`.

| option | default | meaning |
|---|---|---|
| `--short-cut` | 6.0 | structure short-range edge cutoff (Å) |
| `--long-cut` | 10.0 | structure long-range edge cutoff (Å) |
| `--noe-short` | 6.0 | simulated short-mixing NOE cutoff (Å) |
| `--noe-long` | 8.0 | simulated long-mixing NOE cutoff (Å) |
| `--keep-k` | 8 | nearest-K NOE partners kept per methyl (data sparsity) |
| `--labeling` | `A;I;L;M;T;V` | residue types present |
| `--out` | – | write per-peak options TSV |

## Example — maltose-binding protein

```bash
python maus.py examples/mbp/1ANF.pdb examples/mbp/mbp_peaks.tsv \
    --keep-k 8 --out mbp_options.tsv
```

```
methyls(G nodes)=192  peaks=192
unique(1 option)      = 176/192
ambiguous(2-3 options)= 16/192
unassigned            = 0/192
truth in option set   = 192/192 = 100.0%  (error rate 0.0%)
```

Real BMRB 7114 shifts + 1ANF geometry. 176 peaks pinned uniquely; the 16 residual
2–3-option peaks are genuine geminal/local symmetries unresolvable without
stereospecific labeling. See [`examples/mbp/`](examples/mbp) and
[`COMPARISON.md`](COMPARISON.md).
