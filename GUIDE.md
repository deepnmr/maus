# MAUS User Guide

> **Note:** This is a placeholder ("dummy") guide scaffold. Replace the
> bracketed `[TODO]` sections with project-specific content before publishing.

*Methyl Assignments Using Satisfiability* — a clean-room reimplementation of the
SAT-based methyl NMR assignment method of Nerli et al., *Nat. Commun.* 12:691
(2021).

---

## 1. Overview

MAUS assigns methyl NMR resonances by **subgraph isomorphism** rather than by
scoring. Given a protein structure and a set of 2D reference peaks with NOE
connectivity, it returns — for each peak — the *set* of methyls consistent with
every hard constraint. It never emits a single guess unless the data force one,
and it provably never excludes a valid assignment.

| Concept | In this tool |
|---|---|
| Structure graph **G** | one node per methyl carbon; geminal / short / long edges |
| Data graph **H** | 2D peaks (residue type known) + symmetric NOE edges |
| Solver | `python-sat` (Glucose3), assumption-based enumeration |
| Output | per-peak option set (1 / 2–3 / >3 / unassigned) |

---

## 2. Installation

```bash
python3 -m pip install python-sat
git clone https://github.com/deepnmr/maus.git
cd maus
```

Requirements:

- Python ≥ 3.8
- `python-sat` (provides the Glucose3 SAT solver)

---

## 3. Quick start

```bash
python maus.py examples/mbp/1ANF.pdb examples/mbp/mbp_peaks.tsv \
    --keep-k 8 --out mbp_options.tsv
```

Expected summary:

```
methyls(G nodes)=192  peaks=192
unique(1 option)      = 176/192
ambiguous(2-3 options)= 16/192
unassigned            = 0/192
truth in option set   = 192/192 = 100.0%
```

---

## 4. Input formats

### 4.1 PDB structure

Standard `ATOM` records. Only these residue types are parsed by default:
`ALA ILE LEU MET THR VAL`. Chain `A` (or blank) is used.

### 4.2 Peaks TSV

Tab-separated, one peak per line:

```
peak_id <TAB> residue_type <TAB> truth_label
P1        L                 L7CD1
P2        L                 L7CD2
```

- `peak_id` — arbitrary unique string
- `residue_type` — one-letter code (`A I L M T V`)
- `truth_label` — ground-truth methyl label (for benchmarking only)

---

## 5. Command-line options

| Option | Default | Meaning |
|---|---|---|
| `--short-cut` | 6.0 | structure short-range edge cutoff (Å) |
| `--long-cut` | 10.0 | structure long-range edge cutoff (Å) |
| `--noe-short` | 6.0 | simulated short-mixing NOE cutoff (Å) |
| `--noe-long` | 8.0 | simulated long-mixing NOE cutoff (Å) |
| `--keep-k` | 8 | nearest-K NOE partners kept per methyl |
| `--labeling` | `A;I;L;M;T;V` | residue types present |
| `--out` | – | write per-peak options TSV |

---

## 6. Interpreting output

The `--out` TSV columns:

| Column | Meaning |
|---|---|
| `peak_id` | input peak identifier |
| `res_type` | one-letter residue type |
| `n_options` | number of valid methyls |
| `options` | comma-separated methyl labels |
| `truth` | ground-truth label |
| `truth_in_set` | 1 if truth ∈ options, else 0 |

- **1 option** → uniquely assigned.
- **2–3 options** → residual ambiguity, usually geminal/local symmetry.
- **>3 options** → weak local NOE network.
- **0 options** → over-constrained inputs (check cutoffs).

---

## 7. Caveat on the bundled benchmark

> **Important.** In the shipped MBP example the NOE data graph **H** is
> *simulated from the same coordinates* used to build **G**, and peaks are
> indexed identically to methyls. The `truth in option set = 100%` figure is
> therefore a self-consistency check on ideal, noise-free, complete data — not a
> measurement of real-data assignment performance. To reproduce the paper's
> validation, supply an **independent** experimental NOE peak list.

---

## 8. Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `ModuleNotFoundError: pysat` | solver not installed | `pip install python-sat` |
| many `unassigned` peaks | cutoffs too tight | raise `--long-cut` / `--noe-long` |
| everything ambiguous | NOE network too sparse | raise `--keep-k` |

---

## 9. References

- Nerli, De Paula, McShan & Sgourakis. *Backbone-independent NMR resonance
  assignments of methyl probes in large proteins.* **Nat. Commun.** 12:691
  (2021). [doi:10.1038/s41467-021-20984-0](https://doi.org/10.1038/s41467-021-20984-0)
- See [`COMPARISON.md`](COMPARISON.md) for the head-to-head against MAGIC.

---

*[TODO] Add author, license, and contact sections before release.*
