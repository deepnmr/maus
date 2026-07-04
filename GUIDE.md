# MAUS User Guide

> **Note:** This is a placeholder ("dummy") guide scaffold. Replace the
> bracketed `[TODO]` sections with project-specific content before publishing.

*Methyl Assignments Using Satisfiability* — a clean-room reimplementation of the
SAT-based methyl NMR assignment method of Nerli et al., *Nat. Commun.* 12:691
(2021).

---

## 1. Overview

MAUS assigns methyl NMR resonances by **subgraph isomorphism** rather than by
scoring. Given a protein structure and real **peak lists** (HMQC + NOESY), it
returns — for each peak — the *set* of methyls consistent with every hard
constraint. It never emits a single guess unless the data force one, and it
provably never excludes a valid assignment.

| Concept | In this tool |
|---|---|
| Structure graph **G** | one node per methyl carbon; geminal / short / long edges |
| Data graph **H** nodes | HMQC peaks: (¹H, ¹³C) shift + residue type |
| Data graph **H** edges | NOESY cross peaks, endpoints matched to HMQC by frequency |
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
# 1. build peak lists from a BMRB shift file + the PDB
python make_peaklists.py examples/mbp/1ANF.pdb examples/mbp/bmr7114_3.str

# 2. run MAUS on the HMQC + NOESY peak lists
python maus.py examples/mbp/1ANF.pdb examples/mbp/hmqc.tsv examples/mbp/noesy.tsv \
    --tol-h 0.02 --tol-c 0.2 --out mbp_options.tsv
```

Expected summary:

```
methyls(G nodes)=192  HMQC peaks=192  NOESY cross peaks=825
NOE match: firm=508 ambiguous(dropped)=317 unmatched=0
unique(1 option)      = 130/192
ambiguous(2-3 options)= 19/192
ambiguous(>3 options) = 43/192
unassigned            = 0/192
truth in option set   = 192/192 = 100.0%
```

---

## 4. Input formats

### 4.1 PDB structure

Standard `ATOM` records. Only these residue types are parsed by default:
`ALA ILE LEU MET THR VAL`. Chain `A` (or blank) is used.

### 4.2 HMQC peak list (TSV) — data-graph nodes

```
peak_id <TAB> res_type <TAB> H_ppm <TAB> C_ppm <TAB> truth_label
P1        L                 0.828         24.510      L7CD1
```

- `peak_id` — arbitrary unique string
- `res_type` — one-letter code (`A I L M T V`), known from labeling
- `H_ppm`, `C_ppm` — methyl (¹H, ¹³C) chemical shifts
- `truth_label` — ground-truth methyl label (for scoring only)

### 4.3 NOESY peak list (TSV) — data-graph edges

```
peak_id <TAB> H1 <TAB> C1 <TAB> H2 <TAB> C2 <TAB> mix
X1        0.828    24.51   0.712    23.10   short
```

Each row is a methyl-methyl cross peak. During the run both endpoints are
matched back to HMQC peaks by frequency (within `--tol-h`/`--tol-c`); `mix` ∈
`{short, long}` tags the mixing-time distance class.

### 4.4 Generating the peak lists

`make_peaklists.py` builds both lists from a PDB and a BMRB NMR-STAR shift file
(`bmrXXXX_3.str`): the HMQC (¹H,¹³C) coordinates are the real BMRB methyl shifts,
and NOESY cross peaks are emitted for structurally close methyl pairs.

---

## 5. Command-line options

| Option | Default | Meaning |
|---|---|---|
| `--short-cut` | 6.0 | structure short-range edge cutoff (Å) |
| `--long-cut` | 10.0 | structure long-range edge cutoff (Å) |
| `--tol-h` | 0.02 | ¹H NOESY→HMQC match tolerance (ppm) |
| `--tol-c` | 0.20 | ¹³C NOESY→HMQC match tolerance (ppm) |
| `--labeling` | `A;I;L;M;T;V` | residue types present |
| `--out` | – | write per-peak options TSV |

`make_peaklists.py` options: `--noe-short` (6.0), `--noe-long` (8.0), `--keep-k`
(12 nearest NOE partners per methyl).

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

## 7. On the bundled benchmark

MAUS sees only the two peak lists. HMQC (¹H,¹³C) coordinates are the real BMRB
7114 shifts, and each NOESY endpoint is resolved back to an HMQC peak **by
frequency**, not by identity — so the assignment is not circular and shift
degeneracy has real consequences.

- `truth in option set = 100%` is the MAUS **guarantee**: a valid assignment is
  provably never excluded (not a measurement — a property of the exact
  enumeration).
- The residual ambiguity is **real**. NOESY cross peaks whose endpoint matches
  more than one HMQC peak cannot be pinned to a definite methyl pair and are
  dropped, so shift-degenerate methyls keep several options. Tighter
  `--tol-h/--tol-c` recovers more unique calls (170/192 at ±0.01/0.1) as fewer
  NOEs are ambiguous — the honest resolution/degeneracy trade-off.

The NOESY cross peaks themselves are still *generated from the structure* (close
methyl pairs), so this remains a controlled benchmark rather than a picked
experimental NOESY. Swap in a real NOESY peak list to go fully experimental — the
`maus.py` interface does not change.

---

## 8. Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `ModuleNotFoundError: pysat` | solver not installed | `pip install python-sat` |
| many `unassigned` peaks | structure cutoffs too tight | raise `--long-cut` |
| most peaks ambiguous | too many NOEs dropped as ambiguous | tighten `--tol-h`/`--tol-c` |
| high `unmatched` NOE count | tolerance too tight vs shift precision | loosen `--tol-h`/`--tol-c` |

---

## 9. References

- Nerli, De Paula, McShan & Sgourakis. *Backbone-independent NMR resonance
  assignments of methyl probes in large proteins.* **Nat. Commun.** 12:691
  (2021). [doi:10.1038/s41467-021-20984-0](https://doi.org/10.1038/s41467-021-20984-0)
- See [`COMPARISON.md`](COMPARISON.md) for the head-to-head against MAGIC.

---

*[TODO] Add author, license, and contact sections before release.*
