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

# 2. run MAUS on the HMQC + NOESY peak lists (--truth scores against the key)
python maus.py examples/mbp/1ANF.pdb examples/mbp/hmqc.tsv examples/mbp/noesy.tsv \
    --truth examples/mbp/hmqc_true.tsv --tol-h 0.01 --tol-c 0.05 --out mbp_options.tsv
```

Expected summary:

```
methyls(G nodes)=192  HMQC peaks=192  NOESY cross peaks=1650
NOE match: firm=502 ambiguous(dropped)=1148 unmatched=0
unique(1 option)      = 51/192
ambiguous(2-3 options)= 81/192
ambiguous(>3 options) = 60/192
unassigned            = 0/192
truth in option set   = 192/192 = 100.0%
```

---

## 4. Input formats

### 4.1 PDB structure

Standard `ATOM` records. Only these residue types are parsed by default:
`ALA ILE LEU MET THR VAL`. Chain `A` (or blank) is used.

### 4.2 HMQC peak list (TSV) — data-graph nodes (input)

```
label <TAB> H_ppm <TAB> C_ppm <TAB> res_type
P1      0.828        24.510       L
P7      0.340        25.390       L45D2
```

- `label` — anonymous peak id (`P1, P2, …`); leaks nothing about the answer
- `H_ppm`, `C_ppm` — methyl (¹H, ¹³C) chemical shifts
- `res_type` — one-letter code (`A I L M T V`), known from labeling, **or** a
  *tentative assignment* like `L45D2` (Leu45 Cδ2; the leading `C` is optional,
  so `L45D2` = `L45CD2`). A tentative cell pins that peak's candidate methyls to
  the single named one; the constraint then propagates through the NOE network
  to reduce ambiguity elsewhere. A tentative label absent from the structure is
  ignored (falls back to the residue type). Runs print `tentative anchors used`.

Ground truth lives in a **separate** file, never in the input:

```
label <TAB> H_ppm <TAB> C_ppm <TAB> res_type <TAB> True
P1      0.828        24.510       L             L7CD1
```

Pass it with `--truth` to score; omit it to run blind.

### 4.3 NOESY peak list (TSV) — data-graph edges

```
label <TAB> C1 <TAB> C2 <TAB> H2
X1      24.56    25.39   0.340
```

3D `(H)CCH` methyl-methyl cross peak. The **observed** methyl is `(H2, C2)`; the
NOE **partner** contributes carbon `C1` only (a 3D peak has no partner proton).
During the run the observed methyl is matched to an HMQC peak by `(H2,C2)` and
the partner by carbon `C1` alone (within `--tol-h`/`--tol-c`). Single distance
class (structure edge ≤ `long_cut`). Because ¹³C alone is degenerate, most cross
peaks match several candidate partners and drop as ambiguous — carbon tolerance
is the main resolution lever.

### 4.3b HMBC-HMQC peak list (TSV, optional) — geminal links

```
label <TAB> C1 <TAB> C2 <TAB> H2
B1      24.56    25.39   0.340
```

Same layout as the NOESY list. One row per Leu/Val residue: observed methyl
`(H2, C2)`, geminal partner carbon `C1`. Pass with `--hmbc`; MAUS forces the pair
onto a geminal structure edge. The carbon-only partner is degenerate (on MBP
~1/50 links resolve), and even a firm link couples the pair without fixing which
residue — so per-peak option *counts* are essentially unchanged. The
never-exclude guarantee is preserved.

### 4.4 Generating the peak lists

`make_peaklists.py` builds all three files from a PDB and a BMRB NMR-STAR shift
file (`bmrXXXX_3.str`): `hmqc.tsv` (input), `hmqc_true.tsv` (truth key), and
`noesy.tsv`. HMQC (¹H,¹³C) coordinates are the real BMRB methyl shifts; NOESY
cross peaks are emitted for structurally close methyl pairs.

---

## 5. Command-line options

| Option | Default | Meaning |
|---|---|---|
| `--short-cut` | 6.0 | structure short-range edge cutoff (Å) |
| `--long-cut` | 10.0 | structure long-range edge cutoff (Å) |
| `--tol-h` | 0.02 | ¹H NOESY/HMBC→HMQC match tolerance (ppm) |
| `--tol-c` | 0.20 | ¹³C NOESY/HMBC→HMQC match tolerance (ppm) |
| `--hmbc` | – | optional HMBC-HMQC geminal-link peak list |
| `--truth` | – | truth-key TSV for scoring |
| `--labeling` | `A;I;L;M;T;V` | residue types present |
| `--out` | – | write per-peak options TSV |

`make_peaklists.py` options: `--noe-short` (6.0), `--noe-long` (8.0), `--keep-k`
(12 nearest NOE partners per methyl).

---

## 6. Interpreting output

The `--out` TSV columns:

| Column | Meaning |
|---|---|
| `label` | input peak identifier (anonymous) |
| `res_type` | one-letter residue type |
| `n_options` | number of valid methyls |
| `options` | comma-separated methyl labels |
| `truth` | ground-truth label (blank if no `--truth`) |
| `truth_in_set` | 1 if truth ∈ options, else 0 (blank if no `--truth`) |

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
- The residual ambiguity is **real** and dominated by the 3D format: the NOE
  partner is matched on ¹³C only, which is degenerate, so most cross peaks match
  several candidate partners and drop. Carbon tolerance is the main lever
  (unique 29 → 51 → 70 /192 at `--tol-c` 0.1 / 0.05 / 0.02). A 4D experiment
  resolving the partner proton would recover far more.

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
