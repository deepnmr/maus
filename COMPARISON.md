# MAUS vs. MAGIC on maltose-binding protein

This note documents a head-to-head run of two different methyl-assignment
algorithms on the **same** maltose-binding protein (MBP) dataset:

- **MAGIC** (Monneau et al., *J Biomol NMR* 2017) — the scoring / exhaustive-search
  method implemented in the sibling `../magic/` project. Produces a single best
  assignment per peak from a confidence-weighted objective.
- **MAUS** (Nerli, De Paula, McShan & Sgourakis, *Nat. Commun.* 2021, 12:691) —
  the SAT / subgraph-isomorphism method implemented in `maus.py` (clean-room,
  separate program, this directory). Produces, for each peak, the *set* of methyls
  that are consistent with all hard constraints — never a single guess unless the
  data force one.

## Shared input

| item | value |
|---|---|
| protein | *E. coli* maltose-binding protein |
| chemical shifts | **real** experimental methyl ¹³C/¹H from **BMRB entry 7114** |
| structure | PDB **1ANF** |
| methyls / peaks | 192 (A 44, I 22, L 60, M 6, T 20, V 40) |
| NOESY | *simulated* from 1ANF geometry (BMRB has no NOESY peak list); nearest-8 within 8 Å |
| ground truth | BMRB residue label of every peak |

Identical shifts, identical structure, identical simulated NOE network feed both
programs. The only difference is the algorithm.

## Result

| metric | **MAGIC** (scoring) | **MAUS** (SAT) |
|---|---|---|
| assignment style | one methyl per peak | option set per peak |
| **error rate** (truth excluded) | **84.9 %** | **0.0 %** |
| truth present in output | 31/192 = 16.1 % | **192/192 = 100 %** |
| peaks uniquely assigned | 192/192 (commits to all) | 176/192 = 91.7 % |
| of those, correct | 29/192 = 15.1 % | 176/176 = 100 % |
| peaks left as 2–3 options | – | 16/192 |
| peaks unassigned | 0 | 0 |
| wall-clock | 279 s | ~seconds |

Per-type MAGIC recovery: A 9 %, I 23 %, L 10 %, M 67 %, T 25 %, V 12 %.

## Why MAUS wins decisively on this data

The gap is not a bug in MAGIC — it is a structural consequence of how each
method uses the NOESY, amplified by the fact that the NOESY here is *simulated
from the structure*:

1. **MAUS uses the NOE as a hard constraint.** Every symmetric NOE edge must map
   onto a structure-graph edge of a compatible distance class. When the NOE
   network is derived from the structure, each NOE edge corresponds to a real
   contact, so the data graph is (almost) a subgraph of the structure graph and
   the isomorphism is nearly deterministic. MAUS therefore prunes to the true
   methyl for 92 % of peaks and, crucially, **never drops the truth** — the
   remaining 8 % collapse to genuine 2–3-way degeneracies.

2. **MAGIC uses the NOE as a soft score.** Its objective is a near-flat landscape
   over a structure-simulated NOESY (documented in `../magic/VALIDATION.md` §2):
   the corrected Eq. (1) objective does rank the *true* global assignment highest,
   but the runner-ups sit ~1–2 % below, so the bounded search commits to a
   near-optimal-but-wrong single answer for most peaks. Committing to one guess
   with no "I'm not sure" option is what produces the 85 % error rate.

3. **The 16 MAUS ambiguities are real symmetries, not failures.** They are almost
   entirely geminal methyl pairs of the same residue (V8 γ1/γ2, L43 δ1/δ2, L75,
   L121, V240, L311, …) plus two locally-indistinguishable Ala pairs (A84/A141)
   and two Thr pairs (T237/T356). Without stereospecific labeling these swaps are
   genuinely unresolvable from an achiral NOE network — MAUS correctly reports
   both members rather than guessing.

## Takeaway

On MBP with structure-consistent NOE information, the **constraint-satisfaction**
formulation (MAUS) is the right tool: it returns an assignment set that provably
contains the truth for every peak and pins 92 % of them uniquely, in seconds. The
**scoring** formulation (MAGIC) is designed for *experimental* NOESY, where
intensities and crosspeak sparsity carry information a boolean edge cannot; on a
simulated network its single-answer commitment is a liability, not a feature.

The two are complementary, not redundant — MAUS to bound the space with certainty,
MAGIC-style scoring to rank within the residual degeneracies when real NOESY
intensities are available.

## Reproduce

```bash
# MAUS (this project, separate program)
python maus.py examples/mbp/1ANF.pdb examples/mbp/mbp_peaks.tsv \
    --keep-k 8 --out mbp_options.tsv

# MAGIC (same shifts + structure + simulated NOE) — see ../magic/VALIDATION.md
```

`examples/mbp/mbp_options.tsv` holds the per-peak MAUS option sets used above.
