# MAUS example — maltose-binding protein

Input for the SAT assigner `maus.py` (see repo root). MAUS now consumes real
**peak lists**; nothing about the answer is baked into the indexing.

| file | contents |
|---|---|
| `1ANF.pdb` | MBP crystal structure (methyl-carbon coordinates → structure graph G) |
| `hmqc.tsv` | 192 methyl HMQC **input** peaks: `label ⇥ H_ppm ⇥ C_ppm ⇥ res_type` (anonymous `P1…`) |
| `hmqc_true.tsv` | truth key: `label ⇥ H_ppm ⇥ C_ppm ⇥ res_type ⇥ True` (scoring only) |
| `noesy.tsv` | 1650 3D `(H)CCH` NOESY cross peaks: `label ⇥ C1 ⇥ C2 ⇥ H2` (825 pairs ×2 directions) |
| `hmbc.tsv` | 50 optional HMBC-HMQC geminal links (Leu/Val): `label ⇥ C1 ⇥ C2 ⇥ H2` |
| `hmqc_tentative.tsv` | same as `hmqc.tsv` but with 24 tentative anchors in the `res_type` cell |
| `mbp_options.tsv` | reference output: per-peak option sets from the run below |

Both peak lists are generated from **BMRB 7114 chemical shifts** and the 1ANF
geometry by `make_peaklists.py` (repo root). The HMQC (¹H,¹³C) coordinates are the
real methyl shifts; the NOESY endpoints are re-matched back to HMQC peaks **by
frequency** during the run, so shift degeneracy produces genuine ambiguity.

Regenerate the peak lists (optional — they are committed):

```bash
python make_peaklists.py examples/mbp/1ANF.pdb bmr7114_3.str
```

Run MAUS:

```bash
python maus.py examples/mbp/1ANF.pdb examples/mbp/hmqc.tsv examples/mbp/noesy.tsv \
    --truth examples/mbp/hmqc_true.tsv --tol-h 0.01 --tol-c 0.05 \
    --out examples/mbp/mbp_options.tsv
```

Expected (¹H ±0.01 / ¹³C ±0.05 ppm):

```
methyls(G nodes)=192  HMQC peaks=192  NOESY cross peaks=1650
NOE match: firm=502 ambiguous(dropped)=1148 unmatched=0
unique(1 option)      = 51/192
ambiguous(2-3 options)= 81/192
ambiguous(>3 options) = 60/192
unassigned            = 0/192
truth in option set   = 192/192 = 100.0%  (error rate 0.0%)
```

**What the numbers mean.** `truth in option set = 100%` is the MAUS guarantee: a
valid assignment is provably never excluded. The resolution is capped by the
**3D `(H)CCH` format**: each NOESY peak gives the observed methyl's proton but
only the *carbon* of the NOE partner, and ¹³C alone is degenerate — so 1148 of
1650 cross peaks match several candidate partners and drop as ambiguous. Carbon
tolerance is the main lever:

| `--tol-c` | firm NOE edges | unique |
|---|---|---|
| 0.10 | 214 | 29/192 |
| 0.05 | 428 | 51/192 |
| 0.02 | 630 | 70/192 |

A 4D experiment resolving the partner's proton too would recover far more; MAUS
never guesses in either case (geminal Val γ1/γ2 and Leu δ1/δ2 pairs stay
2-option, intrinsically unresolvable).

## Tentative anchors

`hmqc_tentative.tsv` is `hmqc.tsv` with 24 peaks carrying a tentative assignment
in the `res_type` cell (e.g. `L7D1`, `A21B`, `V37G1`) instead of a bare type.
Each pins that peak to one methyl; the constraint propagates through the NOE
network:

```bash
python maus.py examples/mbp/1ANF.pdb examples/mbp/hmqc_tentative.tsv examples/mbp/noesy.tsv \
    --truth examples/mbp/hmqc_true.tsv --tol-h 0.01 --tol-c 0.05
```

```
tentative anchors used = 24
unique(1 option)      = 79/192      (vs 51/192 with no anchors)
ambiguous(2-3 options)= 68/192
ambiguous(>3 options) = 45/192
truth in option set   = 192/192 = 100.0%
```

24 anchors lift unique calls 51 → 79 while the never-exclude guarantee holds.

Compare with MAGIC on the same data in [`../../COMPARISON.md`](../../COMPARISON.md)
(MAGIC lives in the sibling `../magic/` project).
