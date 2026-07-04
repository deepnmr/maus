# MAUS example — maltose-binding protein

Input for the SAT assigner `maus.py` (see repo root). MAUS now consumes real
**peak lists**; nothing about the answer is baked into the indexing.

| file | contents |
|---|---|
| `1ANF.pdb` | MBP crystal structure (methyl-carbon coordinates → structure graph G) |
| `hmqc.tsv` | 192 methyl HMQC **input** peaks: `label ⇥ H_ppm ⇥ C_ppm ⇥ res_type` (anonymous `P1…`) |
| `hmqc_true.tsv` | truth key: `label ⇥ H_ppm ⇥ C_ppm ⇥ res_type ⇥ True` (scoring only) |
| `noesy.tsv` | 825 methyl-methyl NOESY cross peaks: `peak_id ⇥ H1 ⇥ C1 ⇥ H2 ⇥ C2 ⇥ mix` |
| `hmbc.tsv` | 50 optional HMBC-HMQC geminal links (Leu/Val): `label ⇥ C1 ⇥ C2 ⇥ H` |
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
    --truth examples/mbp/hmqc_true.tsv --tol-h 0.02 --tol-c 0.2 \
    --out examples/mbp/mbp_options.tsv
```

Expected (default tolerance ¹H ±0.02 / ¹³C ±0.2 ppm):

```
methyls(G nodes)=192  HMQC peaks=192  NOESY cross peaks=825
NOE match: firm=508 ambiguous(dropped)=317 unmatched=0
unique(1 option)      = 130/192
ambiguous(2-3 options)= 19/192
ambiguous(>3 options) = 43/192
unassigned            = 0/192
truth in option set   = 192/192 = 100.0%  (error rate 0.0%)
```

**What the numbers mean.** `truth in option set = 100%` is the MAUS guarantee: a
valid assignment is provably never excluded. The residual ambiguity is real —
NOESY cross peaks whose endpoints match more than one HMQC peak (317 of 825 at
this tolerance) cannot be assigned to a definite methyl pair and are dropped, so
methyls with degenerate shifts keep several options (e.g. the geminal Val γ1/γ2
and Leu δ1/δ2 pairs, and shift-degenerate Ala/Thr/Leu clusters). Tightening
`--tol-h/--tol-c` recovers more unique calls (170/192 at ±0.01/0.1) as fewer NOEs
are ambiguous — the honest resolution/degeneracy trade-off.

Compare with MAGIC on the same data in [`../../COMPARISON.md`](../../COMPARISON.md)
(MAGIC lives in the sibling `../magic/` project).
