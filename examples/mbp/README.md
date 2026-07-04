# MAUS example — maltose-binding protein

Input for the standalone SAT assigner `maus.py` (see repo root).

| file | contents |
|---|---|
| `1ANF.pdb` | MBP crystal structure (methyl-carbon coordinates → structure graph G) |
| `mbp_peaks.tsv` | 192 peaks: `peak_id ⇥ residue_type ⇥ truth_label` (truth from BMRB 7114) |
| `mbp_options.tsv` | reference output: per-peak option sets from the run below |

Run:

```bash
python maus.py examples/mbp/1ANF.pdb examples/mbp/mbp_peaks.tsv \
    --keep-k 8 --out mbp_options.tsv
```

Expected:

```
unique(1 option)      = 176/192
ambiguous(2-3 options)= 16/192
unassigned            = 0/192
truth in option set   = 192/192 = 100.0%  (error rate 0.0%)
```

The 16 ambiguous peaks are genuine geminal/local symmetries (Val γ1/γ2, Leu δ1/δ2,
two Ala pairs, two Thr pairs) that an achiral NOE network cannot resolve without
stereospecific labeling. MAUS reports both members rather than guessing — hence
0 % error. Compare with MAGIC on the same data in [`../../COMPARISON.md`](../../COMPARISON.md)
(MAGIC lives in the sibling `../magic/` project).
