# TUSCO Paper Update Plan with Line-by-Line Edits

## Data Sources
- **TUSCO Statistics**: 
  - `/data/processed/tusco/hsa/tusco_human_tissue_statistics.tsv`
  - `/data/processed/tusco/mmu/tusco_mouse_tissue_statistics.tsv`
- **Figure Tables**: 
  - `/figs/figure-03/tables/*.tsv`
  - `/figs/figure-04/tables/*.tsv`
  - `/figs/figure-05/tables/*.tsv`
  - `/figs/supp-fig-*/tables/*.tsv`

---

## Figure 1 Updates

### Update tissue-specific TUSCO gene counts
<source>: 
- Tissue ranges: `/data/processed/tusco/hsa/tusco_human_tissue_statistics.tsv` and `/data/processed/tusco/mmu/tusco_mouse_tissue_statistics.tsv`
- Universal gene counts: `/data/processed/tusco/hsa/tusco_human.tsv` and `/data/processed/tusco/mmu/tusco_mouse.tsv`

<old>:
<<<
This process yielded tissue-resolved TUSCO gene sets (human: 60–187 per tissue; mouse: 28–67) and universal cores (human n = 50; mouse n = 32)
>>>

<new>:
<<<
This process yielded tissue-resolved TUSCO gene sets (human: [MIN_FROM_STATS]–[MAX_FROM_STATS] per tissue; mouse: [MIN_FROM_STATS]–[MAX_FROM_STATS]) and universal cores (human n = [UNIVERSAL_COUNT]; mouse n = [UNIVERSAL_COUNT])
>>>

*Note: Extract MIN/MAX from tissue statistics files. Count universal genes from tusco_human.tsv and tusco_mouse.tsv files.*

---

## Figure 3 Updates

### Update cosine similarity values
<source>: `/figs/figure-03/tables/table_s1.csv`

<old>:
<<<
The cDNA PacBio libraries, which were the best-performing experimental option according to LRGASP results, consistently achieved cosim values near unity in both human and mouse samples
>>>

<new>:
<<<
The cDNA PacBio libraries, which were the best-performing experimental option according to LRGASP results, consistently achieved cosim values [EXTRACT_FROM_TABLE_S1] in both human and mouse samples
>>>

### Update cosine similarity range
<source>: `/figs/figure-03/tables/table_s1.csv`

<old>:
<<<
cDNA ONT, dRNA ONT, and hybrid long–short read pipelines (denoted "LS") also exhibited high agreement (cosim: 0.9553–0.9977)
>>>

<new>:
<<<
cDNA ONT, dRNA ONT, and hybrid long–short read pipelines (denoted "LS") also exhibited high agreement (cosim: [MIN_COSIM]–[MAX_COSIM])
>>>

### Update statistical test p-values
<source>: `/figs/figure-03/tables/figure3b-human.tsv` and `/figs/figure-03/tables/figure3b-mouse.tsv`

<old>:
<<<
PTP and FN were significantly higher under TUSCO (one-sided paired t-test: human PTP p = 0.0120, FN p = 2.21 × 10−3; mouse PTP p = 0.0232, FN p = 5.97 × 10−4)
>>>

<new>:
<<<
PTP and FN were significantly higher under TUSCO (one-sided paired t-test: human PTP p = [HUMAN_PTP_P], FN p = [HUMAN_FN_P]; mouse PTP p = [MOUSE_PTP_P], FN p = [MOUSE_FN_P])
>>>

### Update RIN correlation statistics
<source>: `/figs/figure-03/tables/figure3c.tsv`

<old>:
<<<
the fraction of fully recovered TUSCO transcripts—calculated as TP / (TP + PTP)—strongly correlated with the RIN value (R = 0.881, p = 1.41 × 10−6), whereas Sequins showed no correlation (R = 0.075, p = 0.77)
>>>

<new>:
<<<
the fraction of fully recovered TUSCO transcripts—calculated as TP / (TP + PTP)—strongly correlated with the RIN value (R = [TUSCO_R], p = [TUSCO_P]), whereas Sequins showed no correlation (R = [SEQUINS_R], p = [SEQUINS_P])
>>>

### Update sequencing depth correlation
<source>: `/figs/supp-fig-03/tables/fig-s3.tsv`

<old>:
<<<
A pronounced negative correlation (Pearson r = −0.778, p = 2.87 × 10−3) was observed
>>>

<new>:
<<<
A pronounced negative correlation (Pearson r = [PEARSON_R], p = [PEARSON_P]) was observed
>>>

---

## Figure 4 Updates

### Update per-read coverage statistics
<source>: `/figs/figure-04/tables/fig-s5.tsv`

<old>:
<<<
In human, TP coverage averaged 94.47% ± 17.24% (n = 5,694), while PTP averaged 62.96% ± 47.03% (n = 320); in mouse, TP averaged 96.66% ± 9.44% (n = 4,549) and PTP 75.58% ± 29.56% (n = 2,555)
>>>

<new>:
<<<
In human, TP coverage averaged [HUMAN_TP_MEAN]% ± [HUMAN_TP_SD]% (n = [HUMAN_TP_N]), while PTP averaged [HUMAN_PTP_MEAN]% ± [HUMAN_PTP_SD]% (n = [HUMAN_PTP_N]); in mouse, TP averaged [MOUSE_TP_MEAN]% ± [MOUSE_TP_SD]% (n = [MOUSE_TP_N]) and PTP [MOUSE_PTP_MEAN]% ± [MOUSE_PTP_SD]% (n = [MOUSE_PTP_N])
>>>

### Update IQR values
<source>: `/figs/figure-04/tables/fig-s5.tsv`

<old>:
<<<
PTP show much higher dispersion than TP (human IQR 0.866 vs. 0.016; mouse IQR 0.399 vs. 0.021)
>>>

<new>:
<<<
PTP show much higher dispersion than TP (human IQR [HUMAN_PTP_IQR] vs. [HUMAN_TP_IQR]; mouse IQR [MOUSE_PTP_IQR] vs. [MOUSE_TP_IQR])
>>>

---

## Figure 5 Updates

### Update single replicate metrics
<source>: `/figs/figure-05/tables/fig-5b-5c.tsv`

<old>:
<<<
brain sensitivity was 81.3% (95% CI, 74.5–88.0) with PDR 93.1% (95% CI, 89.9–96.4) and precision 67.1% (95% CI, 61.4–72.8; FDR 32.9%); kidney sensitivity was 71.9% (95% CI, 66.4–77.4) with PDR 90.6% (95% CI, 85.9–95.4) and precision 65.2% (95% CI, 61.0–69.4; FDR 34.8%)
>>>

<new>:
<<<
brain sensitivity was [BRAIN_SN_1REP]% (95% CI, [BRAIN_SN_CI_LOW]–[BRAIN_SN_CI_HIGH]) with PDR [BRAIN_PDR_1REP]% (95% CI, [BRAIN_PDR_CI_LOW]–[BRAIN_PDR_CI_HIGH]) and precision [BRAIN_PRE_1REP]% (95% CI, [BRAIN_PRE_CI_LOW]–[BRAIN_PRE_CI_HIGH]; FDR [BRAIN_FDR_1REP]%); kidney sensitivity was [KIDNEY_SN_1REP]% (95% CI, [KIDNEY_SN_CI_LOW]–[KIDNEY_SN_CI_HIGH]) with PDR [KIDNEY_PDR_1REP]% (95% CI, [KIDNEY_PDR_CI_LOW]–[KIDNEY_PDR_CI_HIGH]) and precision [KIDNEY_PRE_1REP]% (95% CI, [KIDNEY_PRE_CI_LOW]–[KIDNEY_PRE_CI_HIGH]; FDR [KIDNEY_FDR_1REP]%)
>>>

### Update two replicate metrics
<source>: `/figs/figure-05/tables/fig-5b-5c.tsv`

<old>:
<<<
precision increased to 87.8% (95% CI, 86.5–89.1) in brain and 81.0% (95% CI, 79.2–82.8) in kidney (FDR 12.2% and 19.0%, respectively), while sensitivity was 81.3% (95% CI, 79.4–83.1) and 70.6% (95% CI, 67.8–73.5)
>>>

<new>:
<<<
precision increased to [BRAIN_PRE_2REP]% (95% CI, [BRAIN_PRE_2REP_CI_LOW]–[BRAIN_PRE_2REP_CI_HIGH]) in brain and [KIDNEY_PRE_2REP]% (95% CI, [KIDNEY_PRE_2REP_CI_LOW]–[KIDNEY_PRE_2REP_CI_HIGH]) in kidney (FDR [BRAIN_FDR_2REP]% and [KIDNEY_FDR_2REP]%, respectively), while sensitivity was [BRAIN_SN_2REP]% (95% CI, [BRAIN_SN_2REP_CI_LOW]–[BRAIN_SN_2REP_CI_HIGH]) and [KIDNEY_SN_2REP]% (95% CI, [KIDNEY_SN_2REP_CI_LOW]–[KIDNEY_SN_2REP_CI_HIGH])
>>>

### Update 3-5 replicate metrics
<source>: `/figs/figure-05/tables/fig-5b-5c.tsv`

<old>:
<<<
With three to five replicates, brain precision stabilized near 88.7–89.3%, at stable sensitivity of 83.1-84.4%). In kidney, sensitivity increased with replicate number—from 70.6% (two replicates) to 78.1% with five, while precision plateaued at approximately 79–81%
>>>

<new>:
<<<
With three to five replicates, brain precision stabilized near [BRAIN_PRE_3REP]–[BRAIN_PRE_5REP]%, at stable sensitivity of [BRAIN_SN_3REP]-[BRAIN_SN_5REP]%). In kidney, sensitivity increased with replicate number—from [KIDNEY_SN_2REP]% (two replicates) to [KIDNEY_SN_5REP]% with five, while precision plateaued at approximately [KIDNEY_PRE_RANGE]%
>>>

### Update tissue-specific comparison
<source>: `/figs/figure-05/tables/fig-5b-5c.tsv`

<old>:
<<<
the brain (65 genes) and kidney (46 genes) panels produced benchmarking metrics that were minimally lower than those of the universal set in all cases, with cosine similarity between universal and tissue results of 0.9999 for brain and 0.9992 for kidney
>>>

<new>:
<<<
the brain ([BRAIN_TISSUE_GENES] genes) and kidney ([KIDNEY_TISSUE_GENES] genes) panels produced benchmarking metrics that were minimally lower than those of the universal set in all cases, with cosine similarity between universal and tissue results of [BRAIN_COSIM] for brain and [KIDNEY_COSIM] for kidney
>>>

---

## Supplementary Figure Updates

### Figure S2: AlphaGenome scores
<source>: `/figs/supp-fig-02/tables/fig-s2-log.tsv`

*Extract median expression and splicing scores for TUSCO vs other single-isoform genes*

### Figure S4: TUSCO-novel performance
<source>: `/figs/supp-fig-04/tables/fig-s4.tsv`, `fig-s4_panels.tsv`

*Extract TP, PTP, FP, FN counts per tool and sample*

### Figure S5: Coverage statistics
<source>: `/figs/figure-04/tables/fig-s5.tsv`

<old>:
<<<
Human: TP 94.5% (median 99.8%, n=5,694); PTP 63.0% (median 19.5%, n=320). Mouse: TP 96.7% (median 98.7%, n=4,549); PTP 75.6% (median 83.8%, n=2,555)
>>>

<new>:
<<<
Human: TP [HUMAN_TP_MEAN]% (median [HUMAN_TP_MEDIAN]%, n=[HUMAN_TP_N]); PTP [HUMAN_PTP_MEAN]% (median [HUMAN_PTP_MEDIAN]%, n=[HUMAN_PTP_N]). Mouse: TP [MOUSE_TP_MEAN]% (median [MOUSE_TP_MEDIAN]%, n=[MOUSE_TP_N]); PTP [MOUSE_PTP_MEAN]% (median [MOUSE_PTP_MEDIAN]%, n=[MOUSE_PTP_N])
>>>

### Figure S6: SIRV spike-in metrics
<source>: `/figs/figure-05/tables/fig-s6.tsv`, `fig-s6_bars.tsv`

<old>:
<<<
brain sensitivity was 86.1% (83.7–88.5) at 1 replicate and ~85.5% at 3–5; kidney sensitivity ranged 80.3% (78.4–82.1; 2 replicates) to 87.0% (5 replicates). F1 remained stable (brain ~92.2%; kidney ~89.0–93.0%)
>>>

<new>:
<<<
brain sensitivity was [BRAIN_SIRV_SN_1REP]% ([BRAIN_SIRV_CI]) at 1 replicate and ~[BRAIN_SIRV_SN_3-5]% at 3–5; kidney sensitivity ranged [KIDNEY_SIRV_SN_2REP]% ([KIDNEY_SIRV_CI_2REP]) to [KIDNEY_SIRV_SN_5REP]% (5 replicates). F1 remained stable (brain ~[BRAIN_F1]%; kidney ~[KIDNEY_F1_RANGE]%)
>>>

---

## Supplementary Tables

### Table S1: Cosine similarity values
<source>: `/figs/figure-03/tables/table_s1.csv`

*Update all cosine similarity values for each pipeline*

### Table S2: Read counts
<source>: Extract from `/figs/figure-05/tables/` or pipeline outputs

*Update FLNC HiFi read counts for samples K31-K35 and B31-B35*

---

## Implementation Notes

1. **Parse each table file** to extract exact numerical values
2. **Replace placeholders** (e.g., [HUMAN_TP_MEAN]) with actual values from tables
3. **Maintain precision** - use same decimal places as in table files
4. **Cross-validate** - ensure consistency between main text and supplementary materials
5. **Document source** - note which table file provided each value

## Validation Checklist

- [ ] All tissue gene counts updated from statistics files
- [ ] All performance metrics updated from figure tables
- [ ] All statistical test p-values updated
- [ ] All correlation coefficients updated
- [ ] All confidence intervals updated
- [ ] Supplementary figure values match main text where referenced
- [ ] Table S1 and S2 fully updated