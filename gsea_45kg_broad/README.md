# GSEA 45kg Liver — Broad Institute GSEA (Java desktop)

Comparison: DLYL45 vs TFBL45 (n=4/group)
Gene sets: KEGG pathways for Sus scrofa (ssc.gene.kegg.gmt)
Method: Preranked GSEA, gene-level ranking by signal-to-noise ratio
Run date: 2026-05-11

## Key files

- `DLYL45-vs-TFBL45.ssc.gene.kegg.Gsea.xls` — Full results table (330 KEGG pathways)
- `DLYL45-vs-TFBL45.ssc.gene.kegg.Gsea.all.xls.gz` — All results including non-KEGG
- `ssc.gene.kegg.gmt` — Gene set database used
- `run_log.txt` — Run parameters and settings
- `all.cls` / `all.sample.cls` — Phenotype class labels

## Interpretation

- Comparison: DLYL45-vs-TFBL45, DLY = positive class
- Positive NES = enriched in DLY (higher expression in DLY liver)
- Negative NES = enriched in TFB (higher expression in TFB liver)
- FDR: Benjamini-Hochberg correction

## Key results (FDR < 0.05)

17 pathways FDR<0.05 (7 DLY-enriched + 10 TFB-enriched)

DLY-enriched (AA biosynthesis):
- Alanine, aspartate and glutamate metabolism (KO00250) NES=+2.01 FDR=0.002
- Tyrosine metabolism (KO00350) NES=+2.10 FDR=0.002
- Arginine biosynthesis (KO00220) NES=+1.96 FDR=0.002
- Cysteine and methionine metabolism (KO00270) NES=+2.08 FDR=0.002
- Biosynthesis of amino acids (KO01230) NES=+1.80 FDR=0.012
- Proteasome (KO03050) NES=+1.69 FDR=0.019
- Arginine and proline metabolism (KO00330) NES=+1.57 FDR=0.046

TFB-enriched (signaling/immune):
- Platelet activation, TCR signaling, PI3K signaling, etc.

## Comparison with previous external GSEA

This analysis replaces the clusterProfiler/KOBAS results previously used.
The NES values are similar but the Group labeling convention differs.
This analysis uses the Broad GSEA standard: DLYL45-vs-TFBL45 → positive NES = DLY-enriched.
