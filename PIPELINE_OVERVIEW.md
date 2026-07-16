# Rare Disease Variant Interpretation Pipeline

**DM1 · MERRF · DMD · Haemophilia A**

A targeted, reproducible Google Colab pipeline (Nextflow-style logic, run
notebook-first) for annotating, prioritizing, and clinically reporting
variants in four rare diseases. It combines standard short/long-read SNV
annotation with disease-specific structural variant, repeat-expansion, and
mitochondrial heteroplasmy workflows, and ends in clinic-ready HTML/PDF
reports.

| Disease | Gene | Variant class this pipeline emphasizes |
|---|---|---|
| Myotonic Dystrophy Type 1 (DM1) | *DMPK* | CTG repeat expansion |
| MERRF | *MT-TK* | Mitochondrial point mutation / heteroplasmy |
| Duchenne Muscular Dystrophy (DMD) | *DMD* | Exon-level deletion/duplication (SV) |
| Haemophilia A | *F8* | SNV + large deletion + intron 22/1 inversion |

**Reference build:** GRCh38/hg38 unless noted. Gene-region coordinates are
placeholders for these four genes — confirm against your own BED files
before a production run.

**How to run it:** execute Section 1 (Environment Setup) once per Colab
session, then run sections top-to-bottom. Section 12 (Appendix) generates
synthetic per-disease test VCFs so the whole pipeline can be smoke-tested
before pointing it at real patient data.

---

## Pipeline Sections

### 1. Environment Setup
Installs the full tool chain in the Colab VM: system packages, Ensembl VEP
(targeted config, no full offline cache by default), ISV and InterVar for
ACMG classification, `sniffles` (via conda/condacolab) and the long-read SV
callers `pbsv` / `Truvari`, the repeat-expansion callers `ExpansionHunter`
(short-read) and `TRGT` (PacBio HiFi), and supporting Python packages.
*AnnotSV is intentionally skipped in this run* — its role is filled by
`SvAnna` further down the pipeline.

### 2. Configuration
Three YAML config blocks are written and loaded, so every downstream
threshold is transparent and editable in one place:

- **`disease_rules.yaml`** — per-disease gene, chromosome, region coordinates,
  inheritance pattern, variant class, and disease-specific parameters (e.g.
  DM1's normal/premutation/pathogenic CTG-repeat cutoffs of 34/35–49/≥50;
  MERRF's known hotspot positions and 10% heteroplasmy pathogenicity
  threshold; DMD's exon count; F8's known inversion sites).
- **`acmg_rules.yaml`** — a simplified point-based mapping of the ACMG/AMP
  evidence codes (PVS1=8, PS1–4=4, PM1–6=2, PP1–5=1, and the matching
  benign-side negative weights) into five classification bands, used as a
  transparent complement to InterVar's own rule engine (not a replacement).
- **`thresholds.yaml`** — the prioritization scoring rules: gnomAD
  rare/common cutoffs, functional-impact weights (HIGH=5, MODERATE=3, LOW=1,
  MODIFIER=0), REVEL/AlphaMissense/CADD/SpliceAI pathogenicity cutoffs,
  ClinVar evidence weights, and the gene-match bonus.

### 3. Input Upload
A Colab file-upload widget for bringing in VCF/BAM inputs (falls back
gracefully when not running inside Colab).

### 4. Quality Control
Runs `bcftools stats` on the input SNV/indel VCF and plots QUAL/DP
distributions and PASS rate before anything is annotated.

### 5. Region Filtering
Writes default BED files straight from `disease_rules.yaml` (overridable),
then filters/bgzips/indexes the input VCF down to each disease's target
region — keeping every downstream step fast and disease-specific.

### 6. Functional Annotation
Normalizes variants, then annotates with **VEP** and **SpliceAI**. (SnpEff,
ClinVar/SnpSift, gnomAD, REVEL, AlphaMissense, and CADD are layered in here
and in the disease-specific branches to build the full annotation set seen
in the final tables.)

### 7. Disease-Specific Branches
Each disease gets its own targeted sub-pipeline:

- **7a. DM1** — runs **ExpansionHunter** against a DMPK repeat-catalog JSON
  to size the CTG repeat from short-read BAM data.
- **7b. MERRF** — checks calls against the known MT-TK hotspot positions and
  reports heteroplasmy fraction.
- **7c. DMD** — runs SV calling and maps calls onto a GRCh38 DMD exon
  coordinate map (79 exons) to determine which exons are deleted and whether
  the deletion is in-frame or frameshift.
- **7d. Haemophilia A** — standard SNV extraction (reusing Section 6's
  annotation functions) plus large-deletion and intron 22/1 inversion
  breakend detection, the latter flagged for manual long-read/inverse-PCR
  confirmation rather than auto-called.

### 8. Variant Prioritization
`score_variant()` computes a transparent, additive priority score per
variant: +5 for gene match, up to +5 for functional impact tier, ±5 for
population rarity/commonness in gnomAD, and small bonuses for each in-silico
predictor that crosses its pathogenicity cutoff (REVEL, AlphaMissense, CADD,
SpliceAI). Every point is logged in a human-readable `score_reasons` string
so the score is auditable.

### 9. ACMG Classification
Runs **InterVar** to apply the ACMG/AMP 2015 guidelines and assign each
variant a five-tier classification (Pathogenic → Benign) with the specific
evidence codes invoked.

### 10. Visualizations
Generates the cohort-level figures: variants-per-disease bar chart,
cross-tool concordance stacked bar chart, DMD exon-deletion map, DM1 repeat-size
bar chart with normal/pathogenic threshold lines, and the ACMG classification
pie chart.

### 11. Clinical Report Generation
Builds the clinic-ready HTML report template, then exports to PDF
(`weasyprint`, with an `fpdf2` fallback).

### 12. Realistic Per-Disease Test Data, Multi-Tool Comparison, and Clinic-Ready Reports
The appendix that makes the whole pipeline reproducible end-to-end without
real patient data:

- **12.1–12.2** — a VCF writer plus one realistic synthetic-variant generator
  per disease, smoke-tested through the same Section-5 region-filtering path
  real data would take.
- **12.3–12.4** — builds the cross-tool comparison table (ClinVar vs
  InterVar vs ClassifyCNV/SvAnna, with a Concordant/Discordant flag and a
  count of how many tools reported a call), then prioritizes each disease
  separately and assembles the combined table.
- **12.5** — the cross-tool concordance figure.
- **12.6–12.7** — the clinic-ready HTML report generator (one report per
  disease, plus the combined master report), and its PDF export.

---

## Output Artifacts

| File | Contents |
|---|---|
| `DM1_clinic_report.html`, `MERRF_clinic_report.html`, `DMD_clinic_report.html`, `HEMA_A_clinic_report.html` | One clinic-ready report per disease: findings summary, prioritized variants table, multi-tool comparison table. |
| `master_clinic_report_all_diseases.html` | All four diseases combined into a single report, plus cohort-level figures. |
| `variant_counts_per_disease.png`, `tool_agreement_by_disease.png`, `dmd_exon_map.png`, `dm1_repeat_summary.png`, `acmg_breakdown.png` | Cohort-level and disease-specific figures embedded in the reports above. |

## Design Notes / Caveats

- Gene-region coordinates in the notebook are **placeholders** — replace
  with confirmed BED files before any production run.
- **AnnotSV is skipped**; SvAnna is used in its place throughout the SV/CNV
  path. If a strict AnnotSV output is required downstream, swap it back in
  at Section 7c/7d.
- The point-based ACMG scoring in `acmg_rules.yaml` is a **transparency
  layer**, not a replacement for InterVar's own rule engine — the two are
  run in parallel and should agree; if they don't, that's exactly the kind
  of discordance the Multi-Tool Comparison table is designed to surface.
- Structural variants and repeat expansions leave REVEL/AlphaMissense/
  CADD/SpliceAI empty (not applicable) — this is expected, not missing data.
- Any variant flagged `tool_agreement = Discordant`, or any SV/inversion
  candidate reported by only one tool (`n_tools_reporting = 1`, e.g. the F8
  intron 22 inversion breakend), should be manually reviewed before clinical
  sign-out.
