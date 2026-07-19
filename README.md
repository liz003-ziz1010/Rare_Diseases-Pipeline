# README.md: Rare Disease Variant Interpretation Pipeline

This notebook implements a comprehensive, reproducible pipeline for the annotation, prioritization, and clinical reporting of variants in four rare diseases: Myotonic Dystrophy Type 1 (DM1), MERRF, Duchenne Muscular Dystrophy (DMD), and Haemophilia A.

It combines standard short-read SNV annotation with disease-specific workflows for structural variants, repeat expansions, and mitochondrial heteroplasmy.

## Table of Contents

1.  [Overview](#overview)
2.  [Diseases Covered](#diseases-covered)
3.  [Key Features](#key-features)
4.  [How to Use](#how-to-use)
5.  [Pipeline Structure](#pipeline-structure)
6.  [Outputs](#outputs)
7.  [Notes and Limitations](#notes-and-limitations)

---

## Overview

This pipeline provides a robust framework for identifying and classifying pathogenic variants in genetically complex rare diseases. It integrates various bioinformatics tools and custom scripts to process VCFs and BAMs, generate comprehensive annotations, prioritize variants based on ACMG guidelines and custom thresholds, and produce clinic-ready reports.

**Reference build:** GRCh38/hg38 is used throughout. Coordinates are placeholders and should be confirmed against your specific BED files for production runs.

## Diseases Covered

| Disease                         | Gene   | Variant Class Emphasized                       |
| :------------------------------ | :----- | :--------------------------------------------- |
| Myotonic Dystrophy Type 1 (DM1) | *DMPK* | CTG repeat expansion                           |
| MERRF                           | *MT-TK*| Mitochondrial point mutation / heteroplasmy    |
| Duchenne Muscular Dystrophy (DMD)| *DMD*  | Exon-level deletion/duplication (SV)           |
| Haemophilia A                   | *F8*   | SNV + large deletion + intron 22/1 inversion   |

## Key Features

*   **Environment Setup**: Automated installation of system tools, annotation engines (VEP, SnpEff), SV/repeat callers (ExpansionHunter, TRGT, pbsv, Sniffles2), and Python packages (cyvcf2, pandas, matplotlib, seaborn, weasyprint, fpdf2).
*   **Configurable Rules**: Disease definitions, ACMG evidence weights, and prioritization thresholds are defined in YAML files for easy auditing and extension.
*   **Comprehensive Annotation**: Integrates functional consequence prediction (VEP, SnpEff), clinical evidence (ClinVar via SnpSift), population frequency (gnomAD), and computational pathogenicity scores (SpliceAI, AlphaMissense, REVEL, CADD).
*   **Disease-Specific Workflows**: Specialized branches for:
    *   **DM1**: Trinucleotide repeat expansion detection and classification.
    *   **MERRF**: Mitochondrial hotspot checking and heteroplasmy estimation from BAM pileups.
    *   **DMD**: Structural variant calling, exon mapping, frame effect prediction, and CNV classification (using SvAnna and ClassifyCNV).
    *   **Haemophilia A**: Standard SNV analysis, large deletion detection, and flagging of intron 22/1 inversion candidates.
*   **Variant Prioritization**: Combines various evidence types into a single interpretable score per variant.
*   **ACMG Classification**: Automated ACMG/AMP evidence code assignment and point-based summarization via InterVar.
*   **Data Visualization**: Generates summary figures for variant counts, ACMG breakdown, repeat size distributions, and exon deletion maps.
*   **Clinical Report Generation**: Produces prioritized variant CSVs, HTML clinical summaries, and PDF exports, suitable for review.
*   **Synthetic Test Data**: Includes a section to generate realistic per-disease synthetic VCFs to smoke-test the entire pipeline.

## How to Use

1.  **Run Environment Setup**: Execute **Section 1** (Environment Setup) once per Colab session. This installs all necessary tools and Python packages.
2.  **Execute Top-to-Bottom**: After setup, run the notebook cells sequentially from top to bottom. The orchestrator (`f3f54a11`) is designed to manage dependencies and trigger subsequent cells, especially after a kernel restart.
3.  **Generate Test Data (Optional)**: If you don't have your own data, proceed to **Section 12** to generate synthetic VCFs for each disease. This allows you to test the full pipeline functionality.
4.  **Upload Input Data (Optional)**: If using real data, utilize **Section 3** (Input Upload) to upload your VCFs, BAMs, and BED files.
5.  **Review Outputs**: Inspect the generated plots, tables, and reports in the respective sections.

## Pipeline Structure

The notebook is organized into the following major sections:

*   **1. Environment Setup**: Installs system tools, bioinformatics software, and Python libraries.
*   **2. Configuration**: Defines disease-specific rules, ACMG evidence weights, and prioritization thresholds.
*   **3. Input Upload**: Provides mechanisms to upload patient VCFs, BAMs, and custom BED files.
*   **4. Quality Control**: Performs basic QC on raw VCF callsets.
*   **5. Region Filtering**: Restricts analysis to gene regions of interest.
*   **6. Functional Annotation**: Applies VEP, SnpEff, ClinVar, gnomAD, SpliceAI, AlphaMissense, REVEL, and CADD annotations.
*   **7. Disease-Specific Branches**: Executes specialized workflows for DM1, MERRF, DMD, and Haemophilia A.
*   **8. Variant Prioritization**: Scores variants based on configured rules.
*   **9. ACMG Classification**: Assigns ACMG classifications using InterVar.
*   **10. Visualizations**: Generates summary plots and figures.
*   **11. Clinical Report Generation**: Creates final HTML, PDF, and CSV reports.
*   **12. Appendix — Synthetic Test Data Generator**: Provides functions to create artificial VCFs for testing the pipeline.

## Outputs

Upon successful execution, the pipeline generates:

*   **Prioritized Variant Tables**: CSV files containing scored and ranked variants (`all_diseases_prioritized_variants.csv` and per-disease CSVs).
*   **Tool Comparison Tables**: A `tool_comparison_all_diseases.csv` file summarizing agreement/disagreement between different annotation tools.
*   **Individual Tool Outputs**: Separate CSVs for each tool's contribution, both combined and per-disease.
*   **QC Plots**: Visualizations of variant quality and depth distributions.
*   **Summary Figures**: Plots showing variant counts, ACMG breakdown, repeat sizes, and exon deletion maps.
*   **Clinical Reports**: Comprehensive HTML (`.html`) and PDF (`.pdf`) reports, combining all findings and figures, plus a raw `.csv` output of prioritized variants.

All outputs are saved within the `DIRS` structure, typically under `/content/RareDiseasePipeline/reports/` and `/content/RareDiseasePipeline/figures/`.

## Notes and Limitations

*   **AnnotSV vs. SvAnna**: Due to the lack of a scriptable API for AnnotSV, the pipeline uses SvAnna (an automatable equivalent) alongside ClassifyCNV for structural variant interpretation.
*   **MAVERICK**: This tool is intentionally left optional, as its functionality overlaps with other pathogenicity predictors included (REVEL, AlphaMissense, CADD, ClinVar, InterVar).
*   **Reference Data**: While the pipeline attempts to download some reference data (e.g., GRCh38 FASTA, VEP cache), users should replace placeholder paths with their own curated BED/JSON resources and ensure all necessary reference files (e.g., ClinVar VCF, gnomAD VCF, BAM files) are available for a full production run.
*   **SnpEff Download Reliability**: SnpEff genome downloads can sometimes be unreliable. The orchestrator includes retry logic, but manual intervention might occasionally be required.
*   **Kernel Restarts**: The pipeline is designed to handle kernel restarts, common in Colab environments, by re-executing necessary setup cells automatically via an orchestrator script.
