# Rare Disease Variant Interpretation Pipeline V2(3): 

A reproducible Google Colab / Jupyter pipeline for annotating, prioritizing, and clinically
reporting genetic variants across **four rare Mendelian/mitochondrial diseases**. It combines
standard SNV annotation with disease-specific workflows for repeat expansions, mitochondrial
heteroplasmy, and structural variants — because none of these four diseases can be correctly
diagnosed with a generic SNV-only pipeline.

**Reference build:** GRCh38/hg38 throughout. Gene coordinates in the notebook are
placeholders for the four target regions — replace with clinically curated coordinates
before any production/clinical use.

---

## 1. What disease data is being generated, and why these four

The notebook targets **four rare diseases chosen specifically because each one defeats a
different part of a "standard" SNV-calling pipeline**. Together they exercise every major
variant class a clinical genomics pipeline needs to handle:

| Disease | Gene | Inheritance | Variant class it depends on | Why a plain SNV caller fails here |
|---|---|---|---|---|
| **Myotonic Dystrophy Type 1 (DM1)** | *DMPK* | Autosomal Dominant | CTG trinucleotide repeat expansion | Repeat expansions aren't "variants" in the SNV sense — VCF callers can't size a repeat tract; you need a specialized repeat-genotyping tool |
| **MERRF** | *MT-TK* (mitochondrial) | Mitochondrial | Point mutation at known hotspots, with **heteroplasmy** (mixed mutant/wild-type mitochondrial copies) | The pathogenic signal is the *fraction* of mitochondrial copies carrying the mutation, not a simple genotype call — needs pileup-level allele-fraction estimation |
| **Duchenne Muscular Dystrophy (DMD)** | *DMD* | X-linked Recessive | Exon-level deletions/duplications (structural variants) across a 79-exon gene | Most pathogenic *DMD* variants are whole-exon CNVs; SNV callers simply don't see them. Whether a deletion is in-frame vs. frameshift is also the key modifier between severe DMD and milder Becker phenotype |
| **Haemophilia A** | *F8* | X-linked Recessive | Full spectrum: SNVs, large deletions, and the recurrent **intron 22/intron 1 inversion** | The intron 22 inversion (found in ~45% of severe cases) arises from intrachromosomal homologous recombination and produces a signature (discordant/split reads) that standard SV callers routinely miss or mis-call as a simple deletion |

Because no single caller handles all of this, the pipeline is architected as **one shared
core (QC → region filtering → functional annotation → prioritization → ACMG classification →
reporting)** plus **four disease-specific branches**, each adding the one assay that actually
finds that disease's pathogenic variant class (Section 7 of the notebook).

The data itself, by default, is **synthetic**: Section 12 generates a small, realistic VCF
per disease (a genuinely pathogenic variant + a VUS + a benign SNP for each) so the whole
pipeline can be smoke-tested end-to-end without real patient data or multi-gigabyte
reference downloads. Real VCF/BAM input can be substituted at Section 3.

---

## 2. Tools and databases used, and why each one is there

### 2.1 Core system / SNV utilities (Section 1.1)
| Tool | Role |
|---|---|
| **bcftools** | VCF stats, normalization (`bcftools norm`), region filtering |
| **samtools** | FASTA indexing/lookup, BAM handling |
| **tabix / bgzip (htslib)** | Indexing and compressing VCFs for random-access lookups (required by nearly every downstream tool) |
| **bedtools** | Interval operations against the per-disease gene BED files |

### 2.2 Functional annotation (Section 6)
| Tool / resource | Role | Why it's needed |
|---|---|---|
| **VEP (Ensembl Variant Effect Predictor)** | Predicts consequence (missense, stop-gained, splice, etc.) per transcript | Industry-standard, transcript-aware consequence caller; needed to know *what* a variant does to the protein/transcript |
| **SnpEff** (v5.2 + SnpSift, GRCh38.99 genome) | A second, independent consequence annotator, and (via SnpSift) ClinVar annotation | Cross-checking VEP against an independently-implemented annotator catches tool-specific edge cases; SnpSift is also used to pull in ClinVar significance |
| **ClinVar** | Database of clinically reported variant classifications | Supplies existing clinical evidence — the single strongest signal for ACMG classification when a variant has been seen before |
| **gnomAD** | Population allele-frequency database | Rare-disease pathogenic variants must be rare in the general population; gnomAD AF is used both as an ACMG evidence input and as a prioritization filter |
| **SpliceAI** | Deep-learning splice-impact predictor | Needed because none of DM1/MERRF/DMD/Haemophilia A's pathology is purely missense — cryptic splice disruption is a recognized mechanism in *DMPK*, *F8*, and *DMD* |
| **AlphaMissense** (precomputed score table) | Missense pathogenicity predictor | Merged by chrom/pos/ref/alt rather than recomputed — provides an orthogonal in-silico pathogenicity estimate for missense calls |
| **REVEL** (precomputed score table) | Ensemble missense pathogenicity predictor | A second, differently-trained missense predictor, used alongside AlphaMissense so no single model's blind spots drive the classification |
| **CADD** (precomputed / tabix lookup) | Genome-wide deleteriousness score | General-purpose severity score used as a prioritization input, independent of variant class |

### 2.3 ACMG/AMP classification (Section 9)
| Tool | Role |
|---|---|
| **InterVar** | Automates ACMG/AMP evidence-code assignment (PVS1, PS1–4, PM1–6, PP1–5, BA1, BS1–4, BP1–7) from VCF input |
| **Custom point-based summarizer** (`acmg_rules.yaml`) | A transparent, auditable re-scoring layer over InterVar's evidence codes — lets a reviewer see exactly which codes drove a classification, as a cross-check on InterVar's own internal logic |

### 2.4 DM1-specific: repeat expansion (Section 7a)
| Tool | Role |
|---|---|
| **ExpansionHunter** | Genotypes short tandem repeats directly from short-read BAM alignments using a variant-catalog JSON for *DMPK*'s CTG repeat |
| **TRGT** | Same purpose, purpose-built for PacBio HiFi long reads, which resolve very large expansions (>150 repeats) that short reads struggle with |

Repeat length is then classified against `disease_rules.yaml` thresholds: ≤34 repeats
(Normal), 35–49 (Premutation), ≥50 (Pathogenic full expansion).

### 2.5 MERRF-specific: mitochondrial hotspot & heteroplasmy (Section 7b)
| Tool | Role |
|---|---|
| **cyvcf2** | Fast VCF parsing to check known *MT-TK* hotspot positions (m.8344A>G, m.8356T>C) |
| **pysam** (BAM pileup) | Direct base-counting at the hotspot position to compute **heteroplasmy fraction** — the proportion of mitochondrial reads carrying the mutant allele, which is the clinically actionable number (compared against a 10% pathogenic threshold) |

### 2.6 DMD-specific: structural variants / exon CNV (Section 7c)
| Tool | Role |
|---|---|
| **pbsv** | PacBio's structural-variant caller, for long-read exon-level deletion/duplication detection |
| **Sniffles2** | A second, independent long-read SV caller, for cross-validation against pbsv |
| **Truvari** | Benchmarks/reconciles SV call sets (e.g. pbsv vs. Sniffles2 concordance, or against a truth set) |
| **DMD exon coordinate map + reading-frame (phase) table** | Custom logic mapping SV breakpoints onto the 79 *DMD* exons and predicting whether a deletion is in-frame (milder, Becker-like) or frameshift (severe, Duchenne) — the single biggest phenotype-determining factor for this gene |
| **ClassifyCNV** | ACMG-based CNV pathogenicity classifier |
| **SvAnna** | Automatable stand-in for AnnotSV (see below) |
| **ClinGen dosage sensitivity** | Gene-level haploinsufficiency/triplosensitivity evidence, feeding CNV interpretation |

**Why AnnotSV was replaced with SvAnna:** AnnotSV is the field-standard SV annotator, but it
has **no scriptable public API** — only a manual web form
(https://www.lbgi.fr/AnnotSV/). Since this pipeline needs to run unattended, **SvAnna** was
installed instead as an automatable, command-line-drivable equivalent, working alongside
ClassifyCNV and ClinGen dosage sensitivity to cover the same interpretive ground. The manual
AnnotSV web tool is still noted as an optional cross-check step for a human reviewer.

### 2.7 Haemophilia A-specific: SNV + deletion + inversion (Section 7d)
| Tool / logic | Role |
|---|---|
| Standard SNV path (Section 6, reused) | Point mutations in *F8*, scoped to the gene region |
| `find_large_deletions_in_gene` (bcftools) | Pulls SV-VCF records overlapping *F8* |
| Custom inversion heuristic | Flags SV records (BND/INV/complex) whose breakpoints fall inside the known intron 22 (int22h-1/2/3) or intron 1 (int1h-1/2) homologous repeat regions — the recombination hotspots that produce the classic *F8* inversions. This is a **candidate flag for manual/long-read confirmation**, not a definitive call, since standard short-read SV callers routinely misrepresent these events |

### 2.8 Explicitly *not* used
- **MAVERICK** — intentionally left out as a default dependency. REVEL, AlphaMissense, CADD,
  ClinVar, and InterVar already cover pathogenicity prediction; MAVERICK's heavier setup is
  left as an optional addition for anyone who wants it.
- **AnnotSV** — see 2.6 above; replaced with SvAnna for scriptability.

### 2.9 Reporting stack (Section 11–12)
| Tool | Role |
|---|---|
| **pandas / numpy** | Variant table manipulation and scoring |
| **matplotlib / seaborn** | QC and summary figures (variant counts, ACMG breakdown, CTG repeat sizes, tool-agreement plots) |
| **Jinja2** | Templates the HTML clinical report |
| **WeasyPrint** and **ReportLab** | Two independent HTML→PDF / native-PDF export paths for the clinic-ready report |
| **BeautifulSoup4** | HTML post-processing for report assembly |

---

## 3. Expected output

Running the full pipeline (real data or the built-in synthetic generator) produces, per
disease and combined:

1. **QC summary** — `bcftools stats`-derived variant counts, Ti/Tv ratio, PASS rate, and
   QUAL/DP distribution plots, computed before any filtering.
2. **Region-filtered VCFs** — one per disease, scoped to the *DMPK*, *MT-TK*, *DMD*, and *F8*
   regions respectively (`{disease}.filtered.vcf.gz` + tabix index).
3. **Annotated variant tables** — every variant merged with VEP/SnpEff consequence, ClinVar
   significance, gnomAD AF, REVEL/AlphaMissense/CADD/SpliceAI scores, and (where relevant)
   InterVar ACMG class and SV/CNV class.
4. **Disease-specific findings**:
   - DM1: CTG repeat count and Normal/Premutation/Pathogenic classification per sample.
   - MERRF: hotspot match (m.8344A>G / m.8356T>C) and heteroplasmy fraction with a
     pathogenic/subthreshold call.
   - DMD: which exons are deleted/duplicated and whether the reading frame is predicted to
     stay in-frame (BMD-like) or shift (DMD-like).
   - Haemophilia A: SNV/deletion calls plus any candidate intron 22/1 inversion flags for
     manual review.
5. **Prioritized variant tables** (`{disease}_prioritized_variants.csv` and
   `all_diseases_prioritized_variants.csv`) — every variant given an interpretable numeric
   priority score built from gene match, functional impact, population rarity, in-silico
   pathogenicity scores, and ClinVar evidence (weights defined in `thresholds.yaml`).
6. **ACMG/AMP classification** — InterVar's evidence codes plus the notebook's own
   point-based re-scoring (`acmg_rules.yaml`), giving Pathogenic / Likely Pathogenic /
   Uncertain Significance / Likely Benign / Benign per variant.
7. **Cross-tool comparison table** (`tool_comparison_all_diseases.csv`) — ClinVar, InterVar,
   and the SV/CNV path (ClassifyCNV + SvAnna + ClinGen dosage) placed side by side per
   variant, with a **Concordant / Discordant** flag — this is the artifact a reviewer
   actually needs before signing off a classification.
8. **Per-tool output files** (`reports/by_tool/`) — each tool's calls saved individually
   (VEP, SnpEff, ClinVar/SnpSift, gnomAD, REVEL, AlphaMissense, CADD, SpliceAI, InterVar,
   ClassifyCNV/SvAnna, ClinGen dosage, MERRF heteroplasmy, DMD exon/frame mapping), both
   combined across diseases and split per disease — mirroring what each tool would emit
   standalone.
9. **Figures** — variant counts per disease, ACMG classification breakdown (pie chart),
   DM1 CTG repeat-size plot with normal/pathogenic threshold lines, tool-agreement summary.
10. **Clinic-ready reports** — one **HTML + PDF** report per disease, plus a combined
    **master report** with the full cross-tool comparison appendix, styled with
    color-coded ACMG classification rows (pathogenic = red, VUS = yellow, benign = green).

All outputs land under `/content/RareDiseasePipeline/` in Colab, organized into `refs/`,
`beds/`, `inputs/`, `work/`, `reports/` (including `reports/by_tool/`), and `figures/`.

---

## 4. How to use

### Quickest path — smoke test with synthetic data (no real files needed)
1. Open the notebook in **Google Colab** (recommended — it relies on `apt-get`, `conda`, and
   Colab-friendly paths; a local Jupyter/WSL environment will need path adjustments).
2. Run **Section 1 (Environment Setup)** top to bottom, once per Colab session. This installs
   bcftools/samtools/bedtools, SnpEff, VEP (Perl/htslib dependencies included), SvAnna,
   ClassifyCNV, InterVar, sniffles2/truvari/pbsv (via conda), and ExpansionHunter/TRGT.
   This step is slow (many downloads/compiles) — expect it to take a while, and re-run
   individual cells if a download flakes rather than restarting the whole section, since most
   install cells are idempotent (they check `if not exists`).
3. Run **Section 2 (Configuration)** — writes `disease_rules.yaml`, `acmg_rules.yaml`, and
   `thresholds.yaml`, and loads them into `DISEASES`, `ACMG_RULES`, `THRESHOLDS`.
4. Skip Section 3 (real upload) and go straight to **Section 12.1** to generate the four
   synthetic per-disease VCFs, then run **Section 12.2** to bgzip/tabix-index them.
5. Continue through **Sections 12.3–12.4** to build the cross-tool comparison table and
   prioritized variant tables from the illustrative annotation values.
6. Run the report-generation cells in **Section 11/12.6** to produce the HTML/PDF outputs.

### Full path — real patient data
1. Complete Section 1 (Environment Setup) as above.
2. In **Section 1.11 (Reference Data Setup)**, replace the placeholder paths with real
   downloads: a GRCh38 FASTA, the full VEP GRCh38 cache, the SnpEff `GRCh38.99` genome
   (`java -jar snpEff/snpEff.jar download GRCh38.99`), a ClinVar GRCh38 VCF, and a gnomAD VCF.
3. In **Section 3**, upload (or place under `DIRS["inputs"]`) your small-variant VCF, an SV
   VCF from a long-read caller, an ExpansionHunter/TRGT repeat VCF, per-disease BED files,
   and — for MERRF heteroplasmy and DMD/Haemophilia A SV calling — a BAM/CRAM.
4. Before trusting results clinically: replace the **illustrative GRCh38 coordinates** in
   `disease_rules.yaml`, curate a real **DMD exon BED (79 exons, NM_004006.3)** with an
   accurate exon-phase table, and supply a real **ExpansionHunter variant-catalog JSON** and
   **TRGT repeat BED** for the *DMPK* CTG repeat.
5. Run Sections 4 → 10 in order (QC → region filtering → functional annotation →
   disease-specific branches → prioritization → ACMG classification → visualizations), then
   Section 11 for the report. Section 6-9 outputs (the real merged annotation dataframe) can
   be substituted directly into Section 12's downstream logic — the scoring, tool-comparison,
   and reporting code all run unchanged against real data.

### Notes and caveats
- The **"Running the Real Analysis Tools on Synthetic VCFs"** and later cells (Sections
  89+) reflect an iterative debugging history (VEP htslib install flag flip-flopping,
  `cyvcf2`/`weasyprint`/`fpdf2` install fixes, orchestrator-cell re-triggers) rather than a
  clean linear script — expect some redundant/superseded cells retained for provenance. The
  canonical linear pipeline is Sections 1–12.
- Full-scale reference downloads (GRCh38 FASTA, VEP cache, ClinVar, gnomAD) are large; the
  synthetic-data smoke test avoids all of them and is the fastest way to confirm the
  pipeline's logic end-to-end before committing to a full run.
- AnnotSV's manual web form (https://www.lbgi.fr/AnnotSV/) remains a recommended optional
  cross-check for DMD/Haemophilia A structural calls, alongside the automated SvAnna path.
