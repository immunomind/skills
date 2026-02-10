---
name: immundata
description: Use this skill for AIRR-seq (Adaptive Immune Receptor Repertoire / VDJ-seq) data analysis with immunarch + immundata in R, including ingestion, receptor schema design, immutable transformations, clonality/diversity/public overlap metrics, and Seurat/AnnData integration.
---

# AIRR data analysis with immundata R package

## Overview

This skill provides a practical workflow for Adaptive Immune Receptor Repertoire (AIRR / VDJ-seq) data analysis with **immunarch** (analysis + visualization) and **immundata** (data ingestion + transformation + schema handling).

Use it to:
- Load AIRR data from one file, many files, glob patterns, or metadata tables.
- Define receptor semantics (chain-agnostic, single-chain, paired-chain).
- Build reproducible immutable pipelines.
- Compute repertoire statistics, clonality, diversity, and public overlap.
- Move annotations between ImmunData and single-cell objects (e.g., Seurat/AnnData).

---

## When to Use This Skill

Use this skill when the user asks to:

- Analyze **bulk** or **single-cell** AIRR/TCR/BCR data.
- Compare repertoires across sample groups (tissue, therapy, cluster, donor, timepoint).
- Compute clonality/diversity/publicity metrics.
- Define or change receptor schema (e.g., `cdr3_aa + v_call`, TRA-only, TRA+TRB).
- Filter receptors by patterns or sequence distance.
- Add/propagate labels between repertoire data and scRNA metadata.
- Convert old immunarch objects to the newer ImmunData pipeline.

---

## Quick Start

### Basic import and first look

```r
library(immunarch)

# Demo data + basic grouping
idata <- get_test_idata() |> agg_repertoires("Therapy")
idata

# Core analyses
airr_stats_genes(idata, gene_col = "v_call") |> vis()
airr_public_jaccard(idata) |> vis()
airr_clonality_prop(idata)
airr_diversity_pielou(idata) |> vis()
````

### Optional: add clonality labels to Seurat metadata

```r
idata <- annotate_clonality_prop(idata)
sdata <- annotate_seurat(idata, sdata, cols = "clonal_prop_bin")
Seurat::DimPlot(sdata, reduction = "umap", group.by = "clonal_prop_bin", shuffle = TRUE)
```

### Ingest AIRR files via immundata directly

```r
library(immundata)

md_path <- system.file("extdata/tsv", "metadata.tsv", package = "immundata")
files <- c(
  system.file("extdata/tsv", "sample_0_1k.tsv", package = "immundata"),
  system.file("extdata/tsv", "sample_1k_2k.tsv", package = "immundata")
)
md <- read_metadata(md_path)

idata <- read_repertoires(
  path     = files,
  schema   = c("cdr3_aa", "v_call"),
  metadata = md
)
```

---

## Typical User Intake (what to ask/assume)

Before coding, identify:

1. **Modality**: bulk vs single-cell AIRR.
2. **Input format**: TSV/CSV/Parquet, gzipped or not, one file vs many.
3. **Schema intent**:

   * chain-agnostic,
   * single-chain (e.g., TRA only),
   * paired-chain (e.g., TRA+TRB, IGH + IGK|IGL).
4. **Grouping variables**: repertoire schema (sample/cluster/tissue/therapy).
5. **Target analyses**: stats, gene usage, clonality, diversity, overlap, annotation transfer.
6. **Scale/performance**: whether snapshots/materialization strategy is needed.

If unknown, default to conservative, reproducible choices and print intermediate summaries.

---

## Standard Analysis Workflow

### 1) Ingest data and metadata

Use `read_metadata()` + `read_repertoires()`.

```r
library(immundata)

inp_files <- "path/to/airr/*.tsv.gz"
md_table  <- read_metadata("path/to/metadata.tsv")

idata <- read_repertoires(
  path              = inp_files,
  schema            = c("cdr3_aa", "v_call"),
  metadata          = md_table,
  repertoire_schema = "Sample"
)
```

### 2) Choose receptor schema explicitly

#### Chain-agnostic (bulk / pre-filtered)

```r
idata <- read_repertoires(
  path   = inp_files,
  schema = c("cdr3_aa", "v_call")
)
```

#### Single-chain (e.g., TRA only)

```r
schema <- make_receptor_schema(
  features = c("cdr3", "v_call"),
  chains   = "TRA"
)

idata <- read_repertoires(
  path        = "path/to/single_cell.csv.gz",
  schema      = schema,
  barcode_col = "barcode",
  locus_col   = "locus"
)
```

#### Paired-chain (e.g., TRA + TRB)

```r
schema <- make_receptor_schema(
  features = c("cdr3", "v_call"),
  chains   = c("TRA", "TRB")
)

idata <- read_repertoires(
  path        = "path/to/single_cell.csv.gz",
  schema      = schema,
  barcode_col = "barcode",
  locus_col   = "locus",
  umi_col     = "umis"
)
```

#### Paired with alternative light chain (e.g., IGH + IGK|IGL)

```r
schema <- make_receptor_schema(
  features = c("cdr3", "v_call"),
  chains   = c("IGH", "IGK|IGL")
)
```

### 3) Transform immutably (filter / annotate / mutate)

```r
# Filtering
idata_f <- idata |>
  filter(v_call == "TRBV2") |>
  filter(imd_proportion >= 0.0002)

# Sequence-distance filter
idata_seq <- idata |>
  filter(seq_options = make_seq_options(
    patterns  = "CASSELAGYRGEQYF",
    query_col = "cdr3",
    method    = "lev",
    max_dist  = 3
  ))

# Annotation join
idata_ann <- annotate(
  idata       = idata,
  annotations = cells[c("barcode", "ident")],
  by          = c(imd_barcode = "barcode"),
  keep_repertoires = FALSE
)

# Mutations
idata_mut <- idata |>
  mutate(big_chains = umis >= 10) |>
  mutate(dist_to_pattern = dd$levenshtein(cdr3, "CASSSVSGNSPLHF"))
```

### 4) Aggregate repertoires for reporting strata

```r
idata_grp <- idata |>
  agg_repertoires("Tissue")
```

### 5) Compute immune repertoire statistics

```r
# Chain and gene statistics
chains <- airr_stats_chains(idata_grp)
genes  <- airr_stats_genes(idata_grp, gene_col = "v_call", level = "receptor")

chains |> vis()
genes  |> vis()
```

### 6) Clonality analysis

```r
cl_line <- airr_clonality_line(idata_grp)
cl_prop <- airr_clonality_prop(idata_grp)
cl_rank <- airr_clonality_rank(idata_grp, bins = c(10, 100, 1000))

cl_prop |> vis()
cl_rank |> vis()
```

### 7) Diversity analysis

```r
d50      <- airr_diversity_dxx(idata_grp, perc = 50)
chao1    <- airr_diversity_chao1(idata_grp)
shannon  <- airr_diversity_shannon(idata_grp)
pielou   <- airr_diversity_pielou(idata_grp)
hill1    <- airr_diversity_index(idata_grp)
hillprof <- airr_diversity_hill(idata_grp, q = c(0, 1, 2))

pielou |> vis()
```

### 8) Public overlap

```r
m_pub <- airr_public_intersection(idata_grp)
m_jac <- airr_public_jaccard(idata_grp)

m_pub |> vis()
m_jac |> vis()
```

### 9) Snapshot expensive steps

```r
# Save intermediate immutable snapshot to avoid recomputing expensive transforms
idata_cached <- immundata::write_immundata(idata_mut, "path/to/snapshot_folder")
```

---

## Common Pitfalls and Best Practices

1. **Changing receptor definition mid-analysis**

   * Receptor schema is foundational. If receptor definition changes, re-read data into a new `ImmunData`.

2. **Using wrong mode for single-cell**

   * Paired analysis requires `barcode_col`, `locus_col`, and typically `umi_col`.
   * Single-chain mode does not necessarily collapse multiple chains per barcode.

3. **Skipping metadata strategy**

   * Prefer metadata-driven ingestion (`path = "<metadata>"`) when file provenance and sample mapping matter.

4. **Direct internal-table edits**

   * Avoid low-level manual edits. Use high-level verbs (`filter`, `annotate`, `mutate`, `agg_repertoires`).

5. **Ignoring immutable pipeline behavior**

   * Each step returns a new object; persist expensive steps with snapshots.

6. **Non-canonical columns**

   * In scripts, rely on canonical `imd_*` columns.
   * In package code, prefer schema keys/aliases (e.g., via `imd_schema()`).

7. **Distance-heavy transforms on large data**

   * Pattern/Levenshtein operations can be expensive; run once, snapshot, then reuse.

---

## Bundled Resources

TBD

---

## Additional Resources

* immunomind docs home: [https://immunomind.github.io/docs/](https://immunomind.github.io/docs/)
* Quick Start: [https://immunomind.github.io/docs/intro/quick_start/](https://immunomind.github.io/docs/intro/quick_start/)
* Reading repertoire files: [https://immunomind.github.io/docs/guides/io/ingesting/](https://immunomind.github.io/docs/guides/io/ingesting/)
* Reading single-/paired-chain data: [https://immunomind.github.io/docs/guides/io/modes/](https://immunomind.github.io/docs/guides/io/modes/)
* Filter: [https://immunomind.github.io/docs/guides/transform/filter/](https://immunomind.github.io/docs/guides/transform/filter/)
* Annotate: [https://immunomind.github.io/docs/guides/transform/annotate/](https://immunomind.github.io/docs/guides/transform/annotate/)
* Mutate: [https://immunomind.github.io/docs/guides/transform/mutate/](https://immunomind.github.io/docs/guides/transform/mutate/)
* Workflow phase 1 (ingestion): [https://immunomind.github.io/docs/concepts/workflow/phase_ingestion/](https://immunomind.github.io/docs/concepts/workflow/phase_ingestion/)
* Workflow phase 2 (transformation): [https://immunomind.github.io/docs/concepts/workflow/phase_transformation/](https://immunomind.github.io/docs/concepts/workflow/phase_transformation/)
* ImmunData structure: [https://immunomind.github.io/docs/guides/data_schema/](https://immunomind.github.io/docs/guides/data_schema/)
* Single-cell end-to-end tutorial: [https://immunomind.github.io/docs/tutorials/single_cell/](https://immunomind.github.io/docs/tutorials/single_cell/)

---

## Tips for Effective Analysis

* Start with a **small subset** and verify schema + grouping before scaling up.
* Print object summaries after ingestion and after major transforms.
* Use explicit variable names for stages (`idata_raw`, `idata_qc`, `idata_ann`, `idata_stats`).
* Prefer pipelines that can be re-executed end-to-end from raw inputs.
* Keep biologically meaningful grouping variables in repertoire schema early.
* Use `vis()` early and often for sanity checks before formal interpretation.
