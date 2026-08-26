---
name: immunarch
description: Use this skill for AIRR-seq (Adaptive Immune Receptor Repertoire / VDJ-seq) data analysis with immunarch + immundata in R, including ingestion, receptor schema design, immutable transformations, clonality/diversity/public overlap metrics, and Seurat/AnnData integration.
---

# AIRR-seq analysis with immunarch and immundata

Use immundata to ingest, represent, transform, and persist AIRR-seq data. Use immunarch to calculate and visualise repertoire-level statistics from `ImmunData`.

## Source of truth

Use current package help and examples. For API details, consult these topics before inventing arguments or function names:

```r
?immundata::read_manifest
?immundata::read_repertoires
?immundata::make_receptor_schema
?immundata::annotate_immundata
?immundata::filter_immundata
?immundata::mutate_immundata
?immundata::agg_repertoires
?immundata::write_immundata
?immunarch::airr_desc
?immunarch::airr_clonality
?immunarch::airr_diversity
?immunarch::repsim
?immunarch::annotate_clonality
```

Use only exported current interfaces demonstrated in those docs.

## When to Use This Skill

Use this skill when the user asks to:

- Analyze **bulk** or **single-cell** AIRR/TCR/BCR data.
- Compare repertoires across sample groups (tissue, therapy, cluster, donor, timepoint).
- Compute clonality/diversity/publicity metrics.
- Define or change receptor schema (e.g., `cdr3_aa + v_call`, TRA-only, TRA+TRB).
- Filter receptors by patterns or sequence distance.
- Add/propagate labels between repertoire data and scRNA metadata.
- Convert old immunarch objects to the newer ImmunData pipeline.

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

If unknown, default to conservative, reproducible choices, with data on smaller scale, so snapshotting only after crucial operations, and print intermediate summaries.

## Tips for Effective Analysis

* Start with a **small subset** and verify schema + grouping before scaling up.
* Print object summaries after ingestion and after major transforms.
* Use explicit variable names for stages (`idata_raw`, `idata_qc`, `idata_ann`, `idata_stats`).
* Prefer pipelines that can be re-executed end-to-end from raw inputs.
* Keep biologically meaningful grouping variables in repertoire schema early.
* Use `vis()` early and often for sanity checks before formal interpretation.

## Data model

- `read_repertoires()` ingests raw files into a disk-backed `ImmunData`; use it once for each receptor definition. It writes `annotations.parquet` and `metadata.json`; keep that directory while the object is in use.
- `make_receptor_schema()` defines which features and chains identify one receptor; use it when selecting one chain or pairing TRA and TRB.
- `agg_repertoires()` defines the biological comparison units and recalculates counts and proportions; use it after adding or changing grouping annotations.
- `filter()` subsets rows or whole matching receptors; use it for cohorts, genes, abundance thresholds, or sequence searches.
- `mutate()` adds lazy row- or group-level features; use it for lengths, flags, distances, and derived annotations.
- `annotate_*()` joins external labels without changing receptor identity; choose the helper matching barcodes, receptors, chains, or arbitrary keys.
- `compute()` executes pending work while keeping an `ImmunData`; use it when an expensive lazy result will be reused. `collect()` returns a tibble; use it only for a deliberately small final result.
- `write_immundata()` persists the complete state; use it at meaningful checkpoints and before ending a session.
- Core data and transformations stay in lazy duckplyr queries. Do not convert the main dataset to a data frame or tibble.
- Operations are immutable. Assign the returned `ImmunData` from every transformation and save.
- Use `ImmunData` methods and public accessors. Do not read or edit internal annotation storage directly.
- A receptor schema defines receptor identity. A repertoire schema defines the biological groups within which receptor counts and proportions are calculated. Changing either definition changes the analysis unit.

## Ingestion

Load the two packages:

```r
library(immundata)
library(immunarch)
```

### Bulk files

For bulk data, supply receptor features and the abundance column. Use this mode when rows represent clonotypes or chains rather than cells. Receptor feature names refer to columns after `rename_columns` is applied.

```r
idata <- read_repertoires(
  path = "path/to/sample.tsv.gz",
  schema = c("cdr3_aa", "v_call"),
  count_col = "counts",
  output_folder = "path/to/processed/sample"
)
```

`path` may also be a vector of same-format files or a glob. The default `repertoire_schema = "<auto>"` creates one repertoire per input file.

### Manifest-driven ingestion

A manifest has one row per repertoire file and may contain donor, tissue, treatment, or other per-file annotations. Use it for multi-sample cohorts where file provenance and grouping must stay explicit. `read_manifest()` validates and resolves file paths; `read_repertoires(path = "<manifest>")` ingests the listed files.

```r
manifest <- read_manifest("path/to/manifest.csv")

idata <- read_repertoires(
  path = "<manifest>",
  manifest = manifest,
  schema = c("cdr3_aa", "v_call"),
  count_col = "counts",
  output_folder = "path/to/processed/cohort"
)
```

The default manifest file column is `file`. If it differs, pass the same name to `read_manifest(file_col = ...)` and `read_repertoires(manifest_file_col = ...)`. Relative paths are resolved against a manifest file; a manifest supplied as a data frame must contain absolute paths.

Use an explicit `repertoire_schema`, such as `c("donor", "timepoint")`, when those columns—not every manifest row—define the repertoires. Use `repertoire_schema = NULL` only when repertoires will be defined later.

### Single-cell receptors

Supplying `barcode_col` selects single-cell processing and requires `umi_col`; use it when observations must be connected through cells. It cannot be combined with `count_col`. Supply `locus_col` when selecting or pairing chains. When several chains from the same locus occur in a cell, the chain with the largest UMI/read value is retained.

```r
schema <- make_receptor_schema(
  features = c("junction_aa", "v_call", "j_call"),
  chains = c("TRA", "TRB")
)

idata <- read_repertoires(
  path = "path/to/contigs.tsv.gz",
  schema = schema,
  barcode_col = "cell_id",
  locus_col = "locus",
  umi_col = "umi_count",
  output_folder = "path/to/processed/paired"
)
```

Paired-chain mode retains cells containing both requested chains; use it when full αβ identity matters. For a single-chain analysis, set `chains = "TRB"`; use it when the selected chain alone defines the biological question or paired recovery is insufficient.

By default, ingestion standardises common 10x column names, removes selected technical columns, keeps productive sequences when possible, prefixes barcodes when the manifest has a `Prefix` column, defines repertoires, and writes the dataset. Override `rename_columns`, `preprocess`, `postprocess`, or `repertoire_schema` only when the input requires it.

## Transform and regroup

Use dplyr verbs on `ImmunData`; their methods retain duckplyr execution.

```r
fr_idata <- idata |>
  filter(Response == "FR")

length_idata <- idata |>
  mutate(cdr3_length = dd$length(cdr3_aa))
```

Sequence matching is configured with `make_seq_options()`:

```r
similar_idata <- idata |>
  filter(
    seq_options = make_seq_options(
      query_col = "cdr3_aa",
      patterns = "ASFPVLSPYNEQF",
      method = "lev",
      max_dist = 4
    )
  )
```

Choose the sequence method by question:

- `"exact"` finds identical sequences; use it for known clonotypes.
- `"regex"` finds motif patterns; use it for explicit sequence rules.
- `"lev"` allows substitutions, insertions, and deletions; use it for variable-length similarity.
- `"hamm"` counts substitutions at aligned positions; use it for equal-length sequences.
- In `filter()`, `max_dist` sets the acceptance threshold. In `mutate()`, distance methods annotate the distance instead of filtering.

Filtering recalculates existing repertoire summaries by default. To add cell annotations and then define repertoires from them, drop the previous summaries during annotation and aggregate afterward:

```r
cell_labels <- tibble::tibble(
  barcode = c("S1_1", "S1_2"),
  cell_type = c("CD8 T cell", "CD4 T cell")
)

cell_idata <- idata |>
  annotate_barcodes(
    annotations = cell_labels,
    annot_col = "barcode",
    keep_repertoires = FALSE
  ) |>
  filter(!is.na(cell_type)) |>
  agg_repertoires(schema = "cell_type")
```

Annotation helpers have distinct typical uses:

- `annotate_barcodes()` adds cell types, clusters, or other cell-level labels.
- `annotate_receptors()` adds specificity, cluster, or database-derived labels shared by all rows of a receptor.
- `annotate_chains()` adds contig- or chain-level quality information.
- `annotate()` joins clinical or experimental information through one or more custom keys using `c(immundata_column = "annotation_column")`.

Annotation keys must be unique; duplicate keys can duplicate AIRR rows and corrupt downstream counts. Set `keep_repertoires = FALSE` when the new annotation will define repertoires, then call `agg_repertoires()`.

## Repertoire analysis

Repertoire-level immunarch methods require defined repertoires. `get_test_idata()` provides a small documented example dataset.

```r
idata <- get_test_idata()
```

### Description

```r
chain_stats <- airr_desc_chains(idata)
length_stats <- airr_desc_lengths(idata, seq_col = "cdr3_aa")
gene_stats <- airr_desc_genes(
  idata,
  gene_col = "v_call",
  level = "receptor"
)

vis(chain_stats)
vis(length_stats)
vis(gene_stats)
```

- `airr_desc_chains()` counts chains per repertoire and locus; use it for library-size QC, chain balance, and locus dropout.
- `airr_desc_lengths()` returns sequence-length distributions; use it to detect technical shifts or compare repertoire structure.
- `airr_desc_genes()` measures gene-segment usage; use it for composition, enrichment, and between-group comparisons.

`airr_desc_genes(level = "receptor")` counts unique receptors and answers how many clonotypes use a gene. Use `level = "barcode"` to weight by cell/UMI abundance and answer how much repertoire mass uses it. The `by` argument controls additional grouping: the default uses the canonical locus when available, `NULL` pools loci, and a character vector selects explicit columns.

### Clonality

```r
vis(airr_clonality_line(idata, limit = 1000))
vis(airr_clonality_rank(idata, bins = c(10, 100, 1000)))
vis(airr_clonality_prop(idata))
```

- `airr_clonality_line()` orders receptors by abundance; use it to see whether a few receptors dominate and how quickly abundance decays.
- `airr_clonality_rank()` measures repertoire space occupied by top-N rank ranges; use it to compare concentration in the leading receptors.
- `airr_clonality_prop()` partitions receptors by within-repertoire proportion; use it to compare expansion classes across repertoires of different sizes.

### Diversity

```r
vis(airr_diversity_dxx(idata, perc = 50))
vis(airr_diversity_chao1(idata))
vis(airr_diversity_rarefaction(idata, step = 2000))
vis(airr_diversity_shannon(idata))
vis(airr_diversity_pielou(idata))
vis(airr_diversity_index(idata))

hill <- airr_diversity_hill(idata, q = c(0, 1, 2))
```

- `airr_diversity_dxx()` reports how many top receptors cover a chosen percentage; use it as an intuitive dominance measure. Smaller values mean stronger concentration.
- `airr_diversity_chao1()` estimates unseen richness from rare receptors; use it when sampling is incomplete and richness is the main question.
- `airr_diversity_rarefaction()` traces expected richness across sampling depth; use it to compare unequally sampled repertoires or assess sequencing saturation.
- `airr_diversity_shannon()` combines richness and evenness; use it as a general diversity summary, while remembering it remains richness-sensitive.
- `airr_diversity_pielou()` normalises Shannon entropy by observed richness; use it when the question is evenness rather than repertoire size.
- `airr_diversity_index()` returns the Hill number at `q = 1`; use it when diversity should be expressed as an effective number of equally abundant receptors.
- `airr_diversity_hill()` evaluates several `q` values; use it to see whether conclusions depend on rare receptors (`q = 0`) or increasingly abundant receptors (`q > 1`).

### Repertoire similarity

```r
vis(repsim_intersection(idata))
vis(repsim_jaccard(idata))
vis(repsim_chao_jaccard(idata))
vis(repsim_morisita_horn(idata))
vis(repsim_bray(idata))
```

- `repsim_intersection()` counts shared unique receptors; use it when the absolute number of shared sequences matters. It is sensitive to repertoire size.
- `repsim_jaccard()` compares presence/absence sets as shared over total unique receptors; use it for size-normalised sequence overlap when abundance should not matter.
- `repsim_chao_jaccard()` uses receptor counts to correct Jaccard similarity for incompletely sampled shared receptors; use it when limited depth may hide overlap.
- `repsim_morisita_horn()` compares abundance profiles, using both receptor identity and counts/proportions; use it when expanded shared receptors should contribute more than rare ones.
- `repsim_bray()` compares abundance profiles as dissimilarity; use it when the magnitude of compositional difference is the target. `0` means identical composition and larger values mean less similarity.

## Receptor annotations and single-cell transfer

Repertoire analysis functions return small summary tables or matrices. Annotation functions instead return a new `ImmunData` with receptor-level columns.

```r
idata <- annotate_clonality_prop(idata)
sdata <- annotate_seurat(idata, sdata, cols = "clonal_prop_bin")
Seurat::DimPlot(
  sdata,
  reduction = "umap",
  group.by = "clonal_prop_bin",
  shuffle = TRUE
)
```

- `annotate_clonality_prop()` adds expansion classes based on within-repertoire proportions; use it for cell-level plots of clonal expansion.
- `annotate_clonality_rank()` adds top-N rank bins; use it when relative ordering is easier to interpret than fixed proportions.
- `annotate_public()` adds receptor-level sharing counts and abundance summaries across repertoires; use it to identify and characterise shared receptors after repertoire aggregation.
- `annotate_seurat()` transfers selected barcode-level columns; use it to display repertoire annotations on single-cell embeddings.

## Persist and resume

`read_repertoires()` already writes the ingested dataset. After expensive transformations, save the returned state:

```r
saved_idata <- write_immundata(idata, "path/to/saved-state")
continued_idata <- read_immundata("path/to/saved-state")
```

For versioned snapshots, establish a project home with `write_immundata(..., output_folder = project_dir, rehome = TRUE)`, then call `write_immundata(idata, tag = "analysis-step")`. Read the latest tagged version with `read_immundata(project_dir, tag = "analysis-step")`, or pass `version` for an exact version. Reusing an explicit output folder replaces its existing ImmunData files.

- Use an explicit output folder for a named standalone state or a project home.
- Use tagged managed snapshots for iterative checkpoints where earlier versions must remain available.
- Use `read_immundata()` to resume without repeating ingestion or expensive transformations.
