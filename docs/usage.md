# Usage

## Installing

```bash
curl -s https://get.nextflow.io | bash      # or: conda install -c bioconda nextflow
nextflow pull AaronAOliver/mpCGC
```

Nextflow 24.04+ and Java 17+ are required. Everything else is provided by the
`conda`, `docker` or `singularity` profile.

## Databases

```bash
nextflow run AaronAOliver/mpCGC --step download_db --outdir refs -profile conda
```

This calls dbCAN's own downloader, which fetches, in one step:

| Database | File | Purpose |
|---|---|---|
| CAZy | `CAZy.dmnd` (~2.2 GB) | CAZyme assignment by DIAMOND |
| dbCAN | `dbCAN.hmm` (~130 MB) | CAZyme family HMMs |
| dbCAN-sub | `dbCAN-sub.hmm` (~5 GB) | CAZyme **subfamily** HMMs |
| SulfAtlas | `sulfatlas_db.dmnd` (~100 MB) | sulfatase subfamilies (S1_*) |
| TCDB | `TCDB.dmnd` (~14 MB) | transporters, SusC/SusD |
| MEROPS | `peptidase_db.dmnd` (~330 MB) | peptidases |
| PRODORIC | `TF.dmnd` | prokaryotic transcription factors |
| STP | `STP.hmm` (~11 MB) | signal-transduction proteins |
| dbCAN-PUL | `dbCAN-PUL/` | reference polysaccharide utilisation loci |

Budget ~8 GB of disk. The result is a single directory to pass as `--db`. Add
`--db_from_s3` to pull the pinned S3 release rather than the moving `db_current`
snapshot — worth doing if you need a run to be reproducible months later.

Databases only need downloading once; point every later run at the same `--db`.

## Input

A CSV samplesheet with a header:

```csv
sample,fasta,taxonomy
MAG_001,/data/MAG_001.fna,d__Bacteria;p__Bacteroidota;c__Bacteroidia
MAG_002,/data/MAG_002.fna.gz,d__Bacteria;p__Verrucomicrobiota
```

- `sample` — used to name every output for that genome. Must be unique.
- `fasta` — nucleotide assembly, optionally gzipped. With `--mode protein` this
  is a protein fasta and a matching GFF is required.
- `taxonomy` — optional GTDB string. Used only to assign genomes to lineages for
  the diversity estimates; genomes without one fall into `other`.

For a larger cohort, pass `--metadata` with a TSV carrying a `Bin_id` column and
a column whose name contains `taxonomy`; extra columns such as completeness are
carried through. The samplesheet is still required.

## Running

```bash
nextflow run AaronAOliver/mpCGC \
    --input samplesheet.csv \
    --db refs/db/dbcan_db \
    --outdir results \
    -profile conda \
    -resume
```

Always keep `-resume` in your command: an interrupted run restarts from the last
completed process rather than re-annotating everything.

### One dataflow at a time

```bash
# annotate and call clusters only
nextflow run . --input samplesheet.csv --db refs/db/dbcan_db --step identify

# later, mine and summarise the catalog that run produced
nextflow run . --input samplesheet.csv --db refs/db/dbcan_db --step mine      --outdir results
nextflow run . --input samplesheet.csv --db refs/db/dbcan_db --step summarize --outdir results
```

`--step mine` and `--step summarize` read `<outdir>/catalog/all_cgc_catalog.tsv`,
so point `--outdir` at the directory the identification run wrote.

### On a cluster

```bash
nextflow run AaronAOliver/mpCGC \
    --input samplesheet.csv --db /shared/refs/dbcan_db \
    -profile singularity,slurm \
    --max_cpus 32 --max_memory 128.GB
```

Edit the `slurm` profile in `nextflow.config` to match your queue names. Resource
requests are attached per label in `conf/base.config`; `resourceLimits` caps every
request (including retries, which double memory) at the `--max_*` ceiling.

### Memory

The dbCAN-sub HMM is about 5 GB and PyHMMER preloads it, so the annotation step
carries the `bigmem` label and asks for 32 GB. If that is more than you have:

```bash
--dbcansub false          # skip subfamily HMMs entirely (families only)
--max_memory 24.GB        # or just lower the ceiling; PyHMMER batches adaptively
```

Dropping dbCAN-sub costs subfamily resolution — `GH13` instead of `GH13_3` — which
matters, because the fingerprints and the colocalisation analysis are built on
subfamilies.

## Parameters

### Identification

| Parameter | Default | Meaning |
|---|---|---|
| `--mode` | `prok` | `prok`, `meta` or `protein` |
| `--methods` | `diamond,hmm,dbCANsub` | CAZyme annotation modules |
| `--e_value_cazyme` | `1e-102` | DIAMOND E-value against CAZy |
| `--dbcansub` | `true` | run the subfamily HMMs |
| `--e_value_tc` / `--coverage_tc` | `1e-4` / `35` | TCDB search |
| `--e_value_sulfatase` / `--coverage_sulfatase` | `1e-4` / `35` | SulfAtlas search |
| `--e_value_peptidase` / `--coverage_peptidase` | `1e-4` / `35` | MEROPS search |
| `--e_value_stp` / `--coverage_stp` | `1e-4` / `0.35` | STP HMMs |

### CGC calling

| Parameter | Default | Meaning |
|---|---|---|
| `--min_cluster_genes` | `2` | minimum genes in a cluster |
| `--min_core_cazyme` | `1` | minimum CAZymes in a cluster |
| `--num_null_gene` | `2` | consecutive unannotated genes tolerated inside a cluster |
| `--additional_genes` | `TC` | accessory classes required alongside the CAZyme core |
| `--additional_logic` | `all` | `all` = every listed class required; `any` = at least `--additional_min_categories` |
| `--extend_mode` | `gene` | `none`, `bp` or `gene` — the gene-boundary extension limit |
| `--extend_gene_count` | `2` | flanking genes added each side when `--extend_mode gene` |
| `--extend_bp` | `0` | flanking base pairs when `--extend_mode bp` |
| `--use_distance` | `false` | also require signature genes within `--base_pair_distance` |

The default requires a transporter in every cluster, which is the classic
polysaccharide-utilisation-locus definition and the setting that best reproduces
mpCGCdb. To also keep clusters built from CAZymes alone, add `CAZyme` to the list
and switch to `any` logic:

```bash
--additional_genes CAZyme,TC --additional_logic any
```

A cluster then qualifies on its CAZyme core alone, so a pair of adjacent
glycoside hydrolases with no transporter is still a CGC. This finds more loci but
widens the boundaries of some existing ones, so the two settings are not directly
comparable — see the sweep in `docs/validation.md`.

### Mining

| Parameter | Default | Meaning |
|---|---|---|
| `--mine_families` | all | comma-separated subset, e.g. `GH16_16,PL7_5` |
| `--min_family_seqs` | `4` | families smaller than this are skipped |
| `--ssn_evalue` | `1e-30` | DIAMOND all-vs-all cutoff |
| `--leiden_resolution` | `1.0` | higher splits communities more finely |
| `--muscle_super5_min` | `1000` | sequence count at which MUSCLE switches to `--super5` |
| `--make_trees` | `true` | MUSCLE + FastTree |
| `--make_tree_figures` | `false` | additionally render trees with ggtree |

Mining every family in a large catalog is the most expensive part of the
pipeline. Start with `--mine_families` on the families you care about.

### Summary

| Parameter | Default | Meaning |
|---|---|---|
| `--collapse_subsets` | `true` | count only maximal fingerprints |
| `--ambiguous_subset` | `most_frequent` | how to resolve a subset with several maximal supersets |
| `--group_by` | `phylum` | metadata column defining lineages |
| `--rarefaction_knots` | `40` | points on each curve |
| `--rarefaction_extrapolate` | `2.0` | extrapolate to this multiple of the observed genome count |
| `--rarefaction_nboot` | `0` | `0` uses the analytic Chao2 interval; >0 bootstraps the curve |
| `--cooccurrence_min_cgc` | `100` | families below this are left out of the matrix |

Diversity estimates need a reasonable number of genomes per lineage; lineages
with fewer than five are skipped rather than reported with a meaningless
asymptote.

## Troubleshooting

**`Process requirement exceeds available memory`** — lower `--max_memory` to just
under what the machine has. `resourceLimits` then caps every request.

**`conda: command not found`** — Nextflow needs `conda` on `PATH` to build
environments. Either put it there, or use `-profile local_envs` with
`--local_env_dbcan`, `--local_env_pydata`, `--local_env_network` and
`--local_env_phylo` pointing at environments you have already created.

**A family fails in `LEIDEN_COMMUNITIES`** — families whose all-vs-all search
produced no edge above `--ssn_evalue` are written out as singletons rather than
failing. If a run does fail there, check that `python-igraph` and `leidenalg` are
present in the `network` environment.

**Empty `diamond.out.peptidase`** — normal for genomes with no MEROPS hits above
the threshold; the file is created regardless.

**Tree figures appear as PNG but not SVG** — ggplot2 4 needs the `svglite` package
to write SVG. The conda environment installs it; the stock ggtree container does
not carry it, so the pipeline falls back to the cairo device and, failing that,
keeps the PNG and logs `SVG not written`. Nothing fails over a missing SVG.
