# mpCGC

**marine polysaccharide CAZyme gene clusterer**

A Nextflow pipeline that identifies, mines and summarises carbohydrate-active
enzyme gene clusters (CGCs) in marine microbial genomes and metagenome-assembled
genomes. mpCGC extends [dbCAN](https://github.com/bcb-unl/run_dbcan) — which
supplies the gene calling, the CAZyme/sulfatase/transporter/peptidase searches
and the cluster calling — with a cross-genome catalog, an enzyme mining layer and
a CGC diversity layer. It is the pipeline behind [mpCGCdb](https://mpcgcdb.com),
an atlas of 289,962 CGCs across 22,607 marine MAGs.

```
                    ┌──────────────────────────────────────────────┐
   genome fasta ───►│  identify   gene calls ─► CAZyme, sulfatase,  │──► CGC catalog
                    │             transporter, peptidase and TF     │    (csv)
                    │             annotation ─► CGC determination   │
                    └───────────────────┬──────────────────────────┘
                                        │  CGC catalog
                      ┌─────────────────┴────────────────┐
                      ▼                                  ▼
      ┌───────────────────────────────┐   ┌──────────────────────────────────┐
      │  mine                         │   │  summarize                       │
      │  enzyme colocalization        │   │  cluster encoding                │
      │  all-vs-all DIAMOND (SSN)     │   │  hierarchical clustering         │
      │  Leiden communities           │   │  CGC diversity (iNEXT / Chao2)   │
      │  MUSCLE + FastTree + ggtree   │   │                                  │
      └───────────────────────────────┘   └──────────────────────────────────┘
              nwk trees, json networks           csv summary metrics
```

---

# Quick start

```bash
# 1. fetch the reference databases (~8 GB, one time)
nextflow run AaronAOliver/mpCGC --step download_db --outdir refs -profile conda

# 2. run the pipeline
nextflow run AaronAOliver/mpCGC \
    --input samplesheet.csv \
    --db refs/db/dbcan_db \
    --outdir results \
    -profile conda \
    -resume
```

---

# Inputs

There is exactly one required input file — a **samplesheet** — plus an optional
**metadata table**. Everything else is a parameter.

## The samplesheet

A CSV with a header row, one genome per row.

| Column | Required | Meaning |
|---|---|---|
| `sample` | **yes** | Name for this genome. Must be unique; it names every output for that genome and prefixes its CGC identifiers. |
| `fasta` | **yes** | The sequence file. What it must contain depends on `--mode` — see below. May be gzipped. |
| `gff` | only for `--mode protein` | Gene coordinates matching the proteins in `fasta`. |
| `taxonomy` | no | GTDB taxonomy string. Drives lineage grouping in the diversity estimates and the colouring of tree figures. |

Paths may be absolute or relative to where you launch Nextflow. Every path is
checked before the run starts, so a typo fails immediately rather than an hour in.

### Which sequence file? `--mode` decides

mpCGC accepts three kinds of input, and this is the choice people most often get
wrong.

#### `--mode prok` (default) — nucleotide assembly, one organism

`fasta` is a **nucleotide** assembly: a finished genome, a draft, or a MAG. Genes
are called for you with Pyrodigal.

```csv
sample,fasta,taxonomy
MAG_001,/data/MAG_001.fna,d__Bacteria;p__Bacteroidota;c__Bacteroidia
MAG_002,/data/MAG_002.fna.gz,d__Bacteria;p__Verrucomicrobiota
```

This is the right mode for MAGs and isolate genomes, and it is what the mpCGCdb
validation run used.

#### `--mode meta` — nucleotide assembly, mixed community

`fasta` is a **nucleotide** metagenome assembly that has not been binned. Same
columns as `prok`; only the gene-calling model differs, using settings suited to
fragmented contigs from many organisms.

Use this for raw assemblies. If you have binned your assembly into MAGs, run those
as separate `prok` rows instead — you will get per-genome CGC counts and usable
diversity estimates, neither of which an unbinned assembly can give you.

#### `--mode protein` — proteins **plus** a GFF

`fasta` is an **amino-acid** FASTA and the `gff` column is **required**.

```csv
sample,fasta,gff,taxonomy
MAG_001,/data/MAG_001.faa,/data/MAG_001.gff,d__Bacteria;p__Bacteroidota
```

Why the GFF is not optional: a CGC is defined by genes sitting next to each other
on a contig. A protein FASTA carries no coordinates, so without a GFF the pipeline
can annotate each protein but cannot tell which proteins are neighbours, and no
cluster can be called. The GFF supplies contig, start, stop and strand.

Requirements for the pair:

- Protein identifiers in the FASTA header (first whitespace-separated token) must
  match the identifiers in the GFF attributes.
- The GFF should be Prodigal-style. Set `--gff_type NCBI_prok` if it came from an
  NCBI annotation instead.
- Both must describe the same gene set. Proteins missing from the GFF cannot be
  placed and are silently dropped from clustering.

Use this mode when you already have an annotation you want to keep — for example
proteins from a curated genome — and want mpCGC's clusters on top of it. If you
have the nucleotide assembly, prefer `prok`: the gene calls then come from the
same Pyrodigal version as everything else, which is one fewer thing to reconcile.

### One genome or many?

The samplesheet scales from one row to tens of thousands, but what you can learn
depends on how many genomes you give it, and this matters most for the diversity
layer.

| Rows | What works |
|---|---|
| 1 | `identify` fully. `mine` works but every sequence-similarity network is drawn from one genome. `summarize` produces fingerprints and clustering, but **no rarefaction** |
| 2–4 per lineage | as above; diversity still skipped |
| ≥ 5 per lineage | rarefaction, extrapolation and Chao2 asymptotes are produced for that lineage |
| hundreds per lineage | the regime the estimates are designed for |

**The rarefaction step needs a suite of genomes, not one.** It treats genomes as
sampling units and distinct CGC fingerprints as species, then asks how richness
accumulates as genomes are added. With one genome there is one sampling unit and
no accumulation curve to fit, so `CGC_RICHNESS` skips that lineage with a message
rather than printing a meaningless asymptote:

```
[richness] skipping Bacteroidota (1 genomes < --min-genomes)
```

The threshold is `--min_genomes` (default 5). Five is a floor, not a
recommendation: Chao2 leans on how many fingerprints are seen in exactly one and
exactly two genomes, so with a handful of genomes the interval will be enormous
and honest rather than small and wrong. Lineages are formed by `--group_by`
(default `phylum`, using the mpCGC lineage groups), so "five genomes" means five
*within a lineage*, not five overall.

If you are running a single genome, you can skip the diversity work entirely:

```bash
nextflow run . --input one_genome.csv --db refs/db/dbcan_db --step identify
```

## The metadata table — taxonomy for grouping and colouring

Taxonomy can be supplied two ways. Both feed the same two consumers: lineage
grouping in the diversity estimates, and colouring of the ggtree figures.

### Option A — the `taxonomy` column in the samplesheet

Simplest, and enough for most runs:

```csv
sample,fasta,taxonomy
MAG_001,/data/MAG_001.fna,d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales
MAG_002,/data/MAG_002.fna,d__Archaea;p__Thermoproteota;c__Nitrososphaeria
```

### Option B — a separate metadata file via `--metadata`

Better for large cohorts, and the right choice if you already have a GTDB-Tk
summary or a MAG statistics table:

```bash
nextflow run . --input samplesheet.csv --metadata mag_metadata.tsv --db refs/db/dbcan_db
```

A **tab-separated** file with a header. Two columns are used:

| Column | Accepted names | Meaning |
|---|---|---|
| genome id | `Bin_id`, `sample`, `MAG` or `genome` | must match the `sample` column of the samplesheet |
| taxonomy | any header containing `taxonomy` (case-insensitive) | GTDB string |

```
Bin_id	Completeness (%)	Contamination (%)	GTDB_taxonomy (v207)
MAG_001	98.4	0.7	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Flavobacteriaceae;g__Polaribacter
MAG_002	95.1	1.2	d__Archaea;p__Thermoproteota;c__Nitrososphaeria
```

Extra columns are ignored, so a GTDB-Tk `gtdbtk.bac120.summary.tsv` works as-is,
as does an existing MAG metadata sheet. `--metadata` overrides the samplesheet's
`taxonomy` column when both are present. Genomes absent from the file, or with an
unparseable string, fall into the `other` group rather than failing the run.

### How taxonomy reaches the tree figures

Tree tips are protein identifiers, which carry no taxonomy, so the colouring is
resolved by two joins the pipeline performs for you:

```
protein  ──(catalog/proteins/protein_map.tsv, written by CGC_PROTEINS)──►  genome
genome   ──(--metadata, or the samplesheet taxonomy column)────────────►  lineage
```

Turn the figures on and pick what to colour by:

```bash
nextflow run . \
    --input samplesheet.csv \
    --metadata mag_metadata.tsv \
    --db refs/db/dbcan_db \
    --make_tree_figures true \
    --tree_color_by taxonomy \
    --tree_color_rank phylum
```

| Parameter | Default | Meaning |
|---|---|---|
| `--make_tree_figures` | `false` | render the trees at all; needs the `rggtree` environment |
| `--tree_color_by` | `taxonomy` | `taxonomy` or `community` (the Leiden community from the sequence similarity network) |
| `--tree_color_rank` | `phylum` | `phylum`, `class`, `order`, `family`, `genus` |

At `--tree_color_rank phylum` the colours are the 11-group mpCGC palette — the
eight dominant phyla, Alpha- and Gammaproteobacteria split out because the phylum
is too coarse to be informative, Archaea pooled, everything else grey — so tree
figures match the other figures in the manuscript. At finer ranks the observed
values are coloured in order of abundance, with everything past
`--max-groups` pooled into grey.

If taxonomy is unavailable the script says so and falls back to colouring by
Leiden community rather than failing:

```
[ggtree] no protein map or taxonomy table given; falling back to community colouring
```

Trees are written for every family with at least `--min_family_seqs` sequences
(default 4). Restrict the work with `--mine_families GH16_16,PL7_5` while you are
finding your feet — mining every family in a large catalog is the most expensive
part of the pipeline.

---

# The three dataflows

Each can be run alone with `--step`; all three run by default.

## `--step identify`

Genes are called with Pyrodigal, then every protein is searched against CAZy
(DIAMOND), the dbCAN and dbCAN-sub HMMs (PyHMMER), SulfAtlas for sulfatases, TCDB
for transporters including SusC/SusD pairs, MEROPS for peptidases, PRODORIC for
transcription factors, and the signal-transduction HMMs. CGCs are then called by a
seed-and-extend procedure: clusters are seeded on co-localised signature genes,
extended across at most `--num_null_gene` unannotated genes, and widened by a
user-provided gene-boundary extension limit (`--extend_mode`).

Output: one annotated directory per genome, plus `catalog/all_cgc_catalog.tsv`.

## `--step mine`

The catalog becomes an enzyme colocalisation network, where families are nodes and
edges are weighted by the number of CGCs containing both. Independently, the
proteins of each enzyme family are compared all-vs-all with DIAMOND
(`-k 0 -e 1e-30`) to build a sequence similarity network, partitioned into
communities with the Leiden algorithm at `--leiden_resolution`. Families are then
aligned with MUSCLE v5 and a maximum-likelihood phylogeny inferred with FastTree,
optionally rendered with ggtree.

Output: `mine/colocalization/*.json`, `mine/communities/*.tsv`, `mine/trees/*.nwk`.

## `--step summarize`

Every CGC is encoded as a binary fingerprint over S1, GH, PL and CE families (one
token per gene, priority GH > PL > CE > S1). Fingerprints are clustered on
pairwise Hamming distance, and two clusters are treated as overlapping in function
when one contains a superset of all annotations of the other. CGC diversity per
lineage is estimated by incidence-based rarefaction and extrapolation in the iNEXT
framework, with Chao2 asymptotic richness.

By default `--collapse_subsets` folds each fingerprint into its maximal superset
before counting, so a cluster nested inside a larger one is not counted as a
separate function. Given the four clusters `GH12`, `(GH12, S1_7)`, `GH23` and
`S1_7`, the collapsed count is two: `(GH12, S1_7)` and `GH23`.

Output: `summarize/cgc_fingerprints.tsv`, `summarize/richness/*`,
`summarize/mpcgc_summary.csv`.

Running the steps separately is fine — `mine` and `summarize` read the catalog
from `<outdir>/catalog/`, so point `--outdir` at the completed identification run:

```bash
nextflow run . --input samplesheet.csv --db refs/db/dbcan_db --step identify  --outdir results
nextflow run . --input samplesheet.csv --db refs/db/dbcan_db --step summarize --outdir results
```

---

# Profiles

| Profile | Use |
|---|---|
| `conda` / `mamba` | per-process conda environments from `env/*.yml` |
| `docker` / `singularity` | one image for every process; build it with `docker build -t mpcgc:1.0.0 .` and select it with `--container` |
| `slurm` | submit to a Slurm cluster |
| `local_envs` | reuse conda environments already on the host (`--local_env_dbcan`, `--local_env_pydata`, `--local_env_network`, `--local_env_phylo`) |
| `test` | one genome, identification only |
| `test_full` | two genomes, all three dataflows |

# Key parameters

```
--step                  identify | mine | summarize | all | download_db
--mode                  prok | meta | protein

CGC calling
--min_core_cazyme       minimum CAZymes per cluster            [1]
--num_null_gene         max intervening unannotated genes      [2]
--additional_genes      accessory classes required with the CAZyme core [TC]
--additional_logic      all | any                              [all]
--extend_mode           none | bp | gene                       [gene]
--extend_gene_count     flanking genes when --extend_mode gene [2]

Mining
--mine_families         restrict to a comma-separated subset   [all]
--min_family_seqs       skip families smaller than this        [4]
--ssn_evalue            DIAMOND all-vs-all E-value             [1e-30]
--leiden_resolution     Leiden resolution                      [1.0]
--make_tree_figures     render trees with ggtree               [false]
--tree_color_by         taxonomy | community                   [taxonomy]
--tree_color_rank       phylum | class | order | family | genus [phylum]

Summary
--collapse_subsets      count maximal fingerprints only        [true]
--group_by              metadata column defining lineages      [phylum]
--min_genomes           lineages smaller than this skip rarefaction [5]
--cooccurrence_min_cgc  min CGCs for a family to be plotted    [100]
```

`nextflow run AaronAOliver/mpCGC --help` prints the full list.
`docs/usage.md` documents every parameter; `docs/output.md` describes every file.

# Requirements

- Nextflow 24.04 or newer, and Java 17+
- conda/mamba, Docker or Singularity
- ~8 GB for the reference databases, and ~32 GB RAM for the dbCAN-sub HMM step
  (lower it with `--dbcansub false`, at the cost of subfamily resolution)

# Tests

```bash
nextflow lint .              # scripts and configs
python test/test_inext.py    # diversity estimators against iNEXT 3.0.2 in R
```

`test/test_inext.py` needs only numpy and scipy: the expected values were
generated once with `test/chao_ref.R` and are stored in the test. All 21 compared
quantities agree with R to within 0.0005.

To check a run against mpCGCdb, `test/reference/` ships the published catalog rows
for both validation genomes, so this works straight from a clone:

```bash
nextflow run . -profile test,conda --db refs/db/dbcan_db --outdir results
python test/validate_against_mpcgcdb.py \
    --run results/catalog/all_cgc_catalog.tsv \
    --reference test/reference/mpcgcdb_GCA_000325705.1_ASM32570v1_genomic.tsv.gz \
    --genome GCA_000325705.1_ASM32570v1_genomic
```

# Reproducing mpCGCdb

The defaults reproduce the settings used to build mpCGCdb, including the
gene-boundary extension of two genes and the requirement that each cluster carry a
transporter. `docs/validation.md` records per-gene comparisons against the
published records for two genomes chosen to bracket the range:

| | *Echinicola vietnamensis*<br>`GCA_000325705.1` | *Thermococcus* sp. AM4<br>`GCA_000151205.2` |
|---|---|---|
| | Bacteroidota, 5.6 Mb | Archaea, 2.1 Mb |
| Clusters | 94 vs 82 published | **16 vs 16** |
| Reference genes matched | **915 / 915 (100 %)** | **163 / 163 (100 %)** |
| Genes missing from the run | **0** | **0** |
| Identical enzyme families | 97.7 % | 96.3 % |
| Clusters identical gene for gene | 90.2 % | 87.5 % |
| Family inventory recovered | 99.3 % | **100 %** |

Every remaining difference is an addition from reference databases that have grown
since — new CAZy families, extra CBM modules on genes mpCGCdb already had, one
reclassification. Nothing present in mpCGCdb is missing from either run. That
document also records one open discrepancy about the transporter requirement, and
what the tests do and do not cover.

# Related

- [mpCGCdb](https://mpcgcdb.com) — the database this pipeline produced
- [KiritimatiallesSCoNe](https://github.com/AaronAOliver/KiritimatiallesSCoNe) —
  sequence similarity / colocalisation network construction for the
  Kiritimatiellales

# Citing

See `CITATION.cff`. mpCGC builds on [dbCAN](https://github.com/bcb-unl/run_dbcan)
for annotation and CGC calling; please cite dbCAN, CAZy, SulfAtlas, TCDB, MEROPS,
DIAMOND, HMMER, Pyrodigal, MUSCLE, FastTree, iNEXT and Leiden as appropriate. A
`versions.yml` recording the exact version of every tool is written for each run.

# Licence

MIT (code). The mpCGCdb dataset is CC-BY 4.0.
