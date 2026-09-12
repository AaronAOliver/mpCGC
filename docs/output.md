# Output

```
results/
├── identify/
│   ├── <sample>/                        one directory per genome
│   │   ├── uniInput.faa                 Pyrodigal proteins
│   │   ├── uniInput.gff                 Pyrodigal gene calls
│   │   ├── diamond.out                  CAZy hits
│   │   ├── diamond.out.sulfatase        SulfAtlas hits
│   │   ├── diamond.out.tc               TCDB hits
│   │   ├── diamond.out.tf               PRODORIC hits
│   │   ├── diamond.out.peptidase        MEROPS hits
│   │   ├── dbCAN_hmm_results.tsv        dbCAN family HMM hits
│   │   ├── dbCANsub_hmm_results.tsv     dbCAN-sub subfamily HMM hits
│   │   ├── STP_hmm_results.tsv          signal-transduction hits
│   │   └── overview.tsv                 consensus per-protein annotation
│   ├── <sample>.faa                     copy of uniInput.faa, named per sample
│   ├── <sample>.overview.tsv            copy of overview.tsv, named per sample
│   └── cgc/
│       ├── <sample>.cgc.gff             annotated GFF used for CGC calling
│       ├── <sample>.cgc_standard_out.tsv   CGCs for that genome
│       └── <sample>.total_cgc_info.tsv
├── catalog/
│   ├── all_cgc_catalog.tsv              every CGC gene across every genome
│   ├── cgc_per_genome.tsv               per-genome counts by gene type
│   └── proteins/
│       ├── <family>.faa                 CGC member proteins, one file per family
│       └── family_counts.tsv
├── mine/
│   ├── colocalization/
│   │   ├── colocalization_edges.tsv     family pairs with CGC and genome support
│   │   └── colocalization_network.json  Cytoscape graph
│   ├── ssn/<family>.ssn.tsv.gz          all-vs-all DIAMOND output
│   ├── communities/
│   │   ├── <family>.communities.tsv     protein to Leiden community
│   │   └── <family>.network.json        Cytoscape graph of the SSN
│   ├── alignments/<family>.aln.faa      MUSCLE alignment
│   ├── trees/<family>.nwk               FastTree phylogeny
│   └── tree_figures/<family>.tree.png   ggtree rendering (optional)
├── summarize/
│   ├── cgc_fingerprints.tsv             per-CGC enzyme fingerprint
│   ├── fingerprint_collapse.tsv         subset to maximal-superset map
│   ├── clustering/
│   │   ├── cgc_clusters.tsv             flat clusters + containment counts
│   │   ├── cgc_linkage.npy              linkage matrix
│   │   └── cgc_dendrogram.png
│   ├── cooccurrence/
│   │   ├── cooccurrence_matrix.tsv      row-conditioned colocalisation matrix
│   │   ├── cooccurrence_order.txt       clustered family order
│   │   └── cooccurrence_heatmap.png
│   ├── richness/
│   │   ├── rarefaction_curves.tsv       interpolated and extrapolated curves
│   │   ├── richness_asymptotes.tsv      Chao2 asymptotes with intervals
│   │   └── rarefaction.png
│   ├── mpcgc_summary.csv                per-genome metrics
│   └── mpcgc_summary_lineage.csv        per-lineage metrics
└── pipeline_info/
    ├── report.html, timeline.html, trace.txt, dag.html
    └── versions.yml                     exact version of every tool used
```

## Key files

### `catalog/all_cgc_catalog.tsv`

One row per gene inside a CGC. This is the table every later step reads, and it
matches the schema published on mpcgcdb.com.

| Column | Notes |
|---|---|
| `CGC#` | `<sample>_CGC<n>`, unique across the run |
| `MAG` | genome name from the samplesheet |
| `Gene Type` | `CAZyme`, `TC`, `TF`, `STP`, `Sulfatase`, `Peptidase` or `null` |
| `Contig ID`, `Protein ID` | |
| `Gene Start`, `Gene Stop`, `Gene Strand` | |
| `Gene Annotation` | e.g. `CAZyme\|GH13_3+CBM48`, `Sulfatase\|S1_7`, `TC\|1.B.14.6.1` |

The `<sample>.faa` and `<sample>.overview.tsv` copies carry the same content as
the files inside `<sample>/`. They exist because every genome's dbCAN run writes
`uniInput.faa`, and collecting many genomes into one process needs distinct
filenames.

`null` rows are the unannotated genes that fall inside a cluster. They are part of
the locus and are deliberately kept: a CGC is a stretch of DNA, not just the
genes that happened to match a database.

`Gene Type` is one label per gene even when a gene matches several databases; the
full set of matches stays in `Gene Annotation`.

### `summarize/cgc_fingerprints.tsv`

| Column | Notes |
|---|---|
| `CGC`, `MAG`, `lineage` | |
| `n_tokens` | families in the fingerprint |
| `fingerprint` | the CGC's own sorted token set |
| `counted_as` | the maximal superset it is counted as; equals `fingerprint` when the CGC is itself maximal |

Count distinct functions with `counted_as`, not `fingerprint`. With
`--collapse_subsets false` the two columns are always identical.

### `summarize/richness/richness_asymptotes.tsv`

| Column | Notes |
|---|---|
| `genomes` | sampling units (T) |
| `S_obs` | fingerprints observed |
| `Q1`, `Q2` | fingerprints seen in exactly one and exactly two genomes |
| `Q0_chao2` | estimated undetected fingerprints |
| `S_chao2` | `S_obs + Q0_chao2` |
| `SE` | analytic standard error of `S_chao2` |
| `ci_lower`, `ci_upper` | log-transformed 95% interval, asymmetric about the estimate |
| `pct_recovered` | `S_obs / S_chao2` |

The interval is asymmetric by construction and never falls below `S_obs`, which is
why it is not reported as `S_chao2 ± 1.96·SE`.

### `summarize/cooccurrence/cooccurrence_matrix.tsv`

Row-conditioned, so the matrix is **not** symmetric:

```
M[anchor, colocalized] = P(colocalized family present | anchor family present)
```

Read along a row to ask what a given family travels with. A rare family that
always accompanies a common one scores high in its own row and low in the common
family's row — the asymmetry is the useful part, and averaging it away loses the
direction of the association.

### `pipeline_info/versions.yml`

Collected per process, recording dbCAN, Pyrodigal, DIAMOND, PyHMMER, MUSCLE,
FastTree, igraph, scipy and numpy versions for the run. Include it with any
manuscript built on a run.
