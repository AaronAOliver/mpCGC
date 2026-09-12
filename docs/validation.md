# Validation

Two checks are recorded here: a gene-by-gene comparison of a pipeline run against
the published mpCGCdb record for one genome, and a numerical check of the
diversity estimators against iNEXT in R.

## 1. Reproducing an mpCGCdb genome

**Genome:** `GCA_000325705.1_ASM32570v1_genomic` — *Echinicola vietnamensis*
(Bacteroidota; Cytophagales), a single complete 5.6 Mb chromosome, 100 %
complete / 0 % contamination in the mpCGCdb metadata. A marine, CAZyme-rich
Bacteroidota, so the comparison exercises glycoside hydrolases, polysaccharide
lyases, sulfatases and SusC/SusD transporters rather than a sparse genome.

**Run:**

```bash
nextflow run . -profile test,local_envs --db ~/dbcan_db --outdir results
python test/validate_against_mpcgcdb.py \
    --run results/catalog/all_cgc_catalog.tsv \
    --reference test/reference/mpcgcdb_GCA_000325705.1_ASM32570v1_genomic.tsv.gz \
    --genome GCA_000325705.1_ASM32570v1_genomic
```

Genes are matched on `(contig, start, stop)` rather than on protein identifier,
because identifiers are only unique within the run that produced them.

### Gene calling

| | Run | mpCGCdb |
|---|---|---|
| Reference genes matched on coordinates | **915 / 915 (100 %)** | |
| Genes only in the reference | **0** | |
| Genes only in the run | 124 | |

Pyrodigal reproduces the gene calls exactly: every gene mpCGCdb places in a CGC
is present at identical coordinates. The 124 extra genes belong to clusters that
only the newer databases find (see below).

### Per-gene annotation, over the 915 shared genes

| | Agreement |
|---|---|
| Identical `Gene Type` | 890 (97.3 %) |
| Identical enzyme families | **894 (97.7 %)** |
| Byte-identical `Gene Annotation` | 867 (94.8 %) |

`Gene Type` differences:

| mpCGCdb | Run | n | Cause |
|---|---|---|---|
| `prodoric` | `TF` | 8 | label change in dbCAN; same PRODORIC hit |
| `null` | `CAZyme` | 11 | newer CAZy/dbCAN-sub finds a family where the older release found none |
| `null` | `STP` | 3 | as above, signal-transduction HMMs |
| `STP` | `CAZyme` | 2 | gene now matches a CAZyme, which outranks STP in the type priority |
| `prodoric` | `STP` | 1 | as above |

All 21 enzyme-family differences are **additions** by the current databases — for
example `GH10` → `CBM4,GH10`, `CE1` → `CBM48,CE1`, and four genes newly assigned
`AA2`, `GH109`, `GT1` or `GT2`. No annotation present in mpCGCdb was lost.

### Cluster boundaries

| | Value |
|---|---|
| mpCGCdb clusters | 82 |
| Run clusters | 94 |
| Reference clusters whose gene set is **identical** | **74 (90.2 %)** |
| Median Jaccard against best match | **1.000** |
| Mean Jaccard against best match | 0.975 |
| Reference clusters with ≥ 0.90 overlap | 76 (92.7 %) |
| Reference clusters with ≥ 0.50 overlap | 81 (98.8 %) |
| Reference clusters split across two run clusters | 5 |

The 12 additional clusters in the run contain **no** mpCGCdb gene at all. Each is
seeded by a family that did not exist, or was not detected, in the CAZy release
mpCGCdb was built from: `GH177` (×3), `GH179`, `GH109`, `GH16`, `CBM50`, `CBM5`,
`AA1` (×2), `GT2`. `GH177` and `GH179` are recent CAZy additions. These are new
findings from a newer database, not disagreements about the same evidence.

### Enzyme family inventory

| | Value |
|---|---|
| Families in mpCGCdb | 137 |
| Families in the run | 147 |
| Families recovered | **136 / 137 (99.3 %)** |
| Shared families with identical gene counts | 127 / 136 |
| Only in mpCGCdb | `GT94` |
| Only in the run | `AA1`, `AA2`, `CBM4`, `CBM48`, `CBM50`, `GH16`, `GH177`, `GH179`, `GH93`, `GT1`, `S1_21` |

The single missing family is a reclassification, not a miss: the gene at
`CP003346.1:5075583-5076641` is `GT94` in mpCGCdb and `GT4` in the current CAZy
release — same gene, same coordinates, reassigned upstream.

### What this run established about the defaults

Three settings were fixed by this comparison, and are now the pipeline defaults.

**`--extend_mode gene --extend_gene_count 2`.** mpCGCdb clusters carry up to two
unannotated genes past their outermost signature gene on each side. With no
extension the run produced 669 gene rows against the reference's 915 and none of
the boundaries matched; with an extension of two genes, 90.2 % of clusters became
byte-identical. This is the manuscript's "user-provided gene-boundary extension
limit".

**`--additional_genes TC --additional_logic all`.** Four signature-class
configurations were run against the same annotated GFF, so only the cluster-calling
rule differed:

| `--additional_genes` | `--additional_logic` | Clusters | Identical to reference | Median Jaccard |
|---|---|---|---|---|
| `TC` (dbCAN default) | `all` | 94 | **74 (90.2 %)** | **1.000** |
| `CAZyme,TC` | `any` | 104 | 74 (90.2 %) | 1.000 |
| `CAZyme,TC,TF,STP,Sulfatase,Peptidase` | `any` | 106 | 13 (15.9 %) | 0.667 |
| `TC,TF,STP,Sulfatase,Peptidase` | `any` | 104 | 13 (15.9 %) | 0.667 |

Requiring a transporter wins clearly. The permissive `any` settings that include
the accessory classes do not merely add clusters — they redraw the boundaries of
clusters that would otherwise match, which is why their exact-match rate collapses
to 15.9 %. Adding `CAZyme` alone preserves every match but contributes ten
clusters the published catalog does not contain.

This corrects an inference drawn earlier from catalog-wide statistics. Across all
of mpCGCdb, 3.6 % of clusters contain no transporter and 5,539 consist only of
CAZymes, neither of which `TC`/`all` can produce — which had suggested a
permissive setting. The per-gene comparison shows the opposite for this genome:
all 82 of its published clusters contain a transporter, and the permissive
settings reproduce them badly. The likeliest explanation is that mpCGCdb was built
with an older dbCAN whose cluster validation differed in how it treated sulfatase,
peptidase and signal-transduction genes; that remains an open discrepancy rather
than something this comparison settles. `--additional_genes` and
`--additional_logic` expose the choice.

**Stripping the dbCAN-sub eCAMI cluster index.** The current dbCAN-sub writes
family names such as `GH140_e33`, `CE4_e749` and `AA1_e`; the trailing `_e<n>` is
eCAMI's internal cluster number, not a CAZy subfamily, and mpCGCdb's annotations
do not carry it. Removing it is also a correctness fix independent of
reproducibility: left in place it would scatter one family across dozens of
pseudo-families in the fingerprints, the colocalisation matrix and the per-family
phylogenies. Family agreement rose from 85.0 % to 97.7 % once it was removed.
Pass `--keep_ecami` to retain it.

### Summary

Gene calls are identical, 97.7 % of shared genes carry the same enzyme families,
99.3 % of the family inventory is recovered, and 90.2 % of clusters are reproduced
gene for gene with a median Jaccard of 1.000. Every remaining difference traces to
a reference database that has grown since mpCGCdb was built: new CAZy families,
extra CBM modules on known genes, and one reclassification. Nothing found in
mpCGCdb is missing from the run.

Reproducing the original catalog byte for byte would require pinning the exact
dbCAN release used at the time; `--db_from_s3` fetches a pinned snapshot rather
than the moving `db_current` and is the right choice when a run must stay
reproducible.

One caveat on scope: these numbers come from two genomes, not the whole catalog.
They were chosen to bracket the range — a CAZyme-rich marine Bacteroidota with 82
clusters across 137 families, and an archaeon a third the size with 16 clusters
across 18 families — and both reproduce with the same defaults. A
genome-to-genome comparison across all 22,607 would still be needed to claim the
defaults are right everywhere, and the transporter-requirement discrepancy above
is direct evidence that they are not right for every cluster in mpCGCdb.

## 1b. A second genome, from the other domain

**Genome:** `GCA_000151205.2_ASM15120v2_genomic` — *Thermococcus* sp. AM4
(Archaea; Thermococci), a single complete 2.09 Mb chromosome. Chosen to be as
unlike the first genome as the catalog allows: a different domain, a third of the
size, an eighth of the CGCs, and a mostly glycogen/trehalose-oriented enzyme set
rather than a marine-polysaccharide one.

| | Run | mpCGCdb |
|---|---|---|
| Clusters | **16** | **16** |
| Reference genes matched on coordinates | **163 / 163 (100 %)** | |
| Genes only in the reference | **0** | |
| Identical `Gene Type` | 157 (96.3 %) | |
| Identical enzyme families | 157 (96.3 %) | |
| Reference clusters with identical gene sets | **14 / 16 (87.5 %)** | |
| Median Jaccard against best match | **1.000** | |
| Family inventory recovered | **18 / 18 (100 %)** | |

The cluster count matches exactly. The six per-gene differences are the same
pattern as the first genome and all additive: three genes the older release called
`null` now carry a CAZyme, two now carry a signal-transduction domain, and one
transporter now also matches a CAZyme. The four families present only in the run
— `CBM20`, `CBM34`, `CBM48`, `CE1` — are carbohydrate-binding modules and an
esterase detected by the current dbCAN-sub on genes mpCGCdb already had.

That a genome this different reproduces with the same defaults, and the same kind
of residual differences, is the main reason to trust that the settings are not
overfitted to one Bacteroidota.

## 2. Diversity estimators against R

`test/test_inext.py` checks `bin/mpcgc_inext.py` against
`iNEXT::ChaoRichness(datatype = "incidence_freq")` from iNEXT 3.0.2 under R 4.5.2,
on three synthetic incidence datasets spanning low to high sample coverage. The
expected values are stored in the test, so R is not needed to run it;
`test/chao_ref.R` regenerates them.

```
python test/test_inext.py
```

| Dataset | T | S_obs | Q1 | Q2 | S_chao2 (mpCGC / R) | SE | 95 % CI |
|---|---|---|---|---|---|---|---|
| low | 40 | 54 | 30 | 12 | 90.5625 / 90.5620 | 18.061 | 68.63 – 145.37 |
| mid | 120 | 124 | 45 | 30 | 157.4688 / 157.4690 | 13.053 | 140.01 – 193.97 |
| high | 500 | 191 | 20 | 25 | 198.9840 / 198.9840 | 4.825 | 193.67 – 214.83 |

Largest deviation across all 21 compared quantities: **0.0005**, which is R's own
printing precision.

The test also asserts properties that must hold for any input: the rarefaction
curve passes through `S_obs` at `t = T`, is monotonic, never exceeds the Chao2
asymptote under extrapolation, the confidence interval contains the estimate and
never falls below `S_obs`, and the interval is asymmetric — it is the
log-transformed form, not `S_chao2 ± 1.96 · SE`.

## 3. What the tests exercise

Being explicit about coverage, since "the pipeline runs" and "the pipeline is
correct" are different claims.

| Process | How it was checked |
|---|---|
| `DBCAN_ANNOTATE` | run on two genomes; output compared gene by gene against mpCGCdb (section 1) |
| `DBCAN_CGC` | same; cluster boundaries compared, and the parameter choices were derived from that comparison |
| `CGC_CATALOG` | same; the compared catalog is its output |
| `CGC_PROTEINS` | run on two genomes — which is what exposed a filename collision when collecting more than one |
| `COLOC_NETWORK` | run on the two-genome catalog |
| `DIAMOND_ALLVSALL`, `LEIDEN_COMMUNITIES` | run on five families |
| `MUSCLE_ALIGN`, `FASTTREE` | run on five families |
| `CGC_FINGERPRINT` | run on a real catalog; subset collapse checked by hand (73 distinct fingerprints reduced to 52 maximal, with `GH177` folding into `GH177,S1_14`) |
| `CGC_CLUSTER`, `CGC_COOCCURRENCE`, `SUMMARY_METRICS` | run on a real catalog |
| `CGC_RICHNESS` | estimators checked numerically against R (section 2); the process itself run on a catalog with too few genomes, confirming it skips rather than reporting a meaningless asymptote |
| `DBCAN_DATABASE` | the underlying `run_dbcan database` call was run to completion to build the databases these tests used; the process wrapper itself is checked only by `-preview` |
| `GGTREE_FIGURE` | run as a Nextflow process against an R 4.5.3 / ggtree 4.0.4 environment, rendering all three mined families (GH10, GT2, S1_7) coloured by taxonomy, PNG and SVG each. The script was separately exercised in all four colouring modes: `--color-by taxonomy` at `--rank phylum` and at `--rank genus`, `--color-by community`, and with no taxonomy supplied, which falls back to community colouring with a message rather than failing |

Diversity estimates were verified against R but have not been checked end to end
against the published mpCGCdb figures, which would need the full 22,607-genome
catalog rather than a single genome.

Running the figure step as a process, rather than only running its R script,
found three defects that the script test could not:

* the process pointed at a `bioconductor-ggtreeextra` container for a script that
  never loads ggtreeExtra, and no biocontainer carries `r-optparse`, so the docker
  profile would have failed. Removing the `optparse` dependency in favour of base-R
  argument parsing let the step use the stock `bioconductor-ggtree` image
* ggplot2 4 requires the `svglite` package to write SVG and no longer falls back to
  the cairo device, so the task wrote its PNG and then aborted. SVG now tries
  svglite, then cairo, and keeps the PNG if neither is available
* **`GGTREE_FIGURE` rendered one figure and reported success.** Its `protein_map`
  and metadata inputs were single-item *queue* channels; a process mixing an
  N-item channel with a one-item queue channel runs once, and the remaining
  families are dropped with no error. Both are now value channels. This could not
  appear with a single family, and produced no failure signal — the run exited 0
  with two thirds of its figures missing

The last of these is the reason the mining dataflow is now tested with three
families rather than one.

## Reproducing these checks

```bash
# databases (~8 GB, once)
nextflow run . --step download_db --outdir refs -profile conda

# identification on the validation genome
nextflow run . -profile test,conda --db refs/db/dbcan_db --outdir results

# comparison, given the mpCGCdb rows for that genome
python test/validate_against_mpcgcdb.py \
    --run results/catalog/all_cgc_catalog.tsv \
    --reference test/reference/mpcgcdb_GCA_000325705.1_ASM32570v1_genomic.tsv.gz \
    --genome GCA_000325705.1_ASM32570v1_genomic \
    --report validation_report.txt

# estimator check (no R required)
python test/test_inext.py
```

The reference rows for both genomes ship with the repository under
`test/reference/` (15 KB and 3 KB), so the comparisons run without downloading
anything; the full catalog is on the downloads page of
[mpcgcdb.com](https://mpcgcdb.com). `test/samplesheet_test.csv`
points at the assembly, which is
[GCA_000325705.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000325705.1/)
at NCBI.
