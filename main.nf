#!/usr/bin/env nextflow
/*
 * mpCGC - marine polysaccharide CAZyme gene clusterer
 *
 * Three dataflows (Figure 1 of the manuscript):
 *   identify   fasta -> gene calls -> CAZyme/accessory annotation -> CGC catalog
 *   mine       CGC catalog -> colocalisation networks, SSN communities, phylogenies
 *   summarize  CGC catalog -> cluster encoding, hierarchical clustering, CGC diversity
 */

nextflow.enable.dsl = 2

include { DOWNLOAD_DB } from './workflows/download_db'
include { IDENTIFY    } from './workflows/identify'
include { MINE        } from './workflows/mine'
include { SUMMARIZE   } from './workflows/summarize'

def helpMessage() {
    log.info """
    ================================================================
     mpCGC  v${workflow.manifest.version}
     marine polysaccharide CAZyme gene clusterer
    ================================================================

    Usage:
      nextflow run AaronAOliver/mpCGC --input samplesheet.csv --db <dbcan_db> -profile conda

    Required:
      --input            Samplesheet CSV: sample,fasta[,taxonomy][,gff]
                         gff is required only when --mode protein
      --metadata         Optional TSV of per-genome taxonomy, overriding the
                         samplesheet column (genome id + a *taxonomy* column)
      --db               dbCAN database directory (build it with --step download_db)
                         --db_from_s3 pulls a pinned release instead of db_current

    Main options:
      --step             identify | mine | summarize | all | download_db  [${params.step}]
      --outdir           Output directory                                  [${params.outdir}]
      --mode             prok | meta | protein                             [${params.mode}]

    CGC calling:
      --min_core_cazyme      Minimum CAZymes per cluster        [${params.min_core_cazyme}]
      --num_null_gene        Max intervening unannotated genes  [${params.num_null_gene}]
      --additional_genes     Signature classes that can satisfy a cluster
                             [${params.additional_genes}]
      --additional_logic     all | any                          [${params.additional_logic}]
      --extend_mode          Gene-boundary extension limit: none | bp | gene [${params.extend_mode}]
      --extend_gene_count    Flanking genes when --extend_mode gene  [${params.extend_gene_count}]
      --extend_bp            Flanking bp when --extend_mode bp       [${params.extend_bp}]
      --keep_ecami           Keep dbCAN-sub eCAMI indices (GH140_e33) [${params.keep_ecami}]

    Mining:
      --ssn_evalue           DIAMOND all-vs-all E-value     [${params.ssn_evalue}]
      --leiden_resolution    Leiden resolution              [${params.leiden_resolution}]
      --make_trees           MUSCLE + FastTree per family   [${params.make_trees}]
      --make_tree_figures    Render trees with ggtree        [${params.make_tree_figures}]
      --tree_color_by        taxonomy | community            [${params.tree_color_by}]
      --tree_color_rank      phylum | class | order | family | genus [${params.tree_color_rank}]

    Summary:
      --collapse_subsets     Count maximal fingerprints only [${params.collapse_subsets}]
      --group_by             Metadata column for lineages    [${params.group_by}]
      --min_genomes          Lineages below this skip rarefaction [${params.min_genomes}]

    Profiles: conda, mamba, docker, singularity, slurm, local_envs, test, test_full

    Full documentation: docs/usage.md
    """.stripIndent()
}

workflow {

    if (params.help) { helpMessage(); return }

    // ---- database-only entry point -------------------------------------
    if (params.step == 'download_db') {
        DOWNLOAD_DB()
        return
    }

    if (!params.db) {
        error "No --db given. Build the databases first:\n" +
              "  nextflow run . --step download_db --outdir <dir>"
    }
    def db_ch = channel.fromPath(params.db, type: 'dir', checkIfExists: true).first()

    // ---- samplesheet ---------------------------------------------------
    if (!params.input) { error "No --input samplesheet given. See docs/usage.md" }

    def genomes = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
            if (!row.fasta)  error "Samplesheet row missing 'fasta': ${row}"
            if (params.mode == 'protein' && !row.gff) {
                error "--mode protein needs a 'gff' column for sample ${row.sample}"
            }
            def meta = [ id: row.sample, taxonomy: row.taxonomy ?: '' ]
            def fa   = file(row.fasta, checkIfExists: true)
            def gff  = row.gff ? file(row.gff, checkIfExists: true) : []
            tuple(meta, fa, gff)
        }

    // ---- metadata table (optional; used for lineage grouping) ----------
    def metadata = params.metadata
        ? channel.fromPath(params.metadata, checkIfExists: true).first()
        : genomes.map { meta, _fa, _gff -> "${meta.id}\t${meta.taxonomy}" }
                 .collectFile(name: 'sample_taxonomy.tsv', newLine: true,
                              seed: "sample\ttaxonomy")

    // ---- dataflows -----------------------------------------------------
    def run_identify  = params.step in ['identify', 'all']
    def run_mine      = params.step in ['mine', 'all']
    def run_summarize = params.step in ['summarize', 'all']

    def catalog
    def proteins
    def protein_map

    if (run_identify) {
        IDENTIFY(genomes, db_ch)
        catalog     = IDENTIFY.out.catalog
        proteins    = IDENTIFY.out.proteins
        protein_map = IDENTIFY.out.protein_map
    } else {
        // resume from a previous run's published catalog
        def prev = file("${params.outdir}/catalog/all_cgc_catalog.tsv")
        if (!prev.exists()) {
            error "--step ${params.step} needs an existing catalog at ${prev}.\n" +
                  "Run --step identify first, or point --outdir at a completed run."
        }
        catalog     = channel.value(prev)
        proteins    = channel.fromPath("${params.outdir}/catalog/proteins/*.faa").collect()
        protein_map = channel.fromPath("${params.outdir}/catalog/proteins/protein_map.tsv",
                                       checkIfExists: true).first()
    }

    if (run_mine)      MINE(catalog, proteins, protein_map, metadata)
    if (run_summarize) SUMMARIZE(catalog, metadata)
}
