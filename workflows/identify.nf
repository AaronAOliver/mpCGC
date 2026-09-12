/*
 * IDENTIFY - orange dataflow of Figure 1.
 *
 *   genome fasta
 *      -> Pyrodigal gene calls                        (DBCAN_ANNOTATE)
 *      -> protein fasta
 *      -> DIAMOND vs CAZy / SulfAtlas / TCDB / PRODORIC / MEROPS
 *         + PyHMMER vs dbCAN, dbCAN-sub, STP          (DBCAN_ANNOTATE)
 *      -> annotated gff3 + seed-and-extend CGC calls  (DBCAN_CGC)
 *      -> CGC catalog (csv)                           (CGC_CATALOG)
 *      -> per-family protein fasta                    (CGC_PROTEINS)
 */

include { DBCAN_ANNOTATE } from '../modules/local/dbcan_annotate'
include { DBCAN_CGC      } from '../modules/local/dbcan_cgc'
include { CGC_CATALOG    } from '../modules/local/cgc_catalog'
include { CGC_PROTEINS   } from '../modules/local/cgc_proteins'

workflow IDENTIFY {

    take:
    genomes     // channel: [ meta, fasta, gff ]  (gff may be [])
    db          // value:   dbCAN database directory

    main:
    ch_versions = channel.empty()

    DBCAN_ANNOTATE(genomes, db)
    ch_versions = ch_versions.mix(DBCAN_ANNOTATE.out.versions.first())

    DBCAN_CGC(DBCAN_ANNOTATE.out.annotated, db)
    ch_versions = ch_versions.mix(DBCAN_CGC.out.versions.first())

    // Each per-genome TSV is named <sample>.cgc_standard_out.tsv, so the
    // catalog builder recovers the genome name from the filename.
    CGC_CATALOG(
        DBCAN_CGC.out.cgc_standard.map { _meta, tsv -> tsv }.collect()
    )

    CGC_PROTEINS(
        CGC_CATALOG.out.catalog,
        DBCAN_ANNOTATE.out.faa.map { _meta, faa -> faa }.collect()
    )

    emit:
    catalog     = CGC_CATALOG.out.catalog
    summary     = CGC_CATALOG.out.per_genome
    proteins    = CGC_PROTEINS.out.families
    protein_map = CGC_PROTEINS.out.protein_map
    faa         = DBCAN_ANNOTATE.out.faa
    versions    = ch_versions
}
