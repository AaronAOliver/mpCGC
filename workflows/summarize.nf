/*
 * SUMMARIZE - blue dataflow of Figure 1.
 *
 *   CGC catalog
 *      -> cluster encoding (binary S1/GH/PL/CE fingerprints)  (CGC_FINGERPRINT)
 *      -> hierarchical clustering of pairwise Hamming distance (CGC_CLUSTER)
 *      -> asymmetric colocalisation matrix                     (CGC_COOCCURRENCE)
 *      -> CGC diversity: iNEXT rarefaction + Chao2 asymptotes  (CGC_RICHNESS)
 *      -> summary metrics (csv)                                (SUMMARY_METRICS)
 */

include { CGC_FINGERPRINT  } from '../modules/local/cgc_fingerprint'
include { CGC_CLUSTER      } from '../modules/local/cgc_cluster'
include { CGC_COOCCURRENCE } from '../modules/local/cgc_cooccurrence'
include { CGC_RICHNESS     } from '../modules/local/cgc_richness'
include { SUMMARY_METRICS  } from '../modules/local/summary_metrics'

workflow SUMMARIZE {

    take:
    catalog     // value: all_cgc_catalog.tsv
    metadata    // value: per-genome metadata TSV

    main:
    ch_versions = channel.empty()

    CGC_FINGERPRINT(catalog, metadata)
    ch_versions = ch_versions.mix(CGC_FINGERPRINT.out.versions)

    CGC_CLUSTER(CGC_FINGERPRINT.out.fingerprints)
    ch_versions = ch_versions.mix(CGC_CLUSTER.out.versions)

    CGC_COOCCURRENCE(CGC_FINGERPRINT.out.fingerprints)
    ch_versions = ch_versions.mix(CGC_COOCCURRENCE.out.versions)

    CGC_RICHNESS(CGC_FINGERPRINT.out.fingerprints)
    ch_versions = ch_versions.mix(CGC_RICHNESS.out.versions)

    SUMMARY_METRICS(
        catalog,
        CGC_FINGERPRINT.out.fingerprints,
        CGC_RICHNESS.out.asymptotes
    )

    emit:
    fingerprints = CGC_FINGERPRINT.out.fingerprints
    richness     = CGC_RICHNESS.out.curves
    asymptotes   = CGC_RICHNESS.out.asymptotes
    summary      = SUMMARY_METRICS.out.summary
    versions     = ch_versions
}
