/*
 * Hierarchical clustering of CGCs by pairwise Hamming distance between their
 * binary fingerprints. Two clusters are treated as overlapping in function when
 * one contains a superset of all annotations of the other.
 */
process CGC_CLUSTER {
    label 'high'
    label 'pydata'
    tag   'hierarchical clustering'

    conda     "${moduleDir}/../../env/pydata.yml"
    container 'quay.io/biocontainers/scikit-learn:1.4.2'

    input:
    path fingerprints

    output:
    path 'cgc_linkage.npy'        , emit: linkage
    path 'cgc_clusters.tsv'       , emit: clusters
    path 'cgc_dendrogram.{png,svg}', emit: figure, optional: true
    path 'versions.yml'           , emit: versions

    script:
    """
    mpcgc_cluster.py \
        --fingerprints ${fingerprints} \
        --out-linkage cgc_linkage.npy \
        --out-clusters cgc_clusters.tsv \
        --figure cgc_dendrogram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
    END_VERSIONS
    """
}
