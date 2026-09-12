/*
 * Row-conditioned asymmetric colocalisation matrix:
 *   M[anchor, colocalized] = P(colocalized family present | anchor present)
 * Families below --cooccurrence_min_cgc distinct CGCs are excluded. Rows and
 * columns are ordered by average-linkage clustering of Jaccard distance.
 */
process CGC_COOCCURRENCE {
    label 'medium'
    label 'pydata'
    tag   'colocalisation matrix'

    conda     "${moduleDir}/../../env/pydata.yml"
    container 'quay.io/biocontainers/scikit-learn:1.4.2'

    input:
    path fingerprints

    output:
    path 'cooccurrence_matrix.tsv'      , emit: matrix
    path 'cooccurrence_order.txt'       , emit: order
    path 'cooccurrence_heatmap.{png,svg}', emit: figure, optional: true
    path 'versions.yml'                 , emit: versions

    script:
    """
    mpcgc_cooccurrence.py \
        --fingerprints ${fingerprints} \
        --min-cgc ${params.cooccurrence_min_cgc} \
        --out-matrix cooccurrence_matrix.tsv \
        --out-order cooccurrence_order.txt \
        --figure cooccurrence_heatmap

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
    END_VERSIONS
    """
}
