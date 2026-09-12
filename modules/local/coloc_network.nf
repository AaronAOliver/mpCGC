/* Enzyme colocalisation network: nodes are families, edges are CGC co-membership. */
process COLOC_NETWORK {
    label 'medium'
    label 'pydata'
    tag   'colocalisation network'

    conda     "${moduleDir}/../../env/pydata.yml"
    container 'quay.io/biocontainers/pandas:2.2.1'

    input:
    path catalog

    output:
    path 'colocalization_network.json', emit: network
    path 'colocalization_edges.tsv'   , emit: edges
    path 'versions.yml'               , emit: versions

    script:
    """
    mpcgc_coloc_network.py \
        --catalog ${catalog} \
        --edges colocalization_edges.tsv \
        --json colocalization_network.json \
        --min-support ${params.coloc_min_support}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
