/* Leiden community detection on the sequence similarity network (weight = bitscore). */
process LEIDEN_COMMUNITIES {
    label 'medium'
    label 'network'
    tag   "${meta.id}"

    conda     "${moduleDir}/../../env/network.yml"
    container 'quay.io/biocontainers/python-igraph:0.11.5--py311h0e8b6ba_0'

    input:
    tuple val(meta), path(ssn)

    output:
    tuple val(meta), path("${meta.id}.communities.tsv"), emit: communities
    tuple val(meta), path("${meta.id}.network.json")   , emit: network
    path 'versions.yml'                                , emit: versions

    script:
    """
    mpcgc_leiden.py \
        --ssn ${ssn} \
        --family ${meta.id} \
        --evalue ${params.ssn_evalue} \
        --resolution ${params.leiden_resolution} \
        --out-communities ${meta.id}.communities.tsv \
        --out-json ${meta.id}.network.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python-igraph: \$(python3 -c "import igraph; print(igraph.__version__)")
    END_VERSIONS
    """
}
