/* Merge per-genome cgc_standard_out.tsv into one catalog keyed by genome. */
process CGC_CATALOG {
    label 'low'
    label 'pydata'
    tag   'CGC catalog'

    conda     "${moduleDir}/../../env/pydata.yml"
    container 'quay.io/biocontainers/pandas:2.2.1'

    input:
    path cgc_tsvs, stageAs: 'cgc/*'

    output:
    path 'all_cgc_catalog.tsv'  , emit: catalog
    path 'cgc_per_genome.tsv'   , emit: per_genome
    path 'versions.yml'         , emit: versions

    script:
    def ecami = params.keep_ecami ? '--keep-ecami' : ''
    """
    mpcgc_catalog.py \
        --inputs cgc/*.cgc_standard_out.tsv \
        --catalog all_cgc_catalog.tsv \
        --per-genome cgc_per_genome.tsv \
        ${ecami}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
