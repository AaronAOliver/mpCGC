/*
 * Cluster encoding: each CGC becomes a binary presence/absence fingerprint over
 * S1 / GH / PL / CE families (one token per gene, priority GH > PL > CE > S1;
 * AA, GT and CBM modules are not encoded). With --collapse_subsets, a
 * fingerprint that is a perfect subset of another is folded into its maximal
 * superset, so nested clusters are not counted as distinct functions.
 */
process CGC_FINGERPRINT {
    label 'medium'
    label 'pydata'
    tag   'cluster encoding'

    conda     "${moduleDir}/../../env/pydata.yml"
    container 'quay.io/biocontainers/pandas:2.2.1'

    input:
    path catalog
    path metadata

    output:
    path 'cgc_fingerprints.tsv'      , emit: fingerprints
    path 'fingerprint_collapse.tsv'  , emit: collapse_map
    path 'versions.yml'              , emit: versions

    script:
    def collapse = params.collapse_subsets ? "--collapse --ambiguous ${params.ambiguous_subset}" : ''
    """
    mpcgc_fingerprint.py \
        --catalog ${catalog} \
        --metadata ${metadata} \
        --group-by ${params.group_by} \
        --out cgc_fingerprints.tsv \
        --collapse-map fingerprint_collapse.tsv \
        ${collapse}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
