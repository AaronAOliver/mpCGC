/* Final per-genome and per-lineage summary table (blue dataflow output). */
process SUMMARY_METRICS {
    label 'low'
    label 'pydata'
    tag   'summary metrics'

    conda     "${moduleDir}/../../env/pydata.yml"
    container 'quay.io/biocontainers/pandas:2.2.1'

    input:
    path catalog
    path fingerprints
    path asymptotes

    output:
    path 'mpcgc_summary.csv'        , emit: summary
    path 'mpcgc_summary_lineage.csv', emit: lineage
    path 'versions.yml'             , emit: versions

    script:
    """
    mpcgc_summary.py \
        --catalog ${catalog} \
        --fingerprints ${fingerprints} \
        --asymptotes ${asymptotes} \
        --out mpcgc_summary.csv \
        --out-lineage mpcgc_summary_lineage.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
