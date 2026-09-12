/* Split CGC member proteins into one fasta per enzyme family. */
process CGC_PROTEINS {
    label 'medium'
    label 'pydata'
    tag   'per-family protein fasta'

    conda     "${moduleDir}/../../env/pydata.yml"
    container 'quay.io/biocontainers/pandas:2.2.1'

    input:
    path catalog
    path faas, stageAs: 'faa/*'

    output:
    path 'families/*.faa', emit: families
    path 'family_counts.tsv', emit: counts
    path 'protein_map.tsv'  , emit: protein_map
    path 'versions.yml'   , emit: versions

    script:
    """
    mkdir -p families
    mpcgc_extract_proteins.py \
        --catalog ${catalog} \
        --faa faa/* \
        --outdir families \
        --counts family_counts.tsv \
        --protein-map protein_map.tsv \
        --min-seqs 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
