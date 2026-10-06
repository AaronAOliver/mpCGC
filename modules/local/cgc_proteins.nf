/* Split CGC member proteins into one fasta per enzyme family. */
process CGC_PROTEINS {
    label 'medium'
    label 'pydata'
    tag   'per-family protein fasta'

    conda     "${moduleDir}/../../env/pydata.yml"

    input:
    path catalog
    path faas, stageAs: 'faa/*'

    output:
    // optional: a run with no CGCs has no families, which is a result, not an error
    path 'families/*.faa', emit: families, optional: true
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

    stub:
    """
    mkdir -p families
    for i in 1 2 3 4; do printf ">p\$i\nMKV\n" >> families/GH16.faa; done
    touch family_counts.tsv protein_map.tsv
    echo '"${task.process}": {}' > versions.yml
    """
}
