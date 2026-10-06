/* All-vs-all DIAMOND within one enzyme family -> sequence similarity network. */
process DIAMOND_ALLVSALL {
    label 'medium'
    label 'dbcan'
    tag   "${meta.id}"

    conda     "${moduleDir}/../../env/dbcan.yml"

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path("${meta.id}.ssn.tsv.gz"), emit: hits
    path 'versions.yml'                           , emit: versions

    script:
    """
    diamond makedb --in ${faa} --db ${meta.id}.dmnd --threads ${task.cpus} --quiet

    diamond blastp \
        --query ${faa} \
        --db ${meta.id}.dmnd \
        --out ${meta.id}.ssn.tsv \
        --outfmt 6 \
        --max-target-seqs ${params.ssn_max_targets} \
        --evalue ${params.ssn_evalue} \
        --threads ${task.cpus} \
        --quiet
    gzip -n ${meta.id}.ssn.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | sed 's/^diamond version //')
    END_VERSIONS
    """

    stub:
    """
    echo -n | gzip > ${meta.id}.ssn.tsv.gz
    echo '"${task.process}": {}' > versions.yml
    """
}
