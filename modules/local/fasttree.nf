/* FastTree maximum-likelihood phylogeny from the family alignment. */
process FASTTREE {
    label 'medium'
    label 'phylo'
    tag   "${meta.id}"

    conda     "${moduleDir}/../../env/phylo.yml"

    input:
    tuple val(meta), path(aln)

    output:
    tuple val(meta), path("${meta.id}.nwk"), emit: tree
    path "${meta.id}.fasttree.log"         , emit: log
    path 'versions.yml'                    , emit: versions

    script:
    def ft = params.fasttree_threads_bin ?: 'FastTreeMP'
    """
    export OMP_NUM_THREADS=${task.cpus}
    ${ft} -nosupport ${params.fasttree_model} ${aln} \
        > ${meta.id}.nwk 2> ${meta.id}.fasttree.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasttree: \$(${ft} -help 2>&1 | head -1 | sed 's/^FastTree //' | cut -d' ' -f1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.nwk ${meta.id}.fasttree.log
    echo '"${task.process}": {}' > versions.yml
    """
}
