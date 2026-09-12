/*
 * MUSCLE v5 alignment of one enzyme family.
 *
 * MUSCLE's default algorithm is accurate but scales poorly, so families at or
 * above --muscle_super5_min sequences are aligned with --super5 instead. The
 * count is taken in the shell rather than in Groovy: countFasta() operates on
 * channel values, not on a path already staged into a task.
 */
process MUSCLE_ALIGN {
    label 'medium'
    label 'phylo'
    tag   "${meta.id}"

    conda     "${moduleDir}/../../env/phylo.yml"
    container 'quay.io/biocontainers/muscle:5.1--h9f5acd7_1'

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path("${meta.id}.aln.faa"), emit: alignment
    path 'versions.yml'                        , emit: versions

    script:
    """
    n_seqs=\$(grep -c '^>' ${faa} || true)
    if [ "\${n_seqs}" -ge ${params.muscle_super5_min} ]; then
        algo=-super5
    else
        algo=-align
    fi
    echo "aligning ${meta.id}: \${n_seqs} sequences with muscle \${algo}" >&2

    muscle \${algo} ${faa} -output ${meta.id}.aln.faa -threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        muscle: \$(muscle -version 2>&1 | sed 's/^muscle //' | cut -d' ' -f1)
        n_sequences: \${n_seqs}
    END_VERSIONS
    """
}
