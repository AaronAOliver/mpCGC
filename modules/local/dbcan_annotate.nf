/*
 * Gene calling + every homology search of the identification dataflow.
 *
 * mpcgc_rename_ids.py          : private copy of the input for dbCAN, with
 *                                short placeholder IDs if the originals are
 *                                too long for dbCAN (see that script)
 * run_dbcan CAZyme_annotation  : Pyrodigal gene calls, DIAMOND vs CAZy,
 *                                PyHMMER vs dbCAN and dbCAN-sub
 * run_dbcan gff_process        : DIAMOND vs TCDB (transporters, SusC/SusD),
 *                                SulfAtlas (sulfatases), MEROPS (peptidases),
 *                                PRODORIC (transcription factors) and PyHMMER
 *                                vs the signal-transduction HMMs; writes the
 *                                annotated cgc.gff used for CGC calling.
 *
 * Outputs still carry the placeholder IDs; DBCAN_CGC restores the originals
 * after cluster calling and is the step that publishes the results.
 */
process DBCAN_ANNOTATE {
    label 'bigmem'
    label 'dbcan'
    tag   "${meta.id}"

    conda     "${moduleDir}/../../env/dbcan.yml"

    input:
    // staged in input/ so no output written here can share a name with (and
    // write through the link to) the user's original files
    tuple val(meta), path(fasta, stageAs: 'input/*'), path(gff, stageAs: 'input/*')
    path  db

    output:
    tuple val(meta), path("${meta.id}.dbcan/"), path("${meta.id}.id_map.tsv"), emit: annotated
    path  'versions.yml'                                                     , emit: versions

    script:
    // true/false given on the command line arrive as text ('false' would count
    // as true), so every boolean parameter is read with toString().toBoolean()
    def methods = params.dbcansub.toString().toBoolean() ? params.methods : params.methods.replaceAll(/,?dbCANsub/, '')
    // in protein mode dbCAN cannot call genes, so an annotation GFF is required
    def protein = params.mode == 'protein'
    def gff_in  = protein ? "--gff ${gff} --out-gff mpcgc_input.gff" : ''
    def gff_arg = protein ? '--input_gff mpcgc_input.gff' : ''
    """
    # dbCAN works on this copy (renamed only if IDs are too long for it), never
    # on the user's files, which Nextflow stages as links
    mpcgc_rename_ids.py \
        --mode ${params.mode} \
        --sample ${meta.id} \
        --fasta ${fasta} ${gff_in} \
        --out-fasta mpcgc_input.fasta \
        --map ${meta.id}.id_map.tsv

    run_dbcan CAZyme_annotation \
        --mode ${params.mode} \
        --input_raw_data mpcgc_input.fasta \
        --output_dir ${meta.id}.dbcan \
        --db_dir ${db} \
        --methods ${methods} \
        --e_value_threshold ${params.e_value_cazyme} \
        --threads ${task.cpus}

    run_dbcan gff_process \
        --output_dir ${meta.id}.dbcan \
        --db_dir ${db} \
        --gff_type ${params.gff_type} \
        --threads ${task.cpus} \
        --e_value_threshold_tc ${params.e_value_tc} \
        --coverage_threshold_tc ${params.coverage_tc} \
        --e_value_threshold_tf_diamond ${params.e_value_tf} \
        --coverage_threshold_tf_diamond ${params.coverage_tf} \
        --e_value_threshold_stp ${params.e_value_stp} \
        --coverage_threshold_stp ${params.coverage_stp} \
        --e_value_threshold_sulfatase ${params.e_value_sulfatase} \
        --coverage_threshold_sulfatase ${params.coverage_sulfatase} \
        --e_value_threshold_peptidase ${params.e_value_peptidase} \
        --coverage_threshold_peptidase ${params.coverage_peptidase} \
        ${gff_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(python -c "import importlib.metadata as m; print(m.version('dbcan'))")
        pyrodigal: \$(python -c "import pyrodigal; print(pyrodigal.__version__)")
        diamond: \$(diamond --version 2>&1 | sed 's/^diamond version //')
        pyhmmer: \$(python -c "import pyhmmer; print(pyhmmer.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}.dbcan
    touch ${meta.id}.dbcan/uniInput.faa ${meta.id}.dbcan/overview.tsv
    printf 'placeholder\\toriginal\\n' > ${meta.id}.id_map.tsv
    echo '"${task.process}": {}' > versions.yml
    """
}
