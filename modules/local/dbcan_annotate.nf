/*
 * Gene calling + every homology search of the identification dataflow.
 *
 * run_dbcan CAZyme_annotation  : Pyrodigal gene calls, DIAMOND vs CAZy,
 *                                PyHMMER vs dbCAN and dbCAN-sub
 * run_dbcan gff_process        : DIAMOND vs TCDB (transporters, SusC/SusD),
 *                                SulfAtlas (sulfatases), MEROPS (peptidases),
 *                                PRODORIC (transcription factors) and PyHMMER
 *                                vs the signal-transduction HMMs; writes the
 *                                annotated cgc.gff used for CGC calling.
 */
process DBCAN_ANNOTATE {
    label 'bigmem'
    label 'dbcan'
    tag   "${meta.id}"

    conda     "${moduleDir}/../../env/dbcan.yml"
    container 'quay.io/biocontainers/dbcan:5.2.9--pyhdfd78af_0'

    input:
    tuple val(meta), path(fasta), path(gff)
    path  db

    output:
    tuple val(meta), path("${meta.id}/")               , emit: annotated
    // named per sample so several genomes can be collected without colliding
    tuple val(meta), path("${meta.id}.faa")            , emit: faa
    tuple val(meta), path("${meta.id}.overview.tsv")   , emit: overview
    path  'versions.yml'                               , emit: versions

    script:
    def methods = params.dbcansub ? params.methods : params.methods.replaceAll(/,?dbCANsub/, '')
    // in protein mode dbCAN cannot call genes, so an annotation GFF is required
    def gff_arg = params.mode == 'protein' ? "--input_gff ${gff}" : ''
    """
    run_dbcan CAZyme_annotation \
        --mode ${params.mode} \
        --input_raw_data ${fasta} \
        --output_dir ${meta.id} \
        --db_dir ${db} \
        --methods ${methods} \
        --e_value_threshold ${params.e_value_cazyme} \
        --threads ${task.cpus}

    run_dbcan gff_process \
        --output_dir ${meta.id} \
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

    cp ${meta.id}/uniInput.faa  ${meta.id}.faa
    cp ${meta.id}/overview.tsv  ${meta.id}.overview.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(python -c "import importlib.metadata as m; print(m.version('dbcan'))")
        pyrodigal: \$(python -c "import pyrodigal; print(pyrodigal.__version__)")
        diamond: \$(diamond --version 2>&1 | sed 's/^diamond version //')
        pyhmmer: \$(python -c "import pyhmmer; print(pyhmmer.__version__)")
    END_VERSIONS
    """
}
