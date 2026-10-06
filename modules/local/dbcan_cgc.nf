/*
 * CGC determination - the seed-and-extend step.
 *
 * Clusters are seeded on co-localised signature genes (CAZymes plus, optionally,
 * transporters / transcription factors / signal-transduction proteins) and
 * extended across at most --num_null_gene unannotated genes, then optionally
 * widened by a user-provided gene-boundary extension limit (--extend_mode).
 *
 * Works on its own copy of DBCAN_ANNOTATE's output, then restores the original
 * sequence IDs in every file (mpcgc_restore_ids.py) if DBCAN_ANNOTATE had to
 * give dbCAN short placeholders. This is the step that publishes the
 * per-genome dbCAN results, so no placeholder reaches the results folder.
 */
process DBCAN_CGC {
    label 'low'
    label 'dbcan'
    tag   "${meta.id}"

    conda     "${moduleDir}/../../env/dbcan.yml"

    input:
    tuple val(meta), path(annot_dir), path(id_map)
    path  db

    output:
    tuple val(meta), path("${meta.id}/")                    , emit: annotated
    // named per sample so several genomes can be collected without colliding
    tuple val(meta), path("${meta.id}.faa")                 , emit: faa
    tuple val(meta), path("${meta.id}.overview.tsv")        , emit: overview
    tuple val(meta), path("${meta.id}.cgc_standard_out.tsv"), emit: cgc_standard
    tuple val(meta), path("${meta.id}.cgc.gff")             , emit: cgc_gff
    tuple val(meta), path("${meta.id}.total_cgc_info.tsv")  , optional: true, emit: cgc_info
    // placeholder -> original ID table, only when IDs had to be shortened
    tuple val(meta), path("${meta.id}.renamed_ids.tsv")     , optional: true, emit: id_map
    path  'versions.yml'                                    , emit: versions

    script:
    // Accessory classes a cluster must carry besides its min_core_cazyme
    // CAZymes. The default, TC with 'all' logic, requires a transporter, which
    // reproduces the mpCGCdb clusters; see --additional_genes in the README.
    def classes = params.additional_genes.tokenize(',').collect { c -> c.trim() }.findAll { c -> c }
    def add = classes
        ? classes.collect { c -> "--additional_genes ${c}" }.join(' ') +
          " --additional_logic ${params.additional_logic}" +
          " --additional_min_categories ${params.additional_min_categories}"
        : ''
    def dist = params.use_distance.toString().toBoolean() ? "--use_distance --base_pair_distance ${params.base_pair_distance}" : ''
    def ext  = params.extend_mode == 'gene' ? "--extend_mode gene --extend_gene_count ${params.extend_gene_count}"
             : params.extend_mode == 'bp'   ? "--extend_mode bp --extend_bp ${params.extend_bp}"
             : '--extend_mode none'
    def nulls = params.num_null_gene > 0 ? "--use_null_genes --num_null_gene ${params.num_null_gene}"
                                         : '--no-use_null_genes'
    """
    # own copy: the input folder belongs to DBCAN_ANNOTATE's task
    cp -rL ${annot_dir} ${meta.id}

    run_dbcan cgc_finder \
        --output_dir ${meta.id} \
        --min_cluster_genes ${params.min_cluster_genes} \
        --min_core_cazyme ${params.min_core_cazyme} \
        ${nulls} \
        ${add} \
        ${ext} \
        ${dist}

    # put the original protein IDs / contig names back in every file
    mpcgc_restore_ids.py --map ${id_map} ${meta.id}

    cp ${meta.id}/uniInput.faa         ${meta.id}.faa
    cp ${meta.id}/overview.tsv         ${meta.id}.overview.tsv
    cp ${meta.id}/cgc_standard_out.tsv ${meta.id}.cgc_standard_out.tsv
    cp ${meta.id}/cgc.gff              ${meta.id}.cgc.gff
    if [ -s ${meta.id}/total_cgc_info.tsv ]; then
        cp ${meta.id}/total_cgc_info.tsv ${meta.id}.total_cgc_info.tsv
    fi
    # publish the ID table only if renaming happened (header line only otherwise)
    if [ \$(wc -l < ${id_map}) -gt 1 ]; then cp ${id_map} ${meta.id}.renamed_ids.tsv; fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(python -c "import importlib.metadata as m; print(m.version('dbcan'))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}
    touch ${meta.id}.faa ${meta.id}.overview.tsv
    printf 'CGC#\tGene Type\tContig ID\tProtein ID\tGene Start\tGene Stop\tGene Strand\tGene Annotation\n' > ${meta.id}.cgc_standard_out.tsv
    touch ${meta.id}.cgc.gff
    echo '"${task.process}": {}' > versions.yml
    """
}
