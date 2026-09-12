/*
 * CGC determination - the seed-and-extend step.
 *
 * Clusters are seeded on co-localised signature genes (CAZymes plus, optionally,
 * transporters / transcription factors / signal-transduction proteins) and
 * extended across at most --num_null_gene unannotated genes, then optionally
 * widened by a user-provided gene-boundary extension limit (--extend_mode).
 */
process DBCAN_CGC {
    label 'low'
    label 'dbcan'
    tag   "${meta.id}"

    conda     "${moduleDir}/../../env/dbcan.yml"
    container 'quay.io/biocontainers/dbcan:5.2.9--pyhdfd78af_0'

    input:
    tuple val(meta), path(annot_dir)
    path  db

    output:
    tuple val(meta), path("${meta.id}.cgc_standard_out.tsv"), emit: cgc_standard
    tuple val(meta), path("${meta.id}.cgc.gff")             , emit: cgc_gff
    tuple val(meta), path("${meta.id}.total_cgc_info.tsv")  , optional: true, emit: cgc_info
    path  'versions.yml'                                    , emit: versions

    script:
    // Signature classes that can satisfy a cluster. mpCGCdb lists CAZyme
    // alongside the accessory classes under 'any' logic, so a cluster needs
    // min_core_cazyme CAZymes plus at least one signature class - which keeps
    // CAZyme-only clusters without requiring a transporter.
    def classes = params.additional_genes.tokenize(',').collect { c -> c.trim() }.findAll { c -> c }
    def add = classes
        ? classes.collect { c -> "--additional_genes ${c}" }.join(' ') +
          " --additional_logic ${params.additional_logic}" +
          " --additional_min_categories ${params.additional_min_categories}"
        : ''
    def dist = params.use_distance ? "--use_distance --base_pair_distance ${params.base_pair_distance}" : ''
    def ext  = params.extend_mode == 'gene' ? "--extend_mode gene --extend_gene_count ${params.extend_gene_count}"
             : params.extend_mode == 'bp'   ? "--extend_mode bp --extend_bp ${params.extend_bp}"
             : '--extend_mode none'
    def nulls = params.num_null_gene > 0 ? "--use_null_genes --num_null_gene ${params.num_null_gene}"
                                         : '--no-use_null_genes'
    """
    run_dbcan cgc_finder \
        --output_dir ${annot_dir} \
        --min_cluster_genes ${params.min_cluster_genes} \
        --min_core_cazyme ${params.min_core_cazyme} \
        ${nulls} \
        ${add} \
        ${ext} \
        ${dist}

    cp ${annot_dir}/cgc_standard_out.tsv ${meta.id}.cgc_standard_out.tsv
    cp ${annot_dir}/cgc.gff              ${meta.id}.cgc.gff
    if [ -s ${annot_dir}/total_cgc_info.tsv ]; then
        cp ${annot_dir}/total_cgc_info.tsv ${meta.id}.total_cgc_info.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(python -c "import importlib.metadata as m; print(m.version('dbcan'))")
    END_VERSIONS
    """
}
