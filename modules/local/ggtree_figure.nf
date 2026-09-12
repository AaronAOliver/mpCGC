/*
 * ggtree rendering of a family phylogeny.
 *
 * Coloured by host lineage when --tree_color_by is 'taxonomy' (the default) and a
 * taxonomy table is available, otherwise by Leiden community. Tree tips are
 * protein identifiers, so the protein map from CGC_PROTEINS supplies the
 * protein -> genome join that taxonomy colouring needs.
 */
process GGTREE_FIGURE {
    label 'low'
    label 'rggtree'
    tag   "${meta.id}"

    conda     "${moduleDir}/../../env/rggtree.yml"
    container 'quay.io/biocontainers/bioconductor-ggtree:3.12.0--r43hdfd78af_0'

    input:
    tuple val(meta), path(nwk), path(communities)
    path  protein_map
    path  taxonomy

    output:
    // svg is absent when the R environment cannot write one; png is always made
    tuple val(meta), path("${meta.id}.tree.png")      , emit: figure
    tuple val(meta), path("${meta.id}.tree.svg")      , optional: true, emit: figure_svg
    path 'versions.yml'                               , emit: versions

    script:
    """
    mpcgc_tree_figure.R \
        --tree ${nwk} \
        --annotation ${communities} \
        --protein-map ${protein_map} \
        --taxonomy ${taxonomy} \
        --color-by ${params.tree_color_by} \
        --rank ${params.tree_color_rank} \
        --max-groups ${params.tree_max_groups} \
        --family ${meta.id} \
        --out-prefix ${meta.id}.tree

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript -e 'cat(strsplit(R.version.string," ")[[1]][3])')
        ggtree: \$(Rscript -e 'cat(as.character(packageVersion("ggtree")))')
    END_VERSIONS
    """
}
