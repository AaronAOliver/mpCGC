/*
 * MINE - green dataflow of Figure 1.
 *
 *   CGC catalog
 *      -> enzyme colocalisation counting              (COLOC_NETWORK)
 *      -> colocalisation network (json)
 *   per-family protein fasta
 *      -> DIAMOND all-vs-all sequence similarity      (DIAMOND_ALLVSALL)
 *      -> Leiden communities on the SSN               (LEIDEN_COMMUNITIES)
 *      -> MUSCLE alignment                            (MUSCLE_ALIGN)
 *      -> FastTree phylogeny (nwk)                    (FASTTREE)
 *      -> ggtree figure                               (GGTREE_FIGURE)
 */

include { COLOC_NETWORK      } from '../modules/local/coloc_network'
include { DIAMOND_ALLVSALL   } from '../modules/local/diamond_allvsall'
include { LEIDEN_COMMUNITIES } from '../modules/local/leiden_communities'
include { MUSCLE_ALIGN       } from '../modules/local/muscle_align'
include { FASTTREE           } from '../modules/local/fasttree'
include { GGTREE_FIGURE      } from '../modules/local/ggtree_figure'

workflow MINE {

    take:
    catalog      // value: all_cgc_catalog.tsv
    proteins     // value: list of per-family faa
    protein_map  // value: protein -> genome map from CGC_PROTEINS
    metadata     // value: per-genome metadata TSV carrying taxonomy

    main:
    ch_versions = channel.empty()

    // Both are single-item channels shared by every family. As queue channels
    // they would be consumed by the first task and the remaining families would
    // be dropped without an error, so make them value channels.
    def protein_map_v = protein_map.first()
    def metadata_v    = metadata.first()

    // ---- enzyme colocalisation network --------------------------------
    COLOC_NETWORK(catalog)
    ch_versions = ch_versions.mix(COLOC_NETWORK.out.versions)

    // ---- per-family sequence similarity networks + phylogenies --------
    def wanted = params.mine_families ? params.mine_families.tokenize(',')*.trim() as Set : null

    def families = proteins
        .flatten()
        .map { faa -> tuple([ id: faa.simpleName ], faa) }
        .filter { meta, _faa -> !wanted || wanted.contains(meta.id) }
        .filter { _meta, faa -> faa.countFasta() >= params.min_family_seqs }

    DIAMOND_ALLVSALL(families)
    ch_versions = ch_versions.mix(DIAMOND_ALLVSALL.out.versions.first())

    LEIDEN_COMMUNITIES(DIAMOND_ALLVSALL.out.hits)
    ch_versions = ch_versions.mix(LEIDEN_COMMUNITIES.out.versions.first())

    if (params.make_trees) {
        MUSCLE_ALIGN(families)
        ch_versions = ch_versions.mix(MUSCLE_ALIGN.out.versions.first())

        FASTTREE(MUSCLE_ALIGN.out.alignment)
        ch_versions = ch_versions.mix(FASTTREE.out.versions.first())

        if (params.make_tree_figures) {
            GGTREE_FIGURE(
                FASTTREE.out.tree.join(LEIDEN_COMMUNITIES.out.communities),
                protein_map_v,
                metadata_v
            )
            ch_versions = ch_versions.mix(GGTREE_FIGURE.out.versions.first())
        }
    }

    emit:
    coloc       = COLOC_NETWORK.out.network
    communities = LEIDEN_COMMUNITIES.out.communities
    trees       = params.make_trees ? FASTTREE.out.tree : channel.empty()
    versions    = ch_versions
}
