/*
 * CGC diversity. Sample-based (incidence) rarefaction and extrapolation of
 * distinct CGC fingerprints per lineage following the iNEXT framework, with
 * Chao2 asymptotic richness and analytic log-transformed confidence intervals.
 * Genomes are the sampling units; a fingerprint is "detected" in a genome if
 * that genome encodes at least one CGC with that fingerprint.
 */
process CGC_RICHNESS {
    label 'medium'
    label 'pydata'
    tag   'CGC diversity'

    conda     "${moduleDir}/../../env/pydata.yml"
    container 'quay.io/biocontainers/scipy:1.13.0'

    input:
    path fingerprints

    output:
    path 'rarefaction_curves.tsv'    , emit: curves
    path 'richness_asymptotes.tsv'   , emit: asymptotes
    path 'rarefaction.{png,svg}'     , emit: figure, optional: true
    path 'versions.yml'              , emit: versions

    script:
    def ci = params.rarefaction_ci ? "--ci --nboot ${params.rarefaction_nboot}" : ''
    """
    mpcgc_richness.py \
        --fingerprints ${fingerprints} \
        --knots ${params.rarefaction_knots} \
        --extrapolate ${params.rarefaction_extrapolate} \
        --min-genomes ${params.min_genomes} \
        --out-curves rarefaction_curves.tsv \
        --out-asymptotes richness_asymptotes.tsv \
        --figure rarefaction \
        ${ci}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
