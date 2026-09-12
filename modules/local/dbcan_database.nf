process DBCAN_DATABASE {
    label 'download'
    label 'dbcan'
    tag   'dbCAN reference databases'

    conda     "${moduleDir}/../../env/dbcan.yml"
    container 'quay.io/biocontainers/dbcan:5.2.9--pyhdfd78af_0'

    output:
    path 'dbcan_db'       , emit: db
    path 'versions.yml'   , emit: versions

    script:
    def s3 = params.db_from_s3 ? '--aws_s3' : ''
    """
    mkdir -p dbcan_db
    run_dbcan database \
        --db_dir dbcan_db \
        --cgc \
        --no-overwrite \
        --retries 5 \
        --timeout 120 \
        ${s3}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(python -c "import importlib.metadata as m; print(m.version('dbcan'))")
    END_VERSIONS
    """
}
