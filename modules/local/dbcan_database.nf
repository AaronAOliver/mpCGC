process DBCAN_DATABASE {
    label 'download'
    label 'dbcan'
    tag   'dbCAN reference databases'

    conda     "${moduleDir}/../../env/dbcan.yml"

    output:
    path 'dbcan_db'       , emit: db
    path 'versions.yml'   , emit: versions

    script:
    def s3 = params.db_from_s3.toString().toBoolean() ? '--aws_s3' : ''
    // The download is ~8 GB and the dbCAN servers occasionally refuse a file or
    // drop a connection. Rerunning in the same directory with --no-overwrite
    // keeps every file already fetched, so a hiccup costs one file rather than
    // a fresh 8 GB task retry.
    """
    mkdir -p dbcan_db
    ok=0
    for attempt in 1 2 3 4 5; do
        if run_dbcan database \
            --db_dir dbcan_db \
            --cgc \
            --no-overwrite \
            --retries 5 \
            --timeout 120 \
            ${s3}; then
            ok=1
            break
        fi
        echo "database download incomplete (attempt \$attempt of 5); retrying in 60 s" >&2
        rm -f dbcan_db/*.part
        sleep 60
    done
    if [ "\$ok" -ne 1 ]; then
        echo "dbCAN database download failed after 5 attempts" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(python -c "import importlib.metadata as m; print(m.version('dbcan'))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dbcan_db
    echo '"${task.process}": {}' > versions.yml
    """
}
