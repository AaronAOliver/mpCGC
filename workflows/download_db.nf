/*
 * download_db - fetch and format every reference database mpCGC needs.
 *
 * dbCAN's own downloader pulls the CAZy DIAMOND database, the dbCAN and
 * dbCAN-sub HMMs, TCDB, the PRODORIC transcription-factor set, the
 * signal-transduction HMMs, SulfAtlas and MEROPS peptidases in one step.
 */

include { DBCAN_DATABASE } from '../modules/local/dbcan_database'

workflow DOWNLOAD_DB {

    main:
    DBCAN_DATABASE()

    emit:
    db       = DBCAN_DATABASE.out.db
    versions = DBCAN_DATABASE.out.versions
}
