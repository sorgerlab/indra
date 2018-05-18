from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from datetime import datetime

import boto3
from sqlalchemy import func

from indra.db.util import get_primary_db


def __report_stat(report_str, fname=None, do_print=True):
    if do_print:
        print(report_str)
    if fname is not None:
        with open(fname, 'a+') as f:
            f.write(report_str + '\n')
    return


def get_text_ref_stats(fname=None, db=None):
    if db is None:
        db = get_primary_db()
    tr_tc_link = db.TextRef.id == db.TextContent.text_ref_id
    tc_rdng_link = db.TextContent.id == db.Readings.text_content_id
    __report_stat("Text ref statistics:", fname)
    __report_stat("--------------------", fname)
    tr_q = db.filter_query(db.TextRef)
    total_refs = tr_q.count()
    __report_stat('Total number of text refs: %d' % total_refs, fname)
    tr_w_cont_q = tr_q.filter(tr_tc_link)
    refs_with_content = tr_w_cont_q.distinct().count()
    __report_stat('Total number of refs with content: %d' % refs_with_content,
                  fname)
    tr_w_fulltext_q = tr_w_cont_q.filter(db.TextContent.text_type == 'fulltext')
    refs_with_fulltext = tr_w_fulltext_q.distinct().count()
    __report_stat('Number of refs with fulltext: %d' % refs_with_fulltext,
                  fname)
    tr_w_abstract_q = tr_w_cont_q.filter(db.TextContent.text_type == 'abstract')
    refs_with_abstract = tr_w_abstract_q.distinct().count()
    __report_stat('Number of refs with abstract: %d' % refs_with_abstract,
                  fname)
    __report_stat(('Number of refs with only abstract: %d'
                   % (refs_with_content-refs_with_fulltext)), fname)
    tr_w_read_content_q = tr_w_cont_q.filter(tc_rdng_link)
    refs_with_reading = tr_w_read_content_q.distinct().count()
    __report_stat('Number of refs that have been read: %d' % refs_with_reading,
                  fname)
    tr_w_fulltext_read_q = tr_w_fulltext_q.filter(tc_rdng_link)
    refs_with_fulltext_read = tr_w_fulltext_read_q.distinct().count()
    __report_stat(('Number of refs with fulltext read: %d'
                   % refs_with_fulltext_read), fname)
    return


def get_text_content_stats(fname=None, db=None):
    if db is None:
        db = get_primary_db()
    tc_rdng_link = db.TextContent.id == db.Readings.text_content_id
    __report_stat("\nText Content statistics:", fname)
    __report_stat('------------------------', fname)
    tc_q = db.filter_query(db.TextContent)
    total_content = tc_q.count()
    __report_stat("Total number of text content entries: %d" % total_content)
    latest_updates = (db.session.query(db.Updates.source,
                                       func.max(db.Updates.datetime))
                      .group_by(db.Updates.source)
                      .all())
    __report_stat(("Latest updates:\n    %s"
                   % '\n    '.join(['%s: %s' % (s, d)
                                    for s, d in latest_updates])),
                  fname
                  )
    tc_w_reading_q = tc_q.filter(tc_rdng_link)
    content_read = tc_w_reading_q.distinct().count()
    __report_stat("Total content read: %d" % content_read, fname)
    tc_fulltext_q = tc_q.filter(db.TextContent.text_type == 'fulltext')
    fulltext_content = tc_fulltext_q.distinct().count()
    __report_stat("Number of fulltext entries: %d" % fulltext_content, fname)
    tc_fulltext_read_q = tc_fulltext_q.filter(tc_rdng_link)
    fulltext_read = tc_fulltext_read_q.distinct().count()
    __report_stat("Number of fulltext entries read: %d" % fulltext_read, fname)
    content_by_source = (db.session.query(db.TextContent.source,
                                          func.count(db.TextContent.id))
                         .distinct()
                         .group_by(db.TextContent.source)
                         .all())
    __report_stat(("Content by source:\n    %s"
                   % '\n    '.join(['%s: %d' % (s, n)
                                    for s, n in content_by_source])),
                  fname
                  )
    content_read_by_source = (db.session.query(db.TextContent.source,
                                               func.count(db.TextContent.id))
                              .filter(tc_rdng_link)
                              .distinct()
                              .group_by(db.TextContent.source)
                              .all())
    __report_stat(("Content read by source:\n    %s"
                   % '\n    '.join(['%s: %d' % (s, n)
                                    for s, n in content_read_by_source])),
                  fname
                  )
    return


def get_readings_stats(fname=None, db=None):
    if db is None:
        db = get_primary_db()

    __report_stat('\nReading statistics:', fname)
    __report_stat('-------------------', fname)
    rdg_q = db.filter_query(db.Readings)
    __report_stat('Total number or readings: %d' % rdg_q.count(), fname)
    # There may be a way to do this more neatly with a group_by clause, however
    # the naive way of doing it leaves us with a miscount due to indistinct.
    reader_versions = (db.session.query(db.Readings.reader_version)
                       .distinct().all())
    sources = db.session.query(db.TextContent.source).distinct().all()
    stats = ''
    for rv, in reader_versions:
        for src, in sources:
            cnt = db.filter_query(
                db.Readings,
                db.TextContent.id == db.Readings.text_content_id,
                db.TextContent.source == src,
                db.Readings.reader_version == rv
            ).distinct().count()
            stats += '    Readings by %s from %s: %d\n' % (rv, src, cnt)
    __report_stat("Readings by reader version and content source:\n%s" % stats,
                  fname)
    return


def get_statements_stats(fname=None, db=None, indra_version=None):
    if db is None:
        db = get_primary_db()
    tc_rdng_link = db.TextContent.id == db.Readings.text_content_id
    stmt_rdng_link = db.Readings.id == db.Statements.reader_ref

    __report_stat('\nStatement Statistics:', fname)
    __report_stat('---------------------', fname)
    stmt_q = db.filter_query(db.Statements)
    if indra_version is not None:
        stmt_q = stmt_q.filter(db.Statements.indra_version == indra_version)
    __report_stat("Total number of statments: %d" % stmt_q.count(), fname)
    readers = db.session.query(db.Readings.reader).distinct().all()
    sources = db.session.query(db.TextContent.source).distinct().all()
    stats = ''
    for reader, in readers:
        for src, in sources:
            cnt = stmt_q.filter(
                stmt_rdng_link,
                tc_rdng_link,
                db.Readings.reader == reader,
                db.TextContent.source == src
            ).distinct().count()
            stats += ('    Statements from %s reading %s: %d\n'
                      % (reader, src, cnt))
    __report_stat("Statements by reader and content source:\n%s" % stats,
                  fname)
    if indra_version is None:
        statements_by_db_source = (
            db.session.query(db.DBInfo.db_name, func.count(db.Statements.id))
                .filter(db.Statements.db_ref == db.DBInfo.id)
                .distinct()
                .group_by(db.DBInfo.db_name)
                .all()
        )
        __report_stat(("Statements by database:\n    %s"
                       % '\n    '.join(['%s: %d' % (s, n)
                                        for s, n in statements_by_db_source])),
                      fname
                      )
        statements_by_indra_version = (
            db.session.query(db.Statements.indra_version,
                             func.count(db.Statements.id))
                .group_by(db.Statements.indra_version)
                .all()
        )
        __report_stat(("Number of statements by indra version:\n    %s"
                       % '\n    '.join(['%s: %d' % (s, n) for s, n
                                        in statements_by_indra_version])),
                      fname
                      )
    return


def get_pa_statement_stats(fname=None, db=None):
    if db is None:
        db = get_primary_db()
    __report_stat('\nStatement Statistics:', fname)
    __report_stat('---------------------', fname)
    stmt_q = db.filter_query(db.PAStatements)
    __report_stat("Total number of statments: %d" % stmt_q.count(), fname)
    statements_produced_by_indra_version = (
        db.session.query(db.PAStatements.indra_version,
                         func.count(db.PAStatements.id))
            .group_by(db.PAStatements.indra_version)
            .all()
    )
    __report_stat(("Number of statements by indra version:\n    %s"
                   % '\n    '.join(['%s: %d' % (s, n) for s, n
                                    in statements_produced_by_indra_version])),
                  fname
                  )
    return


def get_db_statistics(fname=None, db=None, tables=None):
    """Get statistics on the contents of the database"""
    if db is None:
        db = get_primary_db()

    task_dict = {
        'text_ref': get_text_ref_stats,
        'text_content': get_text_content_stats,
        'readings': get_readings_stats,
        'statements': get_statements_stats,
        'pa_statements': get_pa_statement_stats
    }

    # Get the statistics
    if tables is None:
        for stat_meth in task_dict.values():
            stat_meth(fname, db)
    else:
        for table_key in set(tables):
            task_dict[table_key](fname, db)

    return


def get_all_db_statistics_s3():
    utcnow = datetime.utcnow()
    fname = "Primary_Database_Status_Report_%s.txt" % utcnow.strftime("%Y%m%d")
    print("Creating report in: %s." % fname)

    print("\nBegin Report============\n")
    get_db_statistics(fname)
    print("\nEnd Report==============\n")

    print("Saving record to s3.")
    s3 = boto3.client('s3')
    with open(fname, 'rb') as f:
        s3.put_object(Body=f, Bucket='bigmech', Key='indra-db/reports/%s' % fname)
    return


if __name__ == '__main__':
    get_all_db_statistics_s3()
