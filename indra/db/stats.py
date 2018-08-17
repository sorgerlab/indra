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


def _report_groups(db, count_obj, group_obj, fname, *filters):
    """Report on the number of rows by group."""
    q = db.session.query(group_obj, func.count(count_obj))
    if filters:
        q = q.filter(*filters)
    content_by_group = (q.distinct()
                        .group_by(group_obj)
                        .all())
    broad_format = "{table} by {group}:\n    {values}"
    value_strs = ['%s: %d' % (s, n) for s, n in content_by_group]
    value_str = '\n    '.join(value_strs)
    report_str = broad_format.format(group=str(group_obj), values=value_str,
                                     table=str(count_obj))
    __report_stat(report_str, fname)
    return {s: n for s, n in content_by_group}


def get_text_ref_stats(fname=None, db=None):
    if db is None:
        db = get_primary_db()
    tc_rdng_link = db.TextContent.id == db.Reading.text_content_id
    __report_stat("Text ref statistics:", fname)
    __report_stat("--------------------", fname)
    total_refs = db.count(db.TextRef)
    __report_stat('Total number of text refs: %d' % total_refs, fname)
    refs_with_content = db.count(db.TextContent.text_ref_id)
    __report_stat('Total number of refs with content: %d' % refs_with_content,
                  fname)
    refs_by_type = _report_groups(db, db.TextContent.text_ref_id,
                                  db.TextContent.text_type, fname)
    __report_stat(('Number of refs with only abstract: %d'
                   % (refs_with_content-refs_by_type['fulltext'])), fname)
    refs_with_reading = db.count(db.TextContent.text_ref_id,
                                 tc_rdng_link)
    __report_stat('Number of refs that have been read: %d' % refs_with_reading,
                  fname)
    _report_groups(db, db.TextContent.text_ref_id, db.TextContent.text_type,
                   fname, tc_rdng_link)
    return


def get_text_content_stats(fname=None, db=None):
    if db is None:
        db = get_primary_db()
    tc_rdng_link = db.TextContent.id == db.Reading.text_content_id
    __report_stat("\nText Content statistics:", fname)
    __report_stat('------------------------', fname)
    total_content = db.count(db.TextContent)
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
    content_read = db.count(db.Reading.text_content_id)
    __report_stat("Total content read: %d" % content_read, fname)
    fulltext_content = db.count(db.TextContent,
                                db.TextContent.text_type == 'fulltext')
    __report_stat("Number of fulltext entries: %d" % fulltext_content, fname)
    fulltext_read = db.count(db.TextContent,
                             db.TextContent.text_type == 'fulltext',
                             tc_rdng_link)
    __report_stat("Number of fulltext entries read: %d" % fulltext_read, fname)
    _report_groups(db, db.TextContent.id, db.TextContent.source, fname)
    _report_groups(db, db.TextContent.id, db.TextContent.source, fname,
                   tc_rdng_link)
    return


def get_readings_stats(fname=None, db=None):
    if db is None:
        db = get_primary_db()

    __report_stat('\nReading statistics:', fname)
    __report_stat('-------------------', fname)
    total_readings = db.count(db.Reading)
    __report_stat('Total number or readings: %d' % total_readings, fname)
    # There may be a way to do this more neatly with a group_by clause, however
    # the naive way of doing it leaves us with a miscount due to indistinct.
    reader_versions = (db.session.query(db.Reading.reader_version)
                       .distinct().all())
    sources = db.session.query(db.TextContent.source).distinct().all()
    stats = ''
    for rv, in reader_versions:
        for src, in sources:
            cnt = db.count(
                db.Reading,
                db.TextContent.id == db.Reading.text_content_id,
                db.TextContent.source == src,
                db.Reading.reader_version == rv
                )
            stats += '    Readings by %s from %s: %d\n' % (rv, src, cnt)
    __report_stat("Readings by reader version and content source:\n%s" % stats,
                  fname)
    return


def get_statements_stats(fname=None, db=None, indra_version=None):
    if db is None:
        db = get_primary_db()
    tc_rdng_link = db.TextContent.id == db.Reading.text_content_id
    stmt_rdng_link = db.Reading.id == db.RawStatements.reader_ref

    __report_stat('\nStatement Statistics:', fname)
    __report_stat('---------------------', fname)
    if indra_version is not None:
        filters = [db.RawStatements.indra_version == indra_version]
    total_raw_statements = db.count(db.RawStatements, *filters)
    __report_stat("Total number of raw statements: %d" % total_raw_statements,
                  fname)
    readers = db.session.query(db.Reading.reader).distinct().all()
    sources = db.session.query(db.TextContent.source).distinct().all()
    stats = ''
    for reader, in readers:
        for src, in sources:
            cnt = db.count(db.RawStatements, stmt_rdng_link, tc_rdng_link,
                           db.Reading.reader == reader,
                           db.TextContent.source == src, *filters)
            stats += ('    Raw statements from %s reading %s: %d\n'
                      % (reader, src, cnt))
    __report_stat("Statements by reader and content source:\n%s" % stats,
                  fname)
    _report_groups(db, db.RawStatements.id, db.DBInfo.db_name, fname,
                   db.RawStatements.db_info_id == db.DBInfo.id)
    if indra_version is None:
        _report_groups(db, db.RawStatements.id, db.RawStatements.indra_version,
                       fname)
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

    task_order = ['text_ref', 'text_content', 'readings', 'statements',
                  'pa_statements']

    # Get the statistics
    if tables is None:
        for task_name in task_order:
            stat_meth = task_dict[task_name]
            stat_meth(fname, db)
    else:
        table_set = set(tables)
        for task_name in [tn for tn in task_order if tn in table_set]:
            task_dict[task_name](fname, db)

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
