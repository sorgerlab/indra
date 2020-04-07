import pandas

vhn_url = ('http://virhostnet.prabi.fr:9090/psicquic/webservices/current/'\
           'search/query/*')


def process_from_web():
    df = pandas.read_csv(vhn_url)
    return process_df(df)


def process_df():
    vp = VirhostnetProcessor(df)
    vp.extract_statements()
    return vp
