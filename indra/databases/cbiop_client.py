import requests
import pandas
import logging
# Python3
try:
    from io import StringIO
# Python2
except ImportError:
    from StringIO import StringIO

logger = logging.getLogger('cbio')
cbio_url = 'http://www.cbioportal.org/webservice.do'


def send_request(data, skiprows=0):
    '''
    Sends a web service requrest to the cBio portal with arguments given in
    the dictionary data and returns a Pandas data frame on success.

    Parameters
    ----------
    data : dict
        A dict of parameters for the cBioPortal query
    skiprows : int
        Number of rows to skip when reading dataframe. This is useful to align
        headers

    Returns
    -------
    df : Pandas DataFrame
        return the response from cBioPortal as a Pandas DataFrame
    '''
    res = requests.get(cbio_url, params=data)
    status = res.status_code
    if status == 200:
        csv_StringIO = StringIO(res.text)
        df = pandas.read_csv(csv_StringIO, sep='\t', skiprows=skiprows)
        return df
    else:
        logger.error('Request returned with code %d' % status)
