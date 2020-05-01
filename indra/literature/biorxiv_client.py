import json
import logging
import requests


logger = logging.getLogger('biorxiv_client')


covid_json_url = 'https://connect.biorxiv.org/relate/' \
                  'collection_json.php?grp=181'


def get_covid_dois():
    """Get current list of Covid-19 relevant DOIs from biorxiv and medrxiv.

    Browser link at https://connect.biorxiv.org/relate/content/181
    """
    res = requests.get(covid_json_url)
    if res.status_code != 200:
        logger.error('Error querying biorxiv Covid-19 collection: '
                     f'status code {res.status_code}')
        return None
    res_json = res.json()
    dois = []
    for pub in res_json['rels']:
        dois.append(pub['rel_doi'])
    return dois


def get_pdf_from_doi(doi):
    pdf_url = 'https://www.biorxiv.org/content/{doi}.full.pdf'
    res = request.get(pdf_url)


if __name__ == '__main__':
    res = get_covid_dois()

