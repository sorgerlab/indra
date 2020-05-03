import json
import logging
import requests


logger = logging.getLogger(__name__)


# Browser link at https://connect.biorxiv.org/relate/content/181
collection_url = 'https://connect.biorxiv.org/relate/collection_json.php?grp='
covid19_collection_id = '181'
bio_content_url = 'https://www.biorxiv.org/content/'
med_content_url = 'https://www.medrxiv.org/content/'


format_map = {
    'xml': 'source.xml',
    'txt': 'full',
    'pdf': 'full.pdf',
    'abstract': ''
}


def get_collection_dois(collection_id):
    """Get list of DOIs from a biorxiv/medrxiv collection.

    Parameters
    ----------
    collection_id : str
        The identifier of the collection to fetch.

    Returns
    -------
    list of str
        A list of DOIs in the given collection.
    """
    res = requests.get(collection_url + collection_id)
    res.raise_for_status()
    dois = [pub['rel_doi'] for pub in res.json()['rels']]
    return dois


def get_content(doi, bio_or_med='bio', format='xml'):
    url = bio_content_url if bio_or_med == 'bio' else med_content_url
    url += ('%s.%s' % (doi, format))
    res = requests.get(url)
    return res.text


if __name__ == '__main__':
    covid19_dois = get_collection_dois(covid19_collection_id)
    get_content(covid19_dois[0])


