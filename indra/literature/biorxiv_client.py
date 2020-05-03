import re
import logging
import requests


logger = logging.getLogger(__name__)


# Browser link at https://connect.biorxiv.org/relate/content/181
collection_url = 'https://connect.biorxiv.org/relate/collection_json.php?grp='
covid19_collection_id = '181'
bio_content_url = 'https://www.biorxiv.org/content/'
med_content_url = 'https://www.medrxiv.org/content/'


def get_collection_pubs(collection_id):
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
    return res.json()['rels']


def get_content(doi, bio_or_med='bio', format='xml'):
    url = bio_content_url if bio_or_med == 'bio' else med_content_url
    url += ('%s.%s' % (doi, format))
    res = requests.get(url)
    return res.text


def get_pdf_xml_url_base(content):
    match = re.match('(?:.*)"citation_pdf_url" content="([^"]+).full.pdf"',
                     content, re.S)
    if match:
        return match.groups()[0]
    return None


def get_text_url_base(content):
    match = re.match('(?:.*)"citation_html_url" content="([^"]+).full"',
                     content, re.S)
    if match:
        return match.groups()[0]
    return None


def get_pub_content(pub, format):
    # If we're looking for an abstract, that is directly accessible
    # in the pub JSON so we can just return it
    if format == 'abstract':
        return pub.get('rel_abs')

    # The publication JSON does not contain enough information generally
    # to identify the URL for the various formats. Therefore we have to
    # load the landing page for the article and parse out various URLs
    # to reliably get to the desired content.
    landing_page_res = requests.get(pub['rel_link'])

    # The URL for the full PDF and XML is often different in format than
    # the rel_site URL so we need to get the link to it from the content
    # of the landing page. The XML URL doesn't explicitly appear in the
    # page content therefore we work with the citation_pdf_url and get
    # URLs for both the PDF and the XML.
    pdf_xml_url_base = get_pdf_xml_url_base(landing_page_res.text)
    text_url_base = get_text_url_base(landing_page_res.text)
    # For PDFs we return the result in bytes that can then be dumped
    # into a file.
    if format == 'pdf':
        if not pdf_xml_url_base:
            logger.warning('Could not get PDF URL for this content.')
            return None
        url = pdf_xml_url_base + '.full.pdf'
        return requests.get(url).content
    # For xml and text, we return the result as str
    elif format == 'xml':
        if not pdf_xml_url_base:
            logger.warning('Could not get XML URL for this content.')
            return None
        url = pdf_xml_url_base + 'source.xml'
        return requests.get(url).text
    elif format == 'txt':
        url = text_url_base + 'txt'
        return requests.get(url).text
    else:
        logger.warning('Unknown format: %s' % format)
        return None


def get_text_from_rxiv_text(rxiv_text):
    """Return clean text from the raw rxiv text content.

    This function parses out the title, headings and subheadings, and
    the content of sections under headings/subheadings.
    It filters out some irrelevant content e.g., references and
    footnotes.

    Paramteres
    ----------
    rxiv_text : str
        The content of the rxiv full text as obtained from the web.

    Returns
    -------
    str
        The text content stripped out from the raw full text.
    """
    lines = [line.strip() for line in rxiv_text.split('\n') if line.strip()]
    current_section = 'title'
    text = lines[0] + '\n'
    line_idx = 1
    skip_section = {'References', 'Footnotes', 'Acknowledgements',
                    'Supplementary Figures', 'Declaration of Interests',
                    'Author Contributions', 'Code and data availability'}
    for line in lines[line_idx:]:
        line_idx += 1
        match = re.match('## (.+)', line)
        if match:
            current_section = match.groups()[0]
            break
    while line_idx < len(lines):
        for line in lines[line_idx:]:
            line_idx += 1
            match_heading = re.match('## (.+)', line)
            match_subheading = re.match('### (.+)', line)
            if match_heading:
                current_section = match_heading.groups()[0]
                break
            elif current_section in skip_section:
                continue
            elif match_subheading:
                text += (match_subheading.groups()[0] + '\n')
            else:
                text += (line + '\n')
    return text


if __name__ == '__main__':
    covid19_pubs = get_collection_pubs(covid19_collection_id)
    abs = get_pub_content(covid19_pubs[0], 'abstract')
    pdf = get_pub_content(covid19_pubs[0], 'pdf')


