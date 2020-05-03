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


def get_pdf_xml_url_base(content):
    """Return base URL to PDF/XML based on the content of the landing page.

    Parameters
    ----------
    content : str
        The content of the landing page for an rxiv paper.

    Returns
    -------
    str or None
        The base URL if available, otherwise None.
    """
    match = re.match('(?:.*)"citation_pdf_url" content="([^"]+).full.pdf"',
                     content, re.S)
    if match:
        return match.groups()[0]
    return None


def get_text_url_base(content):
    """Return base URL to full text based on the content of the landing page.

    Parameters
    ----------
    content : str
        The content of the landing page for an rxiv paper.

    Returns
    -------
    str or None
        The base URL if available, otherwise None.
    """
    match = re.match('(?:.*)"citation_html_url" content="([^"]+).full"',
                     content, re.S)
    if match:
        return match.groups()[0]
    return None


def get_formats(pub):
    """Return formats available for a publication JSON.

    Parameters
    ----------
    pub : dict
        The JSON dict description a publication.

    Returns
    -------
    dict
        A dict with available formats as its keys (abstract, pdf, xml, txt)
        and either the content (in case of abstract) or the URL
        (in case of pdf, xml, txt) as the value.
    """
    formats = {}
    if 'rel_abs' in pub:
        formats['abstract'] = pub['rel_abs']
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
    if pdf_xml_url_base:
        formats['pdf'] = pdf_xml_url_base + '.full.pdf'
        formats['xml'] = pdf_xml_url_base + '.source.xml'
    text_url_base = get_text_url_base(landing_page_res.text)
    if text_url_base:
        formats['txt'] = text_url_base + 'txt'
    return formats


def get_content_from_pub_json(pub, format):
    """Get text content based on a given format from a publication JSON.

    In the case of abstract, the content is returned
    from the JSON directly. For pdf, the content is returned as bytes
    that can be dumped into a file. For txt and xml, the text is processed
    out of either the raw XML or text content that rxiv provides.

    Parameters
    ----------
    pub : dict
        The JSON dict description a publication.
    format : str
        The format, if available, via which the content should be
        obtained.
    """
    if format == 'abstract':
        return pub.get('rel_abstract')

    formats = get_formats(pub)
    if format not in formats:
        logger.warning('Content not available in format %s' % format)
        return None

    # If we're looking for an abstract, that is directly accessible
    # in the pub JSON so we can just return it
    if format == 'abstract':
        return formats.get('abstract')
    # For PDFs we return the result in bytes that can then be dumped
    # into a file.
    elif format == 'pdf':
        return requests.get(formats[format]).content
    # For xml and text, we return the result as str
    elif format == 'xml':
        return get_text_from_rxiv_xml(requests.get(formats[format]).text)
    elif format == 'txt':
        return get_text_from_rxiv_text(requests.get(formats[format]).text)


def get_text_from_rxiv_xml(rxiv_xml):
    """Return clean text from the raw rxiv xml content.

    Paramteres
    ----------
    rxiv_xml : str
        The content of the rxiv full xml as obtained from the web.

    Returns
    -------
    str
        The text content stripped out from the raw full xml.
    """
    # FIXME: this is a very naive initial solution, we should instead
    # traverse the XML structure properly to get the content.
    text = re.sub('<.*?>', '', rxiv_xml)
    return text


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
    import os
    import json
    fname = 'covid19_pubs.json'
    if os.path.exists(fname):
        with open(fname, 'r') as fh:
            covid19_pubs = json.load(fh)
    else:
        covid19_pubs = get_collection_pubs(covid19_collection_id)
        with open(fname, 'w') as fh:
            json.dump(covid19_pubs, fh)
    contents = {}
    for pub in covid19_pubs:
        doi = pub['rel_doi']
        formats = get_formats(pub)
        if 'txt' in formats:
            print('Getting text for %s' % doi)
            txt = get_content_from_pub_json(pub, 'txt')
        elif 'xml' in formats:
            print('Getting xml for %s' % doi)
            txt = get_content_from_pub_json(pub, 'xml')
        else:
            print('Getting abstract for %s' % doi)
            txt = get_content_from_pub_json(pub, 'abstract')
        contents[doi] = txt
    with open('covid19_contents', 'w') as fh:
        json.dump(contents, fh)

