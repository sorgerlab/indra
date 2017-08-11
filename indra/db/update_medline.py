import zlib
import logging
from io import BytesIO
import xml.etree.ElementTree as ET
from indra import db
from indra.literature import pubmed_client
from indra.util import UnicodeXMLTreeBuilder as UTB
from indra.db.update_content import _get_ftp_connection, ftp_blocksize

logger = logging.getLogger('update_medline')

def download_baseline():
    # Get list of .xml.gz files
    # Unzip each one
    ftp = _get_ftp_connection('/pubmed/baseline')
    # Get the list of .xml.tar.gz files
    xml_files = [f[0] for f in ftp.mlsd() if f[0].endswith('.xml.gz')]
    for xml_file in xml_files[:1]:
        # Download the Gzipped content into a BytesIO
        gzf_bytes = BytesIO()
        print("Downloading %s" % xml_file)
        ftp.retrbinary('RETR %s' % xml_file,
                                callback=lambda b: gzf_bytes.write(b),
                                blocksize=ftp_blocksize)
        # Unzip the BytesIO into uncompressed bytes
        xml_bytes = zlib.decompress(gzf_bytes.getvalue(), 16+zlib.MAX_WBITS)
        # Convert the bytes into an XML ElementTree
        tree = ET.XML(xml_bytes, parser=UTB())
    return tree

def update_text_refs():
    pass

if __name__ == '__main__':
    tree = download_baseline()
    """
    with open('medline17n0001.xml', 'rb') as f:
        content_bytes = f.read()
    tree = ET.XML(content_bytes, parser=UTB())
    article_info = pubmed_client.get_metadata_from_xml_tree(
                            tree, get_abstracts=True, prepend_title=False)
    """
