import urllib, urllib2
import re

trips_url = 'http://trips.ihmc.us/parser/cgi/drum'

def send_query(text, query_args={}):
    query_args['input'] = text
    data = urllib.urlencode(query_args)
    req = urllib2.Request(trips_url,data)
    res = urllib2.urlopen(req)
    html = res.read()
    return html

def get_xml(html):
    ekb = re.findall(r'<ekb .*?>(.*?)</ekb>',html,re.MULTILINE|re.DOTALL)
    if ekb:
        events_terms = ekb[0]
    else:
        events_terms = ''
    header = '<?xml version="1.0" encoding="utf-8" standalone="yes"?><extractions>'
    footer = '</extractions>'
    return header + events_terms + footer

def save_xml(xml,file_name):
    try:
        fh = open(file_name,'wt')
    except IOError:
        print 'Could not open %s for writing.' % file_name
        return
    fh.write(xml)
    fh.close()

if __name__ == '__main__':
    text = 'Active BRAF phosphorylates MEK1 at Ser222.'
    html = send_query(text)
    xml = get_xml(html)
    save_xml(xml,'braf_test.xml')

