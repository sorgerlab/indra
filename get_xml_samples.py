from indra.util import unzip_string
from indra.db import get_primary_db
import pickle
import zlib

db = get_primary_db()

tc_list = db.filter_query(db.TextContent, db.TextContent.text_type == 'fulltext').limit(50).all()

print('type(tc_list[0]) = ',     type(tc_list[0])   )
print('type(tc_list[0].content) = ', type(tc_list[0].content)   )
print('len(tc_list) = ', len(tc_list))

xml_list = [zlib.decompress(tc.content, 16+zlib.MAX_WBITS) for tc in tc_list]

pickle.dump(xml_list, open('xml_samples.p', 'wb'))

