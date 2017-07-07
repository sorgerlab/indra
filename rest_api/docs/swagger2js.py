from collections import OrderedDict
import json
import os

filename = [x for x in os.listdir() if '.json' in x][0]

# with open(filename) as file:
#     swagger = json.load(file)
swagger = json.load(open(filename), object_pairs_hook=OrderedDict)

swagger['host'] = 'localhost:8080'

js_str = 'var api_spec = '
js_str += json.dumps(swagger, indent=2, sort_keys=False)

with open('api_spec.js', 'w') as file:
    file.write(js_str)
