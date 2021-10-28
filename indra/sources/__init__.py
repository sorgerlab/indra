import json
from indra.resources import open_resource_file

with open_resource_file("source_info.json", "r") as f:
    SOURCE_INFO = json.load(f)
