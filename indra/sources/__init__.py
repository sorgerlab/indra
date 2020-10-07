import json
from indra.resources import open_resource_file

with open_resource_file("source_info.json", "r") as f:
    SOURCE_INFO = json.load(f)


# This is a bit hack-y, but the names are so interchanged the alternative is
# a LOT of work.
SOURCE_INFO['trips'] = SOURCE_INFO['drum']
