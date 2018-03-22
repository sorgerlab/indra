# INDRA Database RESTful API

One feature of INDRA is the ability to create and maintain a database of INDRA's pre-assembled statements. Because these statements are a valuable resource in and of themselves, a web api was developed to access that knowledge.

## Implementation

Flask was used to write the python api, as you can see in `api.py`. The means of hosting this api are left to the user. We have had greate success with [Zappa](https://github.com/Miserlou/Zappa) and AWS, and would recommend it for a quick and efficient way to get the API up and running.

## Client Usage

Currently, the only usage for this api is to retrieve preassembled statements from the database. The query parameters available are:

### Agent specifications:
- **`subject`**, **`object`** : Give the HGNC gene symbol of the subject or object of your desired set of statements. For example, if you are looking for statements where "MEK phosphorylates something", you could include `subject=MAP2K1` in your query, and if you want examples of "MEK phosporylates ERK", you could also add `object=MAPK1`. **Note**: only one of each of `subject` and `object` will be accepted per query.
- **`agent`** : If you don't care whether an agent is the subject or object, or that doesn't apply to the type of statements you want, such as Complexes, you can use the `agent` query. For example, if you simply want statements that involve SMAD2, you would add `agent=SMAD2` to your query. **Note**: You can include as many `agent` queries as you like, however you will only get statements that include all agents you query, in addition to those queried for `subject` and `object`.

### Other Options:
- **`type`** : You can also specify what type of statement you are looking for, e.g. Phosphorylation, Activation, Complex, etc. In the case above where you wanted examples of "MEK phosphorylates ERK", you would need to include `type=Phosphorylation` in your query. Note that this field is not case sensitive, so `type=pHoSpHoRyLaTiOn` would have the same result as above.

## Examples

As a web api, there are many ways you can interact with it. The key exception is that access directly via a browser can be complicated by the use of api keys, which we use in our own implementation which we host on AWS. Here are a couple examples of ways getting statements from the web api using `curl` and the `requests` package for `python`. In the examples, let's assume the path to the web api is `https://host.of.api.com/`, and that the api key is `12345`.

### Curl:
`curl` is a very common buildin console command in unix and unix-like OS's, making it a convenient all-purpose tool for making calls to this web api. Using `curl`, a query for "MEK phophorylates ERK" would look like
```bash
curl -i -H "x-api-key:12345" -X GET "https://host.of.api.com/statements/?subject=MAP2K1&object=MAPK1"
```
which will return a json list of statements. Note that if there is no API key, you can simply remove `-H "x-api-key:12345"` from the command.

### Python:
Python is probably the simplest way to use this web api, although if you have access to the entire database itself, this would be circuitous. However, assuming the database is not available directly, you can use the web api via the `requests` package to get the json serialization of the statements:
```python
import json
import requests
resp = requests.get('https://host.of.api.com/statements/', headers={'x-api-key': '12345'}, params={'subject': 'MAP2K1', 'object': 'MAPK1'})
stmts_json = json(resp.content.decode('utf-8'))
```
If you then want to get the actual statements, you simply run `stmts_from_json`:
```python
from indra.statements import stmts_from_json
stmts = stmts_from_json(stmts_json)
```
Again, if there is no API key require, you simply omit that keyword argument.

### Browser:
You can only use a browser if there is no API key required to use the API. Provded there is no API key though, you can simply enter the link: `https://host.of.api.com/statements/?subject=MAP2K1&object=MAPK1` into your browser bar, and you will get a web page full of the JSON serialized statements.
