# INDRA Database RESTful API

One feature of INDRA is the ability to create and maintain a database of INDRA's pre-assembled Statements. Because these statements are a valuable resource in and of themselves, a web api was developed to access that knowledge.

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

As a web api, there are many ways you can interact with it. Virtually every language has tools to send `get` requests to a web API. One limitation to note: access directly via a browser can be complicated by the use of API keys, which we use in our own implementation of the API. Below, we go over a couple specific examples of getting statements using the web API, namely through `curl` and the `requests` package for `python`. In the examples, let's assume the path to the web API is `https://host.of.api.com/`, and that the API key is `12345`.

### Curl:
`curl` is a very common buildin console command in unix and unix-like OS's, making it a convenient all-purpose tool for making calls to this web api. Using `curl`, a query for "MEK phophorylates ERK" would look like
```bash
curl -i -H "x-api-key:12345" -X GET "https://host.of.api.com/statements/?subject=MAP2K1&object=MAPK1&type=phosphorylation"
```
which will return a JSON list of Statements. Note that if there is no API key, you can simply remove `-H "x-api-key:12345"` from the command.

If we instead want to know if there is any support for the proposition there is some kind of interaction between SMURF2 and JNK, your `curl` query would be:
```bash
curl -i -H "x-api-key:12345" -X GET "https://host.of.api.com/statements/?agent=SMURF2&agent=JNK"
```
You would find that indeed there are many types of interaction.

### Python:
Python is probably the most common and rewarding way to use this web API because you then have the full power of INDRA readily available. Of course, if you have access to the actual database from which the Statements are retrieved, using a web API would be circuitous. However, assuming the database is not available to you, using the web api via the `requests` package is almost just as easy. If again we want Statements that are relevant for "MEK phosphorylates ERK", you can get the Statement JSONs as follows:
```python
import json
import requests
resp = requests.get('https://host.of.api.com/statements/', headers={'x-api-key': '12345'}, params={'subject': 'MAP2K1', 'object': 'MAPK1', 'type': 'phosphorylation'})
stmts_json = json.loads(resp.content.decode('utf-8'))
```
You can then get INDRA statements by simply running `stmts_from_json` from `indra.statements`:
```python
from indra.statements import stmts_from_json
stmts = stmts_from_json(stmts_json)
```
As with the `curl` example, if there is no API key required, you simply omit that keyword argument `headers={'x-api-key': '12345'}` from the example above. For those familiar with useing preassembled INDRA Statements, note that the `supports` and `supported_by` lists of the python Statement objects will have many instances of `Unresolved` Statements, which are place-holders with nothing but the UUIDs of the referenced Statements. This is because we can only resolve UUIDs to statements that are included in the original query. Future implementations may alter this behavior.

Now suppose you want to repeat the query for interactions between SMURF2 and JNK shown above using `curl`. You have to be a bit careful here, because the parameters included two agents, which would not work with a Python `dict` as was used in our first Python example, because the keys to Python dicts are of course unique. Thus, you must in this case pass bytes to the `params` keyword argument:
```python
resp = requests.get('https://host.of.api.com/statements/', headers={'x-api-key': '12345'}, params=u'agent=JNK&agent=SMURF2')
```
From you you can continue just as you did in our first Python example to get the INDRA Statement objects from the JSON in the response's content.

### Browser:
You can only use a browser if there is no API key required to use the API. Provded there is no API key though, you can simply enter the link: `https://host.of.api.com/statements/?subject=MAP2K1&object=MAPK1` into your browser bar, and you will get a web page full of the JSON serialized Statements. Future implementations may allow for a more viewer friendly format to be display when a browser is used.
