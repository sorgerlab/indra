# INDRA Database REST API

One feature of INDRA is the ability to create and maintain a database of
INDRA Statements. This web API allows accessing Statements in a database by
searching Statements matching a given set of query paramerters and returning
the Statements in a JSON serialized form.

You need the following information to access a running web service:
- The address of the web service (below shown with the placeholder
host.of.api.com)
- (for some implementations) An API key which needs to be sent in the header of each request to the
service, or any other credentials that are implemented.

If you want to use our implementation of the web API, you can contact us for
the path and the API key.

The service in `api.py` is implemented using the Flask Python package.
The means of hosting this api are left to the user.
We have had success with [Zappa](https://github.com/Miserlou/Zappa) and AWS,
and recommend it for a quick and efficient way to get the API up and running.

## Search parameters

The API accepts requests to retrieve serialized INDRA Statements
from a database according to various search criteria. You can specify
the agent arguments as well as the type of Statement of interest. 
The query parameters available are as follows:
- **`subject`**, **`object`**: The HGNC gene symbol of the subject or
object of the Statement.
**Note**: only one of each of `subject` and `object` will be accepted per
query.
  - Example 1: if looking for Statements where MAP2K1 is a subject
(e.g. "What does MAP2K1 phosphorylate?"), specify
`subject=MAP2K1` as a query parameter
  - Example 2: if looking for Statements where MAP2K1 is the subject and
MAPK1 is the object, add both `subject=MAP2K1` and `object=MAPK1` as
query parameters.
  - Example 3: you can specify the agent id namespace by appending
  `@<namespace>` to the agent id in the parameter, e.g.
  `subject=6871@HGNC`.
- **`agent*`**: This parameter is used if the specific role of the agent
(subject or object) is irrelevant, or the distinction doesn't apply to the
type of Statement of interest (e.g. Complex, Translocation, ActiveForm).
**Note**: You can include as many `agent*` queries as you like, however you
will only get Statements that include all agents you query, in addition to
those queried for `subject` and `object`. Furthermore, to include multiple
agents on our particular implementation, which uses the AWS API Gateway,
you must include a suffix to each agent key, such as `agent0` and `agent1`,
or else all but one agent will be stripped out. Note that you need not use
integers, you can add any suffix you like, e.g. `agentOfDestruction=TP53`
would be entirely valid.
  - Example 1: To obtain Statements that involve SMAD2 in any role, add
  `agent=SMAD2` to the query.
  - Example 2: As with `subject` and `object`, you can specify the
  namespace for an agent by appending `@<namespace>` to the agent's id, e.g.
  `agent=ERK@TEXT`.
  - Example 3: If you wanted to query multiple statements, you could include
  `agent0=MEK@FPLX` and `agent1=ERK@FPLX`. Note that the value of the
   integers has no real bearing on the ordering, and only serves to make the
    agents uniquely keyed. Thus `agent1=MEK@FPLX` and `agent0=ERK@FPLX` will
     give exactly the same result.
- **`type`**: This parameter can be used to specify what type of Statement
of interest (e.g. Phosphorylation, Activation, Complex).
  - Example: To answer the question "Does MAP2K1 phosphorylate MAPK1?"
the parameter `type=Phosphorylation` can be included in your query.
Note that this field is not case sensitive, so `type=phosphorylation` would
give the same result.

## Usage examples

The web service accepts standard GET requests, and any client that can
send such requests can be used to interact with the service. Note, however,
that access directly via a browser can be complicated by the
use of API keys which have to be included in the request header. Here we
provide usage examples with the `curl` command line tool and `python`.

In the examples, let's assume the path to the web API is
`https://host.of.api.com/`, and that the API key is `12345`.

### Curl
`curl` is a command line tool on Linux and Mac, making it a convenient tool
for making calls to this web API.

Example 1: Using `curl`, a query for "MAP2K1 phophorylates MAPK1" can be sent
as:
```bash
curl -i -H "x-api-key:12345" -X GET "https://host.of.api.com/statements/?subject=MAP2K1&object=MAPK1&type=phosphorylation"
```
which will return a JSON list of Statements. Note that if there is no API key,
you can simply remove `-H "x-api-key:12345"` from the command.

Example 2: If we instead want to know if there is any support for the
proposition there is some kind of interaction between SMURF2 and SMAD2,
your `curl` query would be:
```bash
curl -i -H "x-api-key:12345" -X GET "https://host.of.api.com/statements/?
agent=SMURF2&agent=SMAD2"
```
Or, if you are using our implementation of the REST service, you would
need to use
```bash
curl -i -H "x-api-key:12345" -X GET "https://host.of.api.com/statements/?
agent0=SMURF2&agentbondjamesbond=SMAD2"
```


### Python
Python is a convenient way to use this web API and has the important
advantage that Statements returned from the service can directly be used
by INDRA in the same environment. If again we want Statements that are
relevant for "MEK phosphorylates ERK", you can get the Statement JSONs
as follows:
```python
import requests
resp = requests.get('https://host.of.api.com/statements/',
                    headers={'x-api-key': '12345'},
                    params={'subject': 'MAP2K1',
                            'object': 'MAPK1',
                            'type': 'phosphorylation'})
stmts_json = resp.json()
```
You can also instantiate INDRA Statement objects by calling `stmts_from_json`
from `indra.statements` as :
```python
from indra.statements import stmts_from_json
stmts = stmts_from_json(stmts_json)
```
As with the `curl` example, if there is no API key required, you simply omit
that keyword argument `headers={'x-api-key': '12345'}` from the example above.
For those familiar with useing preassembled INDRA Statements, note that the
`supports` and `supported_by` lists of the python Statement objects can
have instances of `Unresolved` Statements, which are placeholders
of referenced Statements that are not included and resolved in detail in
the response.

Now suppose you want to query for interactions between SMURF2 and SMAD2 but
without specifying their specific roles.
This requires specifying two `agent*` parameters with the request which cannot
be represented with a Python `dict` as was used in the previous example.
The agent arguments can be set directly as a string in the `params`
keyword argument:
```python
resp = requests.get('https://host.of.api.com/statements/',
                    headers={'x-api-key': '12345'},
                    params=u'agent0=SMAD2&agent1=SMURF2')
```

### INDRA
Completing the circle of life, you can also access the REST API using a client
implemented in INDRA, specifically `indra.sources.indra_db_rest`. The URL and
API Key (if applicable) are configured in INDRA's config file (usually 
`~/.config/indra/config.ini`), and once added, you can get Statements by
 simply running
```python
from indra.sources import indra_db_rest as dbr
stmts = dbr.get_statements('MEK@FPLX', 'ERK@FPLX', stmt_type='Phosphorylation')
```
This API also handles more complex functionality such as implementing paging to
resolve queries that result in large amounts of content.

### Browser
You can only use a browser if there is no API key required to use the API.
If that is the case, you can simply enter the link:
`https://host.of.api.com/statements/?subject=MAP2K1&object=MAPK1`
into your browser's address bar to ge the JSON response which can
be saved.
