# Visualize causal relationships amongst events

import sys
import json

def shorter_name(key):
    """Finds a shorter name for an id by only taking the last part of the URI,
    after the last / and the last #. Also replaces - and . with _.
    
    Parameters
    ----------
    key: str
        Some URI

    Returns
    -------
    key_short: str
        A shortened, but more ambiguous, identifier
    """
    key_short = key
    for sep in ['#', '/']:
        ind = key_short.rfind(sep)
        if ind is not None:
            key_short = key_short[ind+1:]
        else:
            key_short = key_short
    return key_short.replace('-', '_').replace('.', '_')

def get_shortest_string(string_list):
    """From a list of strings, returns the string with the shortest length.

    Parameters
    ----------
    string_list: list[str]
        A list of strings

    Returns
    -------
    shortest_string: str
        One of the strings in string_list with the shortest length
    """
    shortest_length = None
    shortest_string = None

    for s in string_list:
        if shortest_string is None or len(s) < shortest_length:
            shortest_length = len(s)
            shortest_string = s
    return shortest_string

def type_and_shortest_text_value(entry):
    """Given a JSON-LD entry, returns the abbreviated @type and the 
    text attribute that has the shortest length.

    Parameters
    ----------
    entry: dict
        A JSON-LD entry parsed into a nested python dictionary via the json
        module

    Returns
    -------
    short_type: str
        The shortest type
    short_text: str
        Of the text values, the shortest one
    """
    entry_type = entry['@type']
    entry_type = [shorter_name(t) for t in entry_type]
    entry_type = repr(entry_type)

    text_values = entry[text_attr]
    text_values = [i['@value'] for i in text_values]

    return (entry_type, get_shortest_string(text_values))

if __name__ == '__main__':
    # First argument is the JSON file to visualize
    args = sys.argv
    if len(args) != 2:
        print('Expecting exactly one argument: the name of the JSON-LD file')
        sys.exit(1)
    visualize_file = args[1]

    # Load in the entries of the JSON-LD file and bucket by id
    entries = {}
    with open(visualize_file, 'r') as f:
        contents = f.read()
        db = json.loads(contents)

        for entry in db:
            entry_id = entry['@id']
            assert(entry_id not in entries) #IDs should be unique
            entries[entry_id] = entry

    # Find all the causal assertion events
    causal_assertions = []
    for entry_id in entries:
        entry = entries[entry_id]
        if '@type' in entry:
            entry_type = entry['@type']
            entry_type = entry_type[0] # May have multiple types; check the first
            short_type = shorter_name(entry_type)
            if short_type == 'CausalAssertion':
                causal_assertions.append(entry)

    # Add causal assertions, and arguments, to a graph
    node_labels = {}
    edges = []
    cause_attr = 'http://worldmodelers.com/CauseEffect#has_cause'
    effect_attr = 'http://worldmodelers.com/CauseEffect#has_effect'
    text_attr = 'http://www.ontologyrepository.com/CommonCoreOntologies/has_text_value'
    for assertion in causal_assertions:
        assertion_id = shorter_name(assertion['@id'])

        assertion_text = assertion[text_attr][0]['@value']
        node_labels[assertion_id] = assertion_text


        # If there's an effect, add to the graph
        if cause_attr in assertion:
            cause_uri = assertion[cause_attr][0]['@id']
            cause_entry = entries[cause_uri]

            (cause_type, text) = type_and_shortest_text_value(cause_entry)
            cause_uri = shorter_name(cause_uri)

            node_labels[cause_uri] = text
            edges.append([assertion_id, cause_uri, 'cause'])

            edges.append([cause_uri, cause_uri + '_type', 'type'])
            node_labels[cause_uri + '_type'] = cause_type
        else:
            cause_uri = None

        # If there's a cause, add to the graph
        assertion_id = shorter_name(assertion_id)
        if effect_attr in assertion:
            effect_uri = assertion[effect_attr][0]['@id']
            effect_entry = entries[effect_uri]

            (effect_type, text) = type_and_shortest_text_value(effect_entry)
            effect_uri = shorter_name(effect_uri)

            node_labels[effect_uri] = text
            edges.append([assertion_id, effect_uri, 'effect'])

            edges.append([effect_uri, effect_uri + '_type', 'type'])
            node_labels[effect_uri + '_type'] = effect_type
        else:
            effect_uri = None

    # Write out the graph in graphviz representation
    with open('graph.dot', 'w') as f:
        f.write('DiGraph {\n')
        # Write edges
        for edge in edges:
            f.write('\t%s -> %s [label="%s"]\n' % (
                edge[0], edge[1], edge[2]))

        # Write node labels
        f.write('//Labels\n')
        print('Num labels:', len(node_labels))
        for node in node_labels:
            f.write('\t%s [label="%s"]\n' % (node, node_labels[node]))
        f.write('}\n')

    print('Causal assertions found:', len(causal_assertions))
    #has_cause
    #has_effect
