# Visualize causal relationships amongst events

import sys
import json


def split_long_sentence(sentence, words_per_line):
    """Takes a sentence and adds a newline every "words_per_line" words.

    Parameters
    ----------
    sentence: str
        Sentene to split
    words_per_line: double
        Add a newline every this many words
    """
    words = sentence.split(' ')
    split_sentence = ''
    for i in range(len(words)):
        split_sentence = split_sentence + words[i]
        if (i+1) % words_per_line == 0:
            split_sentence = split_sentence + '\n'
        elif i != len(words) - 1:
            split_sentence = split_sentence + " "
    return split_sentence


def shorter_name(key):
    """Return a shorter name for an id.

    Does this by only taking the last part of the URI,
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


def add_event_property_edges(event_entity, entries):
    """Adds edges to the graph for event properties."""
    do_not_log = ['@type', '@id',
            'http://worldmodelers.com/DataProvenance#sourced_from']

    for prop in event_entity:
        if prop not in do_not_log:
            value = event_entity[prop]
            value_entry = None
            value_str = None
            if '@id' in value[0]:
                value = value[0]['@id']

                if value in entries:
                    value_str = get_entry_compact_text_repr(entries[value],
                                                            entries)
                #get_entry_compact_text_repr(entry, entries)

            if value_str is not None:
                edges.append([shorter_name(event_entity['@id']),
                                           shorter_name(value),
                                           shorter_name(prop)])
                node_labels[shorter_name(value)] = value_str
            #print('\t', prop, value, value_str)


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


def get_shortest_text_value(entry):
    """Given a JSON-LD entry, returns the text attribute that has the
    shortest length.

    Parameters
    ----------
    entry: dict
        A JSON-LD entry parsed into a nested python directionary via the json
        module

    Returns
    -------
    short_text: str
        Of the text values, the shortest one
    """
    text_attr = 'http://www.ontologyrepository.com/CommonCoreOntologies/has_text_value'
    if text_attr in entry:
        text_values = entry[text_attr]
        text_values = [i['@value'] for i in text_values]
        return get_shortest_string(text_values)
    else:
        return None


def get_sourced_from(entry):
    """Get a list of values from the source_from attribute"""
    sourced_from = 'http://worldmodelers.com/DataProvenance#sourced_from'

    if sourced_from in entry:
        values = entry[sourced_from]
        values = [i['@id'] for i in values]
        return values


def get_entry_compact_text_repr(entry, entries):
    """If the entry has a text value, return that.
    If the entry has a source_from value, return the text value of the source.
    Otherwise, return None."""
    text = get_shortest_text_value(entry)
    if text is not None:
        return text
    else:
        sources = get_sourced_from(entry)
        # There are a lot of references to this entity, each of which refer
        # to it by a different text label. For the sake of visualization,
        # let's pick one of these labels (in this case, the shortest one)
        if sources is not None:
            texts = []
            for source in sources:
                source_entry = entries[source]
                texts.append(get_shortest_text_value(source_entry))
            return get_shortest_string(texts)


def get_entity_type(entry):
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
    return entry_type


def load_json(visualize_file):
    # Load in the entries of the JSON-LD file and bucket by id
    entries = {}
    with open(visualize_file, 'r') as f:
        contents = f.read()
        db = json.loads(contents)

        for entry in db:
            entry_id = entry['@id']
            assert(entry_id not in entries) #IDs should be unique
            entries[entry_id] = entry
    return entries


if __name__ == '__main__':
    # First argument is the JSON file to visualize
    args = sys.argv
    if len(args) != 2:
        print('Expecting exactly one argument: the name of the JSON-LD file')
        sys.exit(1)
    visualize_file = args[1]

    entries = load_json(visualize_file)

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
        node_labels[assertion_id] = split_long_sentence(assertion_text, 5)

        # If there's an effect, add to the graph
        if cause_attr in assertion:
            cause_uri = assertion[cause_attr][0]['@id']
            cause_entry = entries[cause_uri]

            cause_type = get_entity_type(cause_entry)
            text = get_shortest_text_value(cause_entry)
            cause_uri = shorter_name(cause_uri)

            node_labels[cause_uri] = text
            edges.append([assertion_id, cause_uri, 'cause'])

            edges.append([cause_uri, cause_uri + '_type', 'type'])
            node_labels[cause_uri + '_type'] = cause_type

            add_event_property_edges(cause_entry, entries)
        else:
            cause_uri = None

        # If there's a cause, add to the graph
        assertion_id = shorter_name(assertion_id)
        if effect_attr in assertion:
            effect_uri = assertion[effect_attr][0]['@id']
            effect_entry = entries[effect_uri]

            effect_type = get_entity_type(effect_entry)
            text = get_shortest_text_value(effect_entry)
            effect_uri = shorter_name(effect_uri)

            node_labels[effect_uri] = text
            edges.append([assertion_id, effect_uri, 'effect'])

            edges.append([effect_uri, effect_uri + '_type', 'type'])
            node_labels[effect_uri + '_type'] = effect_type

            add_event_property_edges(effect_entry, entries)
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
