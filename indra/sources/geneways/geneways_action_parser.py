from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import numpy as np
from os import path
import codecs
import inspect
from indra.sources.geneways.geneways_actionmention_parser \
        import GenewaysActionMentionParser
from indra.sources.geneways.geneways_symbols_parser import GenewaysSymbols


class GenewaysAction(object):
    """Represents a row of data in the Geneways human_action.txt,
    structured so you can access by field."""
    def __init__(self, text_row):
        """Parses a row of text data in human_action.txt and sets the
        field name of this object to the corresponding data."""
        tokens = text_row.split('\t')
        if len(tokens) != 9:
            msg = 'Expected 9 tokens for each line of human_action.txt' + \
                ' but got %d tokens: "%s"' % (len(tokens), text_row)
            raise Exception(msg)

        self.hiid = tokens[0]
        self.up = tokens[1]
        self.dn = tokens[2]
        self.actiontype = tokens[3]
        self.action_count = tokens[4]
        self.actionmention_count = tokens[5]
        self.plo = tokens[6]
        self.max_score = tokens[7]
        self.max_prec = tokens[8]
        self.action_mentions = list() # Initially empty, can be populated later

    def make_annotation(self):
        """Returns a dictionary with all properties of the action
        and each of its action mentions."""
        annotation = dict()

        # Put all properties of the action object into the annotation
        for item in dir(self):
            if len(item) > 0 and item[0] != '_' and \
                    not inspect.ismethod(getattr(self, item)):
                annotation[item] = getattr(self, item)

        # Add properties of each action mention
        annotation['action_mentions'] = list()
        for action_mention in self.action_mentions:
            annotation_mention = action_mention.make_annotation()
            annotation['action_mentions'].append(annotation_mention)

        return annotation

    def __repr__(self):
        r = ''
        first = True
        for item in dir(self):
            if len(item) > 0 and item[0] != '_' and \
                    not inspect.ismethod(getattr(self, item)):

                if not first:
                    r = r + ","

                r = r + item + "=" + repr(getattr(self, item))
                first = False
        return r


class GenewaysActionParser(object):
    """Parses a human_action.txt file, and populates
    a list of GenewaysAction objects with these data."""

    def __init__(self, input_folder):
        """Parses the file and populations the action data"""

        f = 'human_action.txt'
        action_filename = self._search_path(input_folder, f)

        f = 'human_actionmention.txt'
        actionmention_filename = self._search_path(input_folder, f)

        f = 'human_symbols.txt'
        symbols_filename = self._search_path(input_folder, f)

        if action_filename is None or actionmention_filename is None \
            or symbols_filename is None:
            msg = 'Could not find Geneways extracted data: ' + \
                '(human_action.txt, human_actionmention.txt, ' + \
                'human_symbols.txt) in %s' % input_folder
            raise Exception(msg)

        self._init_action_list(action_filename)
        self._link_to_action_mentions(actionmention_filename)
        self._lookup_symbols(symbols_filename)

    def _search_path(self, directory_name, filename):
        """Searches for a given file in the specified directory."""
        full_path = path.join(directory_name, filename)
        if path.exists(full_path):
            return full_path

        # Could not find the requested file in any of the directories
        return None

    def _init_action_list(self, action_filename):
        """Parses the file and populates the data."""

        self.actions = list()
        self.hiid_to_action_index = dict()

        f = codecs.open(action_filename, 'r', encoding='latin-1')
        first_line = True
        for line in f:
            line = line.rstrip()
            if first_line:
                # Ignore the first line
                first_line = False
            else:
                self.actions.append(GenewaysAction(line))

                latestInd = len(self.actions)-1
                hiid = self.actions[latestInd].hiid
                if hiid in self.hiid_to_action_index:
                    raise Exception('action hiid not unique: %d' % hiid)
                self.hiid_to_action_index[hiid] = latestInd

    def _link_to_action_mentions(self, actionmention_filename):
        """Add action mentions"""
        parser = GenewaysActionMentionParser(actionmention_filename)
        self.action_mentions = parser.action_mentions

        for action_mention in self.action_mentions:
            hiid = action_mention.hiid
            if hiid not in self.hiid_to_action_index:
                m1 = 'Parsed action mention has hiid %d, which does not exist'
                m2 = ' in table of action hiids'
                raise Exception((m1 + m2) % hiid)
            else:
                idx = self.hiid_to_action_index[hiid]
                self.actions[idx].action_mentions.append(action_mention)

    def _lookup_symbols(self, symbols_filename):
        """Look up symbols for actions and action mentions"""
        symbol_lookup = GenewaysSymbols(symbols_filename)
        for action in self.actions:
            action.up_symbol = symbol_lookup.id_to_symbol(action.up)
            action.dn_symbol = symbol_lookup.id_to_symbol(action.dn)

    def get_top_n_action_types(self, top_n):
        """Returns the top N actions by count."""
        # Count action types
        action_type_to_counts = dict()
        for action in self.actions:
            actiontype = action.actiontype
            if actiontype not in action_type_to_counts:
                action_type_to_counts[actiontype] = 1
            else:
                action_type_to_counts[actiontype] = \
                        action_type_to_counts[actiontype] + 1

        # Convert the dictionary representation into a pair of lists
        action_types = list()
        counts = list()
        for actiontype in action_type_to_counts.keys():
            action_types.append(actiontype)
            counts.append(action_type_to_counts[actiontype])

        # How many actions in total?
        num_actions = len(self.actions)
        num_actions2 = 0
        for count in counts:
            num_actions2 = num_actions2 + count
        if num_actions != num_actions2:
            raise(Exception('Problem counting everything up!'))

        # Sort action types by count (lowest to highest)
        sorted_inds = np.argsort(counts)
        last_ind = len(sorted_inds)-1

        # Return the top N actions
        top_actions = list()
        if top_n > len(sorted_inds):
            raise Exception('Asked for top %d action types, ' +
                            'but there are only %d action types'
                            % (top_n, len(sorted_inds)))
        for i in range(top_n):
            top_actions.append(action_types[sorted_inds[last_ind-i]])
        return top_actions
