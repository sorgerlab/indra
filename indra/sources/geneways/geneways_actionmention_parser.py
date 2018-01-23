from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import codecs

import inspect

class GenewaysActionMention(object):
    """Represents a row of data in the Geneways human_action.txt, structured
    so you can access by field."""
    def __init__(self, textRow):
        """Parses a row of text data in human_action.txt and sets the
        field name of this object to the corresponding data."""
        tokens = textRow.split('\t')
        if len(tokens) != 11:
            m = 'Expected 11 tokens for each line of human_action.txt but' + \
                    'got %d tokens: "%s"' % (len(tokens), textRow)
            raise Exception(m)

        self.hiid = tokens[0]
        self.actionmentionid = tokens[1]
        self.negative = tokens[2]
        self.upstream = tokens[3]
        self.actiontype = tokens[4]
        self.downstream = tokens[5]
        self.pmid = tokens[6]
        self.isFullText = tokens[7]
        self.sentencenumber = tokens[8]
        self.score = tokens[9]
        self.prec = tokens[10]

    def make_annotation(self):
        """Returns a dictionary with all properties of the action mention."""
        annotation = dict()

        # Put all properties of the action object into the annotation
        for item in dir(self):
            if len(item) > 0 and item[0] != '_' and \
                    not inspect.ismethod(getattr(self, item)):
                annotation[item] = getattr(self, item)

        return annotation

    def __repr__(self):
        r = ''
        first = True
        for item in dir(self):
            if len(item) > 0 and item[0] != '_' and \
                    not inspect.ismethod(getattr(self, item)):
                if not first:
                    r = r + ","

                r = r + item + "=" + getattr(self, item)
                first = False
        return r

class GenewaysActionMentionParser(object):
    """Parses a human_actionmention.txt file, and populates a list of
    GenewaysActionMention objects with these data."""

    def __init__(self, actionmention_filename):
        """Parses the file and populates the data."""

        self.action_mentions = list()

        f = codecs.open(actionmention_filename, 'r', encoding='latin-1')
        first_line = True
        line_number = 1
        for line in f:
            line = line.rstrip()
            if first_line:
                # Ignore the first line
                first_line = False
            else:
                self.action_mentions.append(GenewaysActionMention(line))
            line_number = line_number + 1
