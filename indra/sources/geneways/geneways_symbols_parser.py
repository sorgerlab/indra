from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import codecs

class GenewaysSymbols(object):
    """Parser for the Geneways symbol table. Once parsed, provides a mapping
    from gene symbol names to Entrez IDs"""
    def __init__(self, symbols_filename):
        """Load in the symbol table, row by row, and populate hashmaps
        linking the Geneways ids to Entrez ids and vice versa."""

        self.symbols_to_ids = dict()
        self.ids_to_symbols = dict()
        self.appears_multiple_times = list()

        f = codecs.open(symbols_filename, 'r', encoding='latin-1')
        first_line = True
        for line in f:
            line = line.rstrip()
            if first_line:
                # Don't parse the first line
                first_line = False
            else:
                tokens = line.split()

                if len(tokens) != 2:
                    m = "Apologies! For Geneways symbol table %s, expected' + \
                            ' two tokens per line"
                    m = m % symbols_filename
                    raise Exception(m)

                entrez_id = tokens[0]
                symbol = tokens[1]

                if symbol in self.symbols_to_ids:
                    self.appears_multiple_times.append(symbol)
                else:
                    self.symbols_to_ids[symbol] = list()

                if entrez_id in self.ids_to_symbols:
                    m = 'Apologies! Entrez ID listed multiple times in %s'
                    m = m % symbols_filename
                    raise Exception(m)

                self.symbols_to_ids[symbol].append(entrez_id)
                self.ids_to_symbols[entrez_id] = symbol

        f.close()

    def symbol_to_id(self, symbol):
        """Returns the list of Entrez IDs for a given Geneways symbol
        (there may be more than one)"""

        if symbol not in self.symbols_to_ids:
            m = 'Could not look up Entrez ID for Geneways symbol ' + symbol
            raise Exception(m)
        return self.symbols_to_ids[symbol]

    def id_to_symbol(self, entrez_id):
        """Gives the symbol for a given entrez id)"""

        entrez_id = str(entrez_id)
        if entrez_id not in self.ids_to_symbols:
            m = 'Could not look up symbol for Entrez ID ' + entrez_id
            raise Exception(m)
        return self.ids_to_symbols[entrez_id]

    def symbols_with_multiple_ids(self):
        return self.appears_multiple_times
