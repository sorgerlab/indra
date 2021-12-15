from indra.statements import *


class AcsnProcessor:
    def __init__(self, relations_df, correspondence_df):
        self.relations_df = relations_df
        self.correspondence_df = correspondence_df
        self.statements = []

    def extract_statements(self):
        # This is where we implement the Statement extraction
        pass