class Publication(object):
    def __init__(self, interaction):
        self.pmid = "PMID" + str(interaction['PUBMED_ID'])
        self.modification = interaction['MODIFICATION']
        self.experimental_system = interaction[
            'EXPERIMENTAL_SYSTEM']
        self.experimental_system_type = interaction[
            'EXPERIMENTAL_SYSTEM_TYPE']
        self.throughput = interaction['THROUGHPUT']


    def __str__(self):
        return ("PMID:" + self.pmid)


    def __repr__(self):
        return ("<pmid: %s at 0x%x>" %
                (self.pmid, id(self)))

        
