from indra.sources.medscan.api import MedscanRelation

class MedscanProcessor(object):
    def __init__(self):
        self.statements = []

    def process_relation(self, relation):
        if relation.verb == 'UnknownRegulation-positive':
            pass
        elif relation.verb == 'UnknownRegulation-negative':
            pass
        elif relation.verb == 'UnknownRegulation-unknown':
            pass
        elif relation.verb == 'ExpressionControl-positive':
        pass


