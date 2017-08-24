import sys
import random
from indra.sources import trips
from kqml import KQMLModule, KQMLPerformative, KQMLList

class DrumReader(KQMLModule):
    def __init__(self, **kwargs):
        super(DrumReader, self).__init__(**kwargs)
        self.msg_counter = random.randint(1, 100000)
        self.ready()
        self.extractions = None
        self.read_text('MEK phosphorylates ERK1.')
        self.read_text('BRAF phosphorylates MEK1.')

    def read_text(self, text):
        msg_id = 'RT000%s' % self.msg_counter
        kqml_perf = _get_perf(text, msg_id)
        self.send(kqml_perf)
        self.msg_counter += 1

    def receive_reply(self, msg, content):
        extractions = content.gets(':extractions')
        self.extractions = extractions
        tp = trips.process_xml(self.extractions)
        print(tp.statements)

def _get_perf(text, msg_id):
    text = text.encode('utf-8')
    msg = KQMLPerformative('REQUEST')
    msg.set('receiver', 'DRUM')
    content = KQMLList('run-text')
    content.sets('text', text)
    msg.set('content', content)
    msg.set('reply-with', msg_id)
    return msg

if __name__ == '__main__':
    # NOTE: drum/bin/trips-drum needs to be running
    dr = DrumReader(name='DrumReader')
    dr.start()
