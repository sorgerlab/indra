import sys
import random
from indra.sources import trips
from kqml import KQMLModule, KQMLPerformative, KQMLList

class DrumReader(KQMLModule):
    def __init__(self, **kwargs):
        self.to_read = kwargs.pop('to_read', None)
        super(DrumReader, self).__init__(**kwargs)
        self.msg_counter = random.randint(1, 100000)
        self.ready()
        self.extractions = []
        for text in to_read:
            self.read_text(text)
        self.reply_counter = len(to_read)

    def read_text(self, text):
        print('Reading %s' % text)
        msg_id = 'RT000%s' % self.msg_counter
        kqml_perf = _get_perf(text, msg_id)
        self.send(kqml_perf)
        self.msg_counter += 1

    def receive_reply(self, msg, content):
        extractions = content.gets(':extractions')
        tp = trips.process_xml(extractions)
        self.extractions += tp.statements
        self.reply_counter -= 1
        if self.reply_counter == 0:
            self.exit(0)

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
    to_read = ['MEK phosphorylates ERK1.', 'BRAF phosphorylates MEK1.']
    dr = DrumReader(name='DrumReader', to_read=to_read)
    dr.start()
