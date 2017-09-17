from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import sys
import random
import logging
try:
    from kqml import KQMLModule, KQMLPerformative, KQMLList
    have_kqml = True
except ImportError:
    KQMLModule = object
    have_kqml = False

logger = logging.getLogger('drum_reader')

class DrumReader(KQMLModule):
    def __init__(self, **kwargs):
        if not have_kqml:
            raise ImportError('Install the `pykqml` package to use ' +
                              'the DrumReader')
        self.to_read = kwargs.pop('to_read', None)
        super(DrumReader, self).__init__(name='DrumReader')
        self.msg_counter = random.randint(1, 100000)
        self.ready()
        self.extractions = []
        self.reply_counter = len(self.to_read)
        for text in self.to_read:
            self.read_text(text)

    def read_text(self, text):
        print('Reading %s' % text)
        msg_id = 'RT000%s' % self.msg_counter
        kqml_perf = _get_perf(text, msg_id)
        self.send(kqml_perf)
        self.msg_counter += 1

    def receive_reply(self, msg, content):
        extractions = content.gets(':extractions')
        self.extractions.append(extractions)
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
    dr = DrumReader(to_read=to_read)
    dr.start()
