from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
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
    """Agent which processes text through a local TRIPS/DRUM instance.

    This class is implemented as a communicative agent which sends and receives
    KQML messages through a socket. It sends text (ideally in small blocks
    like one sentence at a time) to the running DRUM instance and receives
    extraction knowledge base (EKB) XML responses asynchronously through
    the socket. To install DRUM and its dependencies locally, follow
    instructions at: https://github.com/wdebeaum/drum
    Once installed, run `drum/bin/trips-drum -nouser` to run DRUM without
    a GUI. Once DRUM is running, this class can be instantiated as
    `dr = DrumReader()`, at which point it attempts to
    connect to DRUM via the socket. You can use `dr.read_text(text)` to
    send text for reading.
    In another usage more, `dr.read_pmc(pmcid)` can be used to read
    a full open-access PMC paper.
    Receiving responses can be started as `dr.start()` which waits for
    responses from the reader and returns when all responses were received.
    Once finished, the list of EKB XML extractions can be accessed via
    `dr.extractions`.

    Attributes
    ----------
    extractions : list[str]
        A list of EKB XML extractions corresponding to the input text list.
    """
    def __init__(self, **kwargs):
        if not have_kqml:
            raise ImportError('Install the `pykqml` package to use ' +
                              'the DrumReader')
        super(DrumReader, self).__init__(name='DrumReader', **kwargs)
        self.msg_counter = random.randint(1, 100000)
        self.ready()
        self.extractions = []
        self.reply_counter = 0

    def read_pmc(self, pmcid):
        """Read a given PMC article.

        Parameters
        ----------
        pmcid : str
            The PMC ID of the article to read. Note that only
            articles in the open-access subset of PMC will work.
        """
        msg = KQMLPerformative('REQUEST')
        msg.set('receiver', 'DRUM')
        content = KQMLList('run-pmcid')
        content.sets('pmcid', pmcid)
        content.set('reply-when-done', 'true')
        msg.set('content', content)
        msg.set('reply-with', 'P-%s' % pmcid)
        self.reply_counter += 1
        self.send(msg)

    def read_text(self, text):
        """Read a given text phrase.

        Parameters
        ----------
        text : str
            The text to read. Typically a sentence or a paragraph.
        """
        logger.info('Reading: "%s"' % text)
        msg_id = 'RT000%s' % self.msg_counter
        kqml_perf = _get_perf(text, msg_id)
        self.reply_counter += 1
        self.msg_counter += 1
        self.send(kqml_perf)

    def receive_reply(self, msg, content):
        """Handle replies with reading results."""
        reply_head = content.head()
        if reply_head == 'error':
            comment = content.gets('comment')
            logger.error('Got error reply: "%s"' % comment)
        else:
            extractions = content.gets('extractions')
            self.extractions.append(extractions)
        self.reply_counter -= 1
        if self.reply_counter == 0:
            self.exit(0)


def _get_perf(text, msg_id):
    """Return a request message for a given text."""
    msg = KQMLPerformative('REQUEST')
    msg.set('receiver', 'DRUM')
    content = KQMLList('run-text')
    content.sets('text', text)
    msg.set('content', content)
    msg.set('reply-with', msg_id)
    return msg
