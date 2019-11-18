import logging

from .core import Reader, ReadingError, ReadingData, get_reader, \
    get_reader_classes, get_reader_class, EmptyReader

from .util import get_dir

from .content import Content

logger = logging.getLogger(__name__)

err_msg = ("Could not load {reader} reader: \"{err}\". {reader} will not be "
           "available.")

# Try loading each reader, and raise a warning if it not available.
try:
    from .isi import IsiReader
except Exception as e:
    logger.warning(err_msg.format(reader="ISI", err=str(e)))

try:
    from .trips import TripsReader
except Exception as e:
    logger.warning(err_msg.format(reader="Trips", err=str(e)))

try:
    from .reach import ReachReader
except Exception as e:
    logger.warning(err_msg.format(reader="REACH", err=str(e)))

try:
    from .sparser import SparserReader
except Exception as e:
    logger.warning(err_msg.format(reader="Sparser", err=str(e)))
