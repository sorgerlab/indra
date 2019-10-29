from os import path, mkdir
from datetime import datetime
from platform import system


class formats:
    JSON = 'json'
    TEXT = 'text'
    XML = 'xml'


def _get_dir(*args):
    dirname = path.join(*args)
    if path.isabs(dirname):
        dirpath = dirname
    elif path.exists(dirname):
        dirpath = path.abspath(dirname)
    else:
        dirpath = path.join(path.dirname(path.abspath(__file__)), dirname)
    if not path.exists(dirpath):
        mkdir(dirpath)
    return dirpath


def _time_stamp():
    return datetime.now().strftime("%Y%m%d%H%M%S")


def _get_mem_total():
    if system() == 'Linux':
        with open('/proc/meminfo', 'r') as f:
            lines = f.readlines()
        tot_entry = [line for line in lines if line.startswith('MemTotal')][0]
        ret = int(tot_entry.split(':')[1].replace('kB', '').strip())/10**6
    else:
        ret = None
    return ret
