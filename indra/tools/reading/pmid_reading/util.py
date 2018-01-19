from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from platform import system


def get_proc_num():
    if system() == 'Linux':
        with open('/proc/cpuinfo', 'r') as f:
            ret = len([
                line for line in f.readlines() if line.startswith('processor')
                ])
    else:
        ret = None
    return ret


def get_mem_total():
    if system() == 'Linux':
        with open('/proc/meminfo', 'r') as f:
            lines = f.readlines()
        tot_entry = [line for line in lines if line.startswith('MemTotal')][0]
        ret = int(tot_entry.split(':')[1].replace('kB', '').strip())/10**6
    else:
        ret = None
    return ret


