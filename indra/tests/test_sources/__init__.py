# -*- coding: utf-8 -*-

"""A submodule for organizing tests for sources."""
import os

RESOURCE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'resources')


def get_resource_file(fname):
    return os.path.join(RESOURCE_PATH, fname)