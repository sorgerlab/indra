"""Handles all imports from jnius to prevent conflicts resulting from attempts
to set JVM options while the VM is already running."""

import os
import warnings
import jnius_config

if '-Xmx4g' not in jnius_config.get_options():
    if not jnius_config.vm_running:
        jnius_config.add_options('-Xmx4g')
    else:
        warnings.warn("Couldn't set memory limit for Java VM because the VM "
                      "is already running.")

path_here = os.path.dirname(os.path.realpath(__file__))
cp = path_here + '/biopax/jars/paxtools.jar'
cp_existing = os.environ.get('CLASSPATH')

if cp_existing is not None:
    os.environ['CLASSPATH'] = cp + ':' + cp_existing
else:
    os.environ['CLASSPATH'] = cp

from jnius import autoclass, JavaException, cast

