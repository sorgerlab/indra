"""Handles all imports from jnius to prevent conflicts resulting from attempts
to set JVM options while the VM is already running."""

import jnius_config

if '-Xmx4g' not in jnius_config.get_options():
    if not jnius_config.vm_running:
        jnius_config.add_options('-Xmx4g')
    else:
        warnings.warn("Couldn't set memory limit for Java VM because the VM "
                      "is already running.")

from jnius import autoclass, JavaException, cast

