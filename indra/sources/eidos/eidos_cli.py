"""
This is a Python based command line interface to Eidos
to complement the Python-Java bridge based interface.
EIDOSPATH (in the INDRA config.ini or as an environmental variable)
needs to be pointing to a fat JAR of the Eidos system.
"""
import os
import glob
import logging
import subprocess
from indra import get_config
from .eidos_api import process_json_ld_file


eip = get_config('EIDOSPATH')
eidos_package = 'org.clulab.wm.eidos'
logger = logging.getLogger('eidos_cli')


def run_eidos(endpoint, *args):
    """Run a given enpoint of Eidos through the command line.

    Parameters
    ----------
    endpoint : str
        The class within the Eidos package to run, for instance
        'apps.ExtractFromDirectory' will run
        'org.clulab.wm.eidos.apps.ExtractFromDirectory'
    *args
        Any further arguments to be passed as inputs to the class
        being run.
    """
    # Make the full path to the class that should be used
    call_class = '%s.%s' % (eidos_package, endpoint)
    # Assemble the command line command and append optonal args
    cmd = ['java', '-Xmx12G', '-cp', eip, call_class] + list(args)
    logger.info('Running Eidos with command "%s"' % (' '.join(cmd)))
    subprocess.call(cmd)


def extract_from_directory(path_in, path_out):
    """Run Eidos on a set of text files in a folder.

    The output is produced in the specified output folder but
    the output files aren't processed by this function.

    Parameters
    ----------
    path_in : str
        Path to an input folder with some text files
    path_out : str
        Path to an output folder in which Eidos places the output
        JSON-LD files
    """
    path_in = os.path.realpath(os.path.expanduser(path_in))
    path_out = os.path.realpath(os.path.expanduser(path_out))
    logger.info('Running Eidos on input folder %s' % path_in)
    run_eidos('apps.ExtractFromDirectory', path_in, path_out)


def extract_and_process(path_in, path_out):
    """Run Eidos on a set of text files and process output with INDRA.

    The output is produced in the specified output folder but
    the output files aren't processed by this function.

    Parameters
    ----------
    path_in : str
        Path to an input folder with some text files
    path_out : str
        Path to an output folder in which Eidos places the output
        JSON-LD files

    Returns
    -------
    stmts : list[indra.statements.Statements]
        A list of INDRA Statements
    """
    path_in = os.path.realpath(os.path.expanduser(path_in))
    path_out = os.path.realpath(os.path.expanduser(path_out))
    extract_from_directory(path_in, path_out)
    jsons = glob.glob(os.path.join(path_out, '*.jsonld'))
    logger.info('Found %d JSON-LD files to process in %s' %
                (len(jsons), path_out))
    stmts = []
    for json in jsons:
        ep = process_json_ld_file(json)
        if ep:
            stmts += ep.statements
    return stmts
