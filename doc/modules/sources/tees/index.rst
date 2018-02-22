TEES (:py:mod:`indra.sources.tees`)
===================================

The TEES processor requires an installaton of TEES. To install TEES:

1. Clone the latest stable version of TEES using

    `git clone https://github.com/jbjorne/TEES.git`

2. Put this TEES cloned repository in one of these three places:
   the same directory as INDRA, your home directory, or ~/Downloads.
   If you put TEES in a location other than one of these three
   places, you will need to pass this directory to
   `indra.sources.tees.tees_api.process_text` each time you call it.

3. Run configure.py within the TEES installation to install TEES dependencies.

TEES API (:py:mod:`indra.sources.tees.tees_api`)
------------------------------------------------

.. automodule:: indra.sources.tees.tees_api
    :members:

TEES Processor (:py:mod:`indra.sources.tees.processor`)
-------------------------------------------------------

.. automodule:: indra.sources.tees.processor
    :members:

