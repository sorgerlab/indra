[![Build Status](https://travis-ci.org/sorgerlab/indra.svg?branch=travis_ci)](https://travis-ci.org/sorgerlab/indra)

INDRA
=====
INDRA (Integrated Dynamical Reasoning Assembler) generates executable models of pathway dynamics from natural language (using the TRIPS parser), BioPAX and BEL sources (including the PathwayCommons database and NDEx).

Usage example
-------------
In this example a model is generated from the natural language description of a mechanism. 

```python
from indra.pysb_assembler import PysbAssembler
from indra.trips import trips_api
pa = PysbAssembler()
# Process a natural language description of a mechanism
tp = trips_api.process_text('MEK2 phosphorylates ERK1 at Thr-202 and Tyr-204')
# Collect extracted mechanisms in PysbAssembler
pa.add_statements(tp.statements)
# Assemble the model
model = pa.make_model(policies='two_step')
```
