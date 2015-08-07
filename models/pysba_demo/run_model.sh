export PYTHONPATH=$PYTHONPATH:../../
export CLASSPATH=$CLASSPATH:../../indra/biopax/biopax-pattern-04801fb13eea.jar:../../indra/biopax/cpath-client-4.2.0-SNAPSHOT-fat-jar.jar:../../indra/biopax/paxtools-4.3.0.jar:jena-arq-2.13.0.jar

python run_model.py --render --simulate

dot reactions.dot -T pdf -o reactions.pdf
ccomps -x species.dot | dot | gvpack -m0 | neato -n2 -T pdf -o species.pdf

open reactions.pdf
open species.pdf

