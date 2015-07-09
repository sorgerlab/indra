export PYTHONPATH=$PYTHONPATH:../../
export CLASSPATH=$CLASSPATH:../../belpy/biopax/biopax-pattern-04801fb13eea.jar:../../belpy/biopax/cpath-client-4.2.0-SNAPSHOT-fat-jar.jar:../../belpy/biopax/paxtools-4.3.0.jar:jena-arq-2.13.0.jar

python run_model.py --offline

dot reactions.dot -T pdf -o reactions.pdf
ccomps -x species.dot | dot | gvpack -m0 | neato -n2 -T pdf -o species.pdf

open reactions.pdf
open species.pdf

# Something with simulation figure
