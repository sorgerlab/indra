#!/bin/bash
model_name=$1
export CLASSPATH="$HOME/src/reach/*"
platform=`uname`
if [[ "$platform" == 'Linux' ]]; then
    export JAVA_HOME="/usr/lib/jvm/java-8-oracle/"
else
    export JAVA_HOME="/Library/Java/JavaVirtualMachines/jdk1.8.0_74.jdk/Contents/Home"
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"
export PYTHONPATH=$PYTHONPATH:../../../indra
echo $PYTHONPATH
log_file="$model_name/log.txt"
config_file="$model_name/config.yaml"
python rasmachine.py $model_name $config_file >> $log_file 2>> $log_file
