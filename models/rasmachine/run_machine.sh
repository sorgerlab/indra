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
export PYTHONPATH=../../../indra:..
twitter_cred="$model_name/twitter_cred.txt"
gmail_cred="$model_name/gmail_cred.txt"
ndex_cred="$model_name/ndex_cred.txt"
log_file="$model_name/log.txt"
err_file="$model_name/err.txt"
belief_file="$model_name/belief.txt"
python rasmachine.py --model $model_name --twitter $twitter_cred --gmail $gmail_cred --ndex $ndex_cred --belief $belief_file 1>> $log_file 2>> $err_file
