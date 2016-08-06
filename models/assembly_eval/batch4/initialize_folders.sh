#!/usr/bin/env bash
declare -a folders=("reach" "trips" "combined")
declare -a subfolders=("index_cards" "other_outputs")

if [ ! -f pmcids.txt ]; then
    touch pmcids.txt
fi

for f in "${folders[@]}"
do
    if [ ! -d "$f" ]; then
        mkdir $f
    fi
    for sf in "${subfolders[@]}"
    do
        if [ ! -d "$f/$sf" ]; then
            mkdir $f/$sf
        fi
    done
done
