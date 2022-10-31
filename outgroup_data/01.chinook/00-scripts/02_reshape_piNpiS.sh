#!/bin/bash

#input consist of the file geneID for each pop.

input=$1

grep -v "fst.clean.fst.clean.fst.sites4foldonly"  $input|\
    sed 's/.fst.clean.fst.clean.fst//g' > "$input"_cleaned
