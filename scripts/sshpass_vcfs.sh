#!/bin/#!/usr/bin/env bash

# Loop to move over all VCF from Clements desktop to yeaman03
clem_ps="sphyrna lewini"
MAIN=/lu213/james.whiting/RepAdapt

# sshpass from volume 1
for species in $(ls -d *)
do
echo "Starting $species VCF"
sshpass -p "$clem_ps" rsync -azpP clementrougeux@10.13.72.139:/Volumes/RepAdapt_3_CR/04_genotypes/$species/03*/*concat* $species/
done

# sshpass for reference genomes

for species in $(ls -d * | grep -v "A.")
do
echo "Starting $species Genome"
sshpass -p "$clem_ps" rsync --exclude intervals/ -azpP clementrougeux@10.13.72.139:/Volumes/RepAdapt_4_CR/05_references/00_references/${species}/* ./${species}/
done
