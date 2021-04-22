#!/bin/bash

# Script to push/pull bulk from Clement's desktop
clem_ps="sphyrna lewini"

# rsync
for species in $(ls -d * | grep -v "A.")
do
echo "Starting $species Genome"
sshpass -p "$clem_ps" rsync --exclude intervals/ -azpP clementrougeux@10.13.72.139:/Volumes/RepAdapt_4_CR/05_references/00_references/${species}/* ./${species}/
done
