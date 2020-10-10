#!/bin/bash
nfiles=$(ls | grep -c root) 
echo "$nfiles jobs completed"
if [[ $nfiles -ge $1 ]]; then
    exit 0
else
    exit 1
fi
