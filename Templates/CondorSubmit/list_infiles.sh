#!/bin/bash
echo "transfer_input_files =" $(ls *.root | grep -v / | tr '\n' ',')
