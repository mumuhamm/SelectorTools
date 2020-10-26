#!/bin/bash
if [[ -f $1 ]]; then
    rm $(find . -regex '.*_.*_[0-9]+\.root')
fi
