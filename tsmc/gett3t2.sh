#!/usr/bin/env bash

[ $# -ge 1 -a -f "$1" ] && input="$1" || input="-"
cat $input | sed '/^\[/!d' | cut -d: -f 2,3,4,5 | 
    sed 's/,([0-9]:/,/g; s/[0-9]://; s/):/,/; s/);//g' | python gett3t2.py
