#!/bin/bash
cat f | awk '{print $1 "\t" $2 "\t" $3 "\t" $1":"$2"-"$3 "\t" "." "\t" "+"}' | tail -n +2
