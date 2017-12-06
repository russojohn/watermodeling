#!/bin/bash

if [ $# -ne 2 ]; then echo "$0 [xyz file] [property file]" ; exit ; fi

head -n2 $1 ; paste <(sed '1,2d' $1) $2
