#!/bin/sh

if [ "$#" -ne "1" ] ; then
    echo "Usage: will_05-paintvspanel1by1.sh <idfile>"
    echo "<idfile>: A file with the IDs of all individuals to be painted against the reference panel, 1 per row"
    exit 0
fi

idfile="$1"
