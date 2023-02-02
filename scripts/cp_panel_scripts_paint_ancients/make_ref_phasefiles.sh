#!/bin/bash

echo "Found $# arguments"
if [ "$#" -lt "3" ] ; then
    echo "Usage:  make_ref_phasefiles.sh <id> <phaselinenumber> <phaselinenumber2>"
    echo "<id>: the individual which is being painted against the reference panel"
    echo "<phaselinenumber>: the  line of the phasefile containing the individual's first haplotype"
    echo "<phaselinenumber2>: the line of the phasefile containing the individual's second haplotype"
    exit 0
fi
set -e

name="$1"
phaselinenumber="$2"
phaselinenumber2="$3"
dir="temp.$name"
chrlist=`seq 1 22`
nhaps=636

mkdir -p ref_phasefiles

for chr in $chrlist ; do
  if [ ! -f  ref_phasefiles/"$chr".merged.phase ]; then
    echo "Making ref phasefile for chromosome $chr at ref_phasefiles/"
    touch ref_phasefiles/"$chr".merged.phase
    head -n 3 phasefiles/filtered."$chr".phase > ref_phasefiles/"$chr".merged.phase
    sed -i "1s/.*/$nhaps/" ref_phasefiles/"$chr".merged.phase
  fi
  done

for chr in $chrlist ; do
  awk "NR>=$phaselinenumber && NR<=$phaselinenumber2" phasefiles/filtered."$chr".phase >> ref_phasefiles/"$chr".merged.phase
  done

