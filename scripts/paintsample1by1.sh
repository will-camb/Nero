#!/bin/bash

echo "Found $# arguments"
if [ "$#" -lt "4" ] ; then
    echo "Usage: paintsample1by1.sh <id> <name> <dir> <refdata.cp> <optional: removename>"
    echo "<id>: the individual which is being painted against the reference panel"
    echo "<n>: the number of the individual in the idfile"
    echo "<phaselinenumber>: the line of the phasefile containing the individual's first haplotype"
    echo "<cp_panel_scripts>: the location of the panel scripts which will be copied to each temp directory"
    exit 0
fi

name="$1"
number="$2"
phaselinenumber="$3"
cp_panel_scripts="$4"
dir="temp.$name"
chrlist=`seq 1 22`
nhaps=638

for chr in $chrlist; do
  if [ ! -f $chr.master_all_copyprobsperlocus.txt ]; then
    echo "Making $chr.master_all_copyprobsperlocus.txt"
    touch $chr.master_all_copyprobsperlocus.txt
  fi
  done

echo "Using temporary directory $dir"
mkdir -p "$dir"
mkdir -p "$dir/phasefiles"
cmd="cp $cp_panel_scripts/* $dir"
echo "Copying scripts from $cp_panel_scripts"
$cmd
echo "Done copying scripts, now making new idfile at $dir"
{ head -n $number ordered_all_pop_ids_mapped | tail -n 1 && tail -n 318 ordered_all_pop_ids_mapped; } > $dir/ordered_all_pop_ids_mapped

phaselinenumber2=$(($phaselinenumber + 1))
for chr in $chrlist ; do
  touch $dir/phasefiles/$chr.merged.phase
  head -n 3 phasefiles/$chr.merged.phase > $dir/phasefiles/$chr.merged.phase
  awk "NR>=$phaselinenumber && NR<=$phaselinenumber2" phasefiles/$chr.merged.phase >> $dir/phasefiles/$chr.merged.phase
  tail -n 636 phasefiles/$chr.merged.phase >> $dir/phasefiles/$chr.merged.phase
  sed -i '' "1s/.*/$nhaps/" $dir/phasefiles/$chr.merged.phase
  echo "Copied lines to $dir/phasefiles/$chr.merged.phase"
  done

echo "Moving to directory $dir to run remaining commands"
cd $dir
cmd="bash will_04-paintvspanel.sh"
echo $cmd
$cmd
wait
cmd="bash will_modernvsancient/painting/commandlist1.txt"
echo $cmd
$cmd
wait
bash will_04-paintvspanel.sh

for chr in $chrlist ; do
  cat $dir/will_modernvsancient/painting/$chr.all_copyprobsperlocus.txt >> ../$chr.master_all_copyprobsperlocus.txt
  rm $chr.all_copyprobsperlocus.txt
  done
rm -r $dir