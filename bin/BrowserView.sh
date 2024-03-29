#!/usr/bin/env zsh

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
     -input_folder) 
     input=$2
     shift
     ;;  
     *)

esac
shift
done

if [[ -z "$input" ]]; then
   echo ""
   echo "Usage:"
   echo "./BrowserView.sh -input_folder [input folder]"
   echo ""
   exit
fi


output=$input"_BrowserView"

mkdir -p $output
ls $input | while read one
do
  echo -n > $output/$one
  echo $one | sed 's/:/\t/g' | awk '$3=="+" {print "browser position " $1 ":" $5 "-" $2}' >> $output/$one
  echo $one | sed 's/:/\t/g' | awk '$3=="-" {print "browser position " $1 ":" $2 "-" $5}' >> $output/$one

  cat $input/$one | awk -F'\t' '{print $4"\t"$0}' | awk 'BEGIN{FS="\t";OFS="\t"}{sub(/:.*/, "", $1); print $0}' | sort -k1,1 > $output/one.sorted.tmp
  
  initCount=0 
  cat $output/one.sorted.tmp | cut -f '1' | uniq | while read OneRead
    do
      CountRead=$(( 1+$initCount ))
      ColorVar1=$(shuf -i 0-255 -n 1)
      ColorVar2=$(shuf -i 0-255 -n 1)
      ColorVar3=$(shuf -i 0-255 -n 1)
      echo "track name=read"$CountRead " description=read"$CountRead  " color="$ColorVar1","$ColorVar2","$ColorVar3 >> $output/$one
      join -t$'\t' $output/one.sorted.tmp <(echo $OneRead) | cut -f'2-' >> $output/$one
      initCount=$CountRead
    done
   
done

rm -rf $output/one.sorted.tmp
