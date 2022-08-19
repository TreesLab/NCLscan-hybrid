#!/usr/bin/env zsh

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
     -input_folder) 
     input=$2
     shift
     ;;  
     -o | --output)
     out=$2
     shift
     ;;
     -L)
     extend_len=$2
     shift
     ;;
     *)

esac
shift
done

if [[ -z "$input" ]]; then
   echo ""
   echo "Usage:"
   echo "./OutOfCircle.sh -input_folder [input folder] -o [output_prefix] -L [extend_len]"
   echo ""
   exit
fi


if [[ -z "$out" ]]; then
   echo ""
   echo "Usage:"
   echo "./OutOfCircle.sh -input_folder [input folder] -o [output_prefix] -L [extend_len]"
   echo ""
   exit
fi

if [[ -z "$extend_len" ]]; then
   extend_len=100
fi


# Out of circle: long read extends more than 100 bp on one side (upstream or downstream)  
mkdir -p OC.events
echo -n > $out\_OC.result  
echo -n > $out\_one.OCreadID.tmp

ls $input | while read one
do
   one_name=$(echo $one | sed 's/\.bed12//g')
   intra_chr=$(echo $one | sed 's/:/\t/g' | awk '{print $1}')
   intra_donor=$(echo $one | sed 's/:/\t/g' | awk '{print $2}')
   intra_acceptor=$(echo $one | sed 's/:/\t/g' | awk '{print $5}')
   upstream_pos=$(echo $intra_donor $intra_acceptor | tr ' ' \\n | sort | head -n 1 | awk '{print $1-L}' L=$extend_len)
   downstream_pos=$(echo $intra_donor $intra_acceptor | tr ' ' \\n | sort | tail -n 1 | awk '{print $1+L}' L=$extend_len)
   long_list=$(cat $input/$one | awk '{print $4}' | sed 's/:/\t/g' | awk '{print $1}' | sort  | uniq)
   num_long=$(echo $long_list | wc -l)

   cat $input/$one | awk -F'\t' '{print $4"\t"$0}' | awk 'BEGIN{FS="\t";OFS="\t"}{sub(/:.*/, "", $1); print $0}' | sort -k1,1 > $out\_one.sorted.tmp
   
   OC=0
   countOC=0
   echo $long_list | while read oneLong 
   do
      pieces_OneLongRead=$(join -t$'\t' $out\_one.sorted.tmp <(echo $oneLong) | wc -l)
      upstream_OneLongRead=$(join -t$'\t' $out\_one.sorted.tmp <(echo $oneLong) | cut -f '2-' | awk '$1==chr' chr=$intra_chr | awk '{print $2 " " $3}' | tr ' ' \\n | sort | head -n 1)
      downstream_OneLongRead=$(join -t$'\t' $out\_one.sorted.tmp <(echo $oneLong) | cut -f '2-' | awk '$1==chr' chr=$intra_chr | awk '{print $2 " " $3}' | tr ' ' \\n | sort | tail -n 1)
      
      if [[ "$pieces_OneLongRead" -eq 2 ]]; then
         if [[ "$upstream_OneLongRead" -lt "$upstream_pos" ]] || [[ "$downstream_OneLongRead" -gt "$downstream_pos" ]]
           then {
            OC=1  countOC=$(( $countOC + 1 ))
            echo $oneLong >>  $out\_one.OCreadID.tmp
           }
          else OC=0  countOC=$(( $countOC ))
         fi
      fi  
    done

   join $out\_one.sorted.tmp <(cat $out\_one.OCreadID.tmp | sort -k1,1) | tr ' ' \\t | cut -d$'\t' -f 2- > $out/OC_events/$one

      

   if [[ "$countOC" -gt 0  ]] 
        then OC=1
        else OC=0 
      fi
      

   echo $one_name   $OC   $countOC >> $out\_OC.result  
done

rm -r -f $out\_one.sorted.tmp
rm -r -f $out\_one.OCreadID.tmp
