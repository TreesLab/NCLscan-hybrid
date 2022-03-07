#!/usr/bin/env zsh

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
     -input_folder) 
     input=$2
     shift
     ;;  
     -c | --config)
     config=$2
     shift
     ;;
     -o | --output)
     out=$2
     shift
     ;;
     *)

esac
shift
done

if [[ -z "$input" ]]; then
   echo ""
   echo "Usage:"
   echo "./WithinCircle.sh -input_folder [input folder] -c [configure file link] -o [output_prefix]"
   echo ""
   exit
fi


if [[ -z "$out" ]]; then
   echo ""
   echo "Usage:"
   echo "./WithinCircle.sh -input_folder [input folder] -c [configure file link] -o [output_prefix]"
   echo ""
   exit
fi

source $config

# Within Circle : obtain more than one pseudo-reference (200 bp NCL sequence) 
echo $input
echo $out
echo -n > $out/$out\_WithinCircle.result  

mkdir -p \
    $out/WithinCircle_tmp \
    $out/WithinCircle_events \
    $out/WithinCircle_tmp/WithinCircle_reads_tmp


cat $input | awk '{print $1 "\t" $5 "\t" $6}' | sort -k3,3 -k1,1 | uniq -c | awk '$1 > 1' | sed 's/:/\t/g' | awk '$4==$7' | awk '{print $2 "\t" $3 "\t" $4":"$5":"$6":"$7":"$8":"$9}' > $out/WithinCircle_tmp/FlankingRead_80_twice.list
cat $out/WithinCircle_tmp/FlankingRead_80_twice.list | cut -f '1' | sort | uniq > $out/WithinCircle_tmp/FlankingRead_80_twice.readID

cat $out/WithinCircle_tmp/FlankingRead_80_twice.list > $out/WithinCircle_tmp/FlankingRead_80_twice_intra.list
cat $out/WithinCircle_tmp/FlankingRead_80_twice_intra.list | awk '{print $1}' | sort | uniq > $out/WithinCircle_tmp/FlankingRead_80_twice_intra.readID

cat $input | sort -k1,1 > $out/WithinCircle_tmp/FlankingRead_80.sorted_by_read_id.paf
join -t$'\t' $out/WithinCircle_tmp/FlankingRead_80.sorted_by_read_id.paf $out/WithinCircle_tmp/FlankingRead_80_twice_intra.readID > $out/WithinCircle_tmp/FlankingRead_80_twice_intra.paf

cat $out/WithinCircle_tmp/FlankingRead_80_twice_intra.list | awk '{print $3}' | sort  | uniq > $out/WithinCircle_tmp/circ_twice_intra.list

echo -n > $out/WithinCircle_tmp/circ.split.bed

cat $out/WithinCircle_tmp/circ_twice_intra.list | while read CircOne
do
  join -t$'\t' -1 3 -2 1 -o 1.1,1.2,1.3 $out/WithinCircle_tmp/FlankingRead_80_twice_intra.list <(echo $CircOne) > $out/WithinCircle_tmp/CircOne_read.list

  CircOneNum=$(cat $out/WithinCircle_tmp/CircOne_read.list | wc -l )
  for readi in $(seq 1 1 $CircOneNum)
  do 
       one=$(cat $out/WithinCircle_tmp/CircOne_read.list | head -n $readi | tail -n 1)
       ReadOne=$(echo $one | awk '{print $1}')
       
       intraChr=$(echo $one | awk '{print $3}' | sed 's/:/\t/g' | awk '{print $1}')
       intraDonor=$(echo $one | awk '{print $3}' | sed 's/:/\t/g' | awk '{print $2}')
       intraAcceptor=$(echo $one | awk '{print $3}'| sed 's/:/\t/g' | awk '{print $5}')
       intraMax=$(echo $intraDonor"\n"$intraAcceptor | sort -r | head -n 1)
       intraMin=$(echo $intraDonor"\n"$intraAcceptor | sort | head -n 1)

       join -t$'\t' $out/WithinCircle_tmp/FlankingRead_80_twice_intra.paf <(echo $ReadOne) | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' | sort -nk3 | awk '{print $0 "\t" int(($4-$3+1)/2)}' | awk '{print $0 "\t" $3+$5}' > $out/WithinCircle_tmp/$CircOne.split.tmp1

       fragmentNum=$(cat $out/WithinCircle_tmp/$CircOne.split.tmp1 | wc -l)
       echo -n > $out/WithinCircle_tmp/$CircOne.split.bed    
     
       #fragmentNum=1
       cat $out/WithinCircle_tmp/$CircOne.split.tmp1 | head -n 1 | awk '{print $1 "\t" 0 "\t" $6}' >> $out/WithinCircle_tmp/$CircOne.split.bed
     
       #fragmentNum=2~(n-1)
       for fragmentN in $(seq 2 1 $fragmentNum)
       do
           fragmentNminus=$(($fragmentN-1))
           start_fragment=$(cat $out/WithinCircle_tmp/$CircOne.split.tmp1 | head -n $fragmentNminus | tail -n 1 | awk '{print $6}')
           cat $out/WithinCircle_tmp/$CircOne.split.tmp1 | head -n $fragmentN | tail -n 1 | awk '{print $1 "\t" sg "\t" $6}' sg=$start_fragment >> $out/WithinCircle_tmp/$CircOne.split.bed
       done
     
       #fragmentNum=n
       start_fragment=$(cat $out/WithinCircle_tmp/$CircOne.split.tmp1 | head -n $fragmentNum | tail -n 1 | awk '{print $6}')
       cat $out/WithinCircle_tmp/$CircOne.split.tmp1 | head -n $fragmentNum | tail -n 1 | awk '{print $1 "\t" sg "\t" $2}' sg=$start_fragment >> $out/WithinCircle_tmp/$CircOne.split.bed

       cat $out/WithinCircle_tmp/$CircOne.split.bed >> $out/WithinCircle_tmp/circ.split.bed
   done
done
   
$seqtk_link subseq $out/tmp/All.fa $out/WithinCircle_tmp/circ.split.bed > $out/WithinCircle_tmp/circ.split.fa

echo "Step1: to align split sequences of a long read against whole genome"
$minimap2_link -t 10 -ax splice $genome_prefix.mmi $out/WithinCircle_tmp/circ.split.fa | $samtools_link view -bS - > $out/WithinCircle_tmp/circ.bam
$bedtools_link bamtobed -bed12 -i $out/WithinCircle_tmp/circ.bam | awk '$1~/^chr[0-9XY]/'| awk '$5==60'| sort -k4  > $out/WithinCircle_tmp/circ.tmp.bed12

cat $out/WithinCircle_tmp/circ.tmp.bed12 | awk -F'\t' '{print $4"\t"$0}' | awk 'BEGIN{FS=OFS="\t"}{sub(/:.*/, "", $1); print $0}' | sort -k1,1 > $out/WithinCircle_tmp/circ.tmp.bed12.with_read_id

echo "Step2: to report a circRNA long read"
cat $out/WithinCircle_tmp/circ_twice_intra.list | while read CircOne
do
    join -t$'\t' -1 3 -2 1 -o 1.1,1.2,1.3 $out/WithinCircle_tmp/FlankingRead_80_twice_intra.list <(echo $CircOne) > $out/WithinCircle_tmp/CircOne_read.list

    CircOneNum=$(cat $out/WithinCircle_tmp/CircOne_read.list | wc -l )
    for readi in $(seq 1 1 $CircOneNum)
    do 
       one=$(cat $out/WithinCircle_tmp/CircOne_read.list | head -n $readi | tail -n 1)
       # echo $one
       ReadOne=$(echo $one | awk '{print $1}')
       
       intraChr=$(echo $one | awk '{print $3}' | sed 's/:/\t/g' | awk '{print $1}')
       intraDonor=$(echo $one | awk '{print $3}' | sed 's/:/\t/g' | awk '{print $2}')
       intraAcceptor=$(echo $one | awk '{print $3}'| sed 's/:/\t/g' | awk '{print $5}')
       intraMax=$(echo $intraDonor"\n"$intraAcceptor | sort -r | head -n 1)
       intraMin=$(echo $intraDonor"\n"$intraAcceptor | sort | head -n 1)

       join -t$'\t' $out/WithinCircle_tmp/circ.tmp.bed12.with_read_id <(echo $ReadOne) | cut -f '2-' > $out/WithinCircle_tmp/$CircOne.tmp.bed12

       bedNum=$(cat $out/WithinCircle_tmp/$CircOne.tmp.bed12 | awk '{print $4}' | wc -l)
       
       chrNum=$(cat $out/WithinCircle_tmp/$CircOne.tmp.bed12 | awk '{print $1}' | uniq | wc -l) 
       chrSplit=$(cat $out/WithinCircle_tmp/$CircOne.tmp.bed12 | awk '{print $1}' | uniq | head -n 1)  
       startSplit=$(cat $out/WithinCircle_tmp/$CircOne.tmp.bed12 | awk '{print $2 "\n" $3}' | sort | head -n 1)
       endSplit=$(cat $out/WithinCircle_tmp/$CircOne.tmp.bed12 | awk '{print $2 "\n" $3}' | sort -r | head -n 1)

       if [ "$bedNum" -eq $(($fragmentNum+1)) ] && [ "$chrNum" -eq 1 ] && [ $(($intraMin-10)) -le "$startSplit" ] && [ "$endSplit" -le $(($intraMax+10)) ] && [ "$intraChr"=="$chrSplit" ]
       then 
             cat $out/WithinCircle_tmp/$CircOne.tmp.bed12 > $out/WithinCircle_tmp/WithinCircle_reads_tmp/$CircOne.$readi.bed12
             cat $out/WithinCircle_tmp/$CircOne.tmp.bed12 >> $out/WithinCircle_events/$CircOne.bed12
       fi 
    done    
done

ls $out/WithinCircle_tmp/WithinCircle_reads_tmp/ | sed 's/\./\t/' | awk '{print $1}' | sort | uniq -c | awk '{print $2  "\t" "1" "\t" $1}' > $out/$out\_WithinCircle.result

