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
     -t | --threads)
     threads=$2
     shift
     ;;
     *)

esac
shift
done

if [[ -z "$input" ]]; then
   echo ""
   echo "Usage:"
   echo "./WithinCircle.sh -input_folder [input folder] -c [configure file link] -o [output_prefix_name] -t [num_threads]"
   echo ""
   exit
fi


if [[ -z "$out" ]]; then
   echo ""
   echo "Usage:"
   echo "./WithinCircle.sh -input_folder [input folder] -c [configure file link] -o [output_prefix_name] -t [num_threads]"
   echo ""
   exit
fi

if [[ -z "$threads" ]]; then
   threads=1
fi

source $config

## Within Circle : obtain more than one pseudo-reference (200 bp NCL sequence) 
echo $input
echo $out
echo -n > $out/$out\_WithinCircle.result  

mkdir -p \
    $out/WithinCircle_tmp \
    $out/WithinCircle_events \
    $out/WithinCircle_tmp/WithinCircle_reads_tmp


cat $input | awk '{print $1 "\t" $5 "\t" $6}' | sort -k3,3 -k1,1 -k2,2 | uniq -c | awk '$1 > 1' | sed 's/:/\t/g' | awk '$4==$7' | awk '{print $2 "\t" $3 "\t" $4":"$5":"$6":"$7":"$8":"$9}' > $out/WithinCircle_tmp/FlankingRead_80_twice.list
cat $out/WithinCircle_tmp/FlankingRead_80_twice.list > $out/WithinCircle_tmp/FlankingRead_80_twice_intra.list
cat $out/WithinCircle_tmp/FlankingRead_80_twice_intra.list | tr '\t' '_' | sort > $out/WithinCircle_tmp/FlankingRead_80_twice_intra.key

cat $input | awk '{print $1"_"$5"_"$6"\t"$0}' | sort -k1,1 | sort -k1,1 > $out/WithinCircle_tmp/FlankingRead_80.paf.with_key
join -t$'\t' $out/WithinCircle_tmp/FlankingRead_80.paf.with_key $out/WithinCircle_tmp/FlankingRead_80_twice_intra.key | cut -f'2-' > $out/WithinCircle_tmp/FlankingRead_80_twice_intra.paf

cat $out/WithinCircle_tmp/FlankingRead_80_twice_intra.list | awk '{print $3}' | sort  | uniq > $out/WithinCircle_tmp/circ_twice_intra.list


cat $out/WithinCircle_tmp/FlankingRead_80_twice_intra.paf \
    | awk -F'\t' '{print $1"\t"$2"\t"$NF"\t"$6}' \
    | sort -k4,4 -k1,1 -k3,3n \
    | awk 'BEGIN{
            FS=OFS="\t";
            eventID="";
            readID="";
            readLen=0;
            start=0;
        }{
            if(!(($1==readID) && ($4==eventID))){
                print readID, start, readLen; 
                eventID=$4;
                readID=$1; 
                readLen=$2; 
                start=0;
            }; 
            print readID, start, $3;
            start=$3
        }
        END{
            print readID, start, readLen
        }' \
    | sed '1d' > $out/WithinCircle_tmp/circ.split.bed

   
$seqtk_link subseq $out/tmp/All.fa $out/WithinCircle_tmp/circ.split.bed > $out/WithinCircle_tmp/circ.split.fa

echo "Step1: to align split sequences of a long read against whole genome"

if [ ! -e "$genome_prefix.mmi" ]; then
   $minimap2_link -d $genome_prefix.mmi $genome_prefix.fa
fi

#$minimap2_link -t $threads -ax splice $genome_prefix.mmi $out/WithinCircle_tmp/circ.split.fa | $samtools_link view -bS - > $out/WithinCircle_tmp/circ.bam
#$bedtools_link bamtobed -bed12 -i $out/WithinCircle_tmp/circ.bam | awk '$1~/^chr[0-9XY]/'| awk '$5==60'| sort -k4,4  > $out/WithinCircle_tmp/circ.tmp.bed12

cat $out/WithinCircle_tmp/circ.tmp.bed12 | awk -F'\t' '{print $4"\t"$0}' | awk 'BEGIN{FS=OFS="\t"}{sub(/:.*/, "", $1); print $0}' | sort -k1,1 > $out/WithinCircle_tmp/circ.tmp.bed12.with_read_id

echo "Step2: to report a circRNA long read"

cat $out/WithinCircle_tmp/FlankingRead_80_twice_intra.list | awk -F'\t' '{print $3"\t"$0}' | sort -k1,1 -k2,2 -k3,3 | $NCLscan_hybrid_bin/split_file_by_first_column.py - -o $out/WithinCircle_tmp/WithinCircle_list_tmp -s .list

join -t$'\t' $out/WithinCircle_tmp/circ.tmp.bed12.with_read_id <(cat $out/WithinCircle_tmp/FlankingRead_80_twice_intra.list | sort -k1,1) | sort -k15,15 -k1,1 -k5,5 | awk -F'\t' '{print $15"\t"$0}' | cut -f '-14' | $NCLscan_hybrid_bin/split_file_by_first_column.py - -o $out/WithinCircle_tmp/WithinCircle_bed_tmp -s .bed12.tmp
touch $(cat $out/WithinCircle_tmp/circ_twice_intra.list | awk '{print dir"/"$1".bed12.tmp"}' dir=$out/WithinCircle_tmp/WithinCircle_bed_tmp)

cat $out/WithinCircle_tmp/FlankingRead_80_twice_intra.paf | cut -f '1,5,6' | sort -k3,3 -k1,1 -k2,2 | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$1}' | awk 'BEGIN{FS=OFS="\t";event=""; idx=0}{if($3!=event){event=$3; idx=1;}; print $0,idx; idx+=1;}' > $out/WithinCircle_tmp/FlankingRead_80_twice_intra.paf.fragmentNum.readi

cat $out/WithinCircle_tmp/FlankingRead_80_twice_intra.paf.fragmentNum.readi  | while read ReadOne strand CircOne fragmentNum readi
do
  #echo $ReadOne
  #echo $strand
  #echo $CircOne
  #echo $fragmentNum
  #echo $readi
  #ReadOne="m130427_012547_42142_c100503071630000001823073108081320_s1_p0/145009/5756_6733"
  #strand="+"
  #CircOne="chr17:16382130:+:chr17:16381902:+"
  #fragmentNum=2
  #readi=2
 
  intraChr=$(echo $CircOne | cut -d':' -f'1')
  intraDonor=$(echo $CircOne | cut -d':' -f'2')
  intraAcceptor=$(echo $CircOne | cut -d':' -f'5')
  intraMax=$(echo $intraDonor"\n"$intraAcceptor | sort -r | head -n 1)
  intraMin=$(echo $intraDonor"\n"$intraAcceptor | sort | head -n 1)
  
  join -t$'\t' $out/WithinCircle_tmp/WithinCircle_bed_tmp/$CircOne.bed12.tmp <(echo $ReadOne) | cut -f '2-' > $out/WithinCircle_tmp/$CircOne.$readi.tmp.bed12
        
  ##intra boundaries10bp
  echo $intraChr"\t"$(($intraMin-10))"\t"$(($intraMin+10))"\t""downstream""\n"$intraChr"\t"$(($intraMax-10))"\t"$(($intraMax+10))"\t""upstream" >  $out/WithinCircle_tmp/intraBoundary_10bp.bed
  if [ -s $out/WithinCircle_tmp/$CircOne.$readi.tmp.bed12 ]
  then
     print "Yes1"
     split_start=$(cat $out/WithinCircle_tmp/$CircOne.$readi.tmp.bed12 | sed  's/:/\t/' | sed 's/-/\t/' | sort -k5,5n | head -n 1) 
     split_end=$(cat $out/WithinCircle_tmp/$CircOne.$readi.tmp.bed12 | sed  's/:/\t/' | sed 's/-/\t/' | sort -k5,5n | tail -n 1) 
     echo $split_start | awk '{print $1 "\t" $2-1 "\t" $2 "\n" $1 "\t" $3-1 "\t" $3}' > $out/WithinCircle_tmp/split_start_points.bed
     echo $split_end | awk '{print $1 "\t" $2-1 "\t" $2 "\n" $1 "\t" $3-1 "\t" $3}' > $out/WithinCircle_tmp/split_end_points.bed
     check_split_start_num=$($bedtools_link intersect -a $out/WithinCircle_tmp/intraBoundary_10bp.bed -b $out/WithinCircle_tmp/split_start_points.bed -wa -wb  | wc -l)
     check_split_end_num=$($bedtools_link intersect -a $out/WithinCircle_tmp/intraBoundary_10bp.bed -b $out/WithinCircle_tmp/split_end_points.bed -wa -wb | wc -l)
     if [ "$check_split_start_num" -gt 0  ] && [ "$check_split_end_num" -gt 0 ]
     then
         print "Yes2"
         check_split_start=$(echo $check_split_start_1 | awk '{print $4}' | sort | uniq | head -n 1)
         check_split_end1=$(echo $check_split_end_1 | awk '{print $4}' | sort | uniq | head -n 1)
         check_split_end2=$(echo $check_split_end_1 | awk '{print $4}' | sort -r | uniq | head -n 1)
         ##one side of first split fragment and the last split fragment locate within intra boundaries 10bp
         if [ "$check_split_start"!="$check_split_end1" ] || [ "$check_split_start"!="$check_split_end2" ] 
         then
             print "Yes3"
             ##two sides of between split fragments locate within intra boundaries 10bp
             split_bewteen=$(cat $out/WithinCircle_tmp/$CircOne.$readi.tmp.bed12 | sed 's/:/\t/' | sed 's/-/\t/' | sort -k5,5n | sed -e '1d' | sed -e '$d')
             split_bewteen_num=$(cat $out/WithinCircle_tmp/$CircOne.$readi.tmp.bed12 | sed 's/:/\t/' | sed 's/-/\t/' | sort -k5,5n | sed -e '1d' | sed -e '$d' | wc -l)
             if [ $split_bewteen_num -gt 0 ]
             then
                   print "Yes4"
                   echo -n > $out/WithinCircle_tmp/split_between_count.txt
                   echo $split_bewteen | while read one
                   do
                       echo $one | awk '{print $1 "\t" $2-1 "\t" $2 "\n" $1 "\t" $3-1 "\t" $3}' > $out/WithinCircle_tmp/split_between_points.bed
                       $bedtools_link intersect -a $out/WithinCircle_tmp/intraBoundary_10bp.bed -b $out/WithinCircle_tmp/split_between_points.bed -wa -wb \
                      | awk '{print $4}' | sort | uniq | wc -l \
                      | awk '{print $1 "\t" one}' one=$one >> $out/WithinCircle_tmp/split_between_count.txt
                   done
                   check_split_between=$(cat $out/WithinCircle_tmp/split_between_count.txt | awk '$1!=2')
                   if [ -s $check_split_bwtween ]
                     then
                     print "Yes5"
                      ##check pseudo refernece number (fragment number) + 1 equal to bed numbers
                      ##check side positions of splits witin intra region with 10 bp buffer 
                      bedNum1=$(cat $out/WithinCircle_tmp/$CircOne.$readi.tmp.bed12 | awk '{print $4}' | sort | uniq | wc -l)
                      bedNum2=$(cat $out/WithinCircle_tmp/$CircOne.$readi.tmp.bed12 | awk '{print $4}' | wc -l)
                      chrNum=$(cat $out/WithinCircle_tmp/$CircOne.$readi.tmp.bed12 | awk '{print $1}' | uniq | wc -l) 
                      chrSplit=$(cat $out/WithinCircle_tmp/$CircOne.$readi.tmp.bed12 | awk '{print $1}' | uniq | head -n 1)  
                      startSplit=$(cat $out/WithinCircle_tmp/$CircOne.$readi.tmp.bed12 | awk '{print $2 "\n" $3}' | sort | head -n 1)
                      endSplit=$(cat $out/WithinCircle_tmp/$CircOne.$readi.tmp.bed12 | awk '{print $2 "\n" $3}' | sort -r | head -n 1)
                       if [ "$bedNum1" -eq "$bedNum2" ] && [ "$bedNum2" -eq $(($fragmentNum+1)) ] && [ "$chrNum" -eq 1 ] && [ $(($intraMin-10)) -le "$startSplit" ] && [ "$endSplit" -le $(($intraMax+10)) ] && [ "$intraChr"=="$chrSplit" ]
                         then 
                          print "Yes6"
                          cat $out/WithinCircle_tmp/$CircOne.$readi.tmp.bed12 > $out/WithinCircle_tmp/WithinCircle_reads_tmp/$CircOne.$readi.bed12
                          cat $out/WithinCircle_tmp/$CircOne.$readi.tmp.bed12 >> $out/WithinCircle_events/$CircOne.bed12
                       fi
                   fi  
               fi   
              
           fi
       fi 
   fi
  
done

ls $out/WithinCircle_tmp/WithinCircle_reads_tmp/ | sed 's/\./\t/' | awk '{print $1}' | sort | uniq -c | awk '{print $2  "\t" "1" "\t" $1}' > $out/$out\_WithinCircle.result