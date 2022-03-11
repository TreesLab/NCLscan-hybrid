#!/usr/bin/env zsh

while [[ $# -gt 1 ]]
do 
key="$1"

case $key in 
     -long | --longread)
     longread=$2
     shift
     ;;
     -long_type | --longread_type)
     long_type=$2
     shift
     ;;
     -nclscan | --NCLscan)
     NCLscan=$2
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
      
if [[ -z "$longread" ]]; then
   echo ""
   echo "Usage:"
   echo "./NCLscan-hybrid.sh -long [input long read fasta file] -long_type [pb or ont] -nclscan [NCLscan result file] -c [configure file link] -o [output_prefix]"
   echo ""
   exit
fi


source $config
# BASEDIR=$(pwd)


rm -rf $out

mkdir -p \
    $out \
    $out/tmp \
    $out/pass1 \
    $out/fail1 \
    $out/pass2_intra \
    $out/pass2_inter 


echo -n > $out/$out\_long_intra.result
echo -n > $out/$out\_long_inter.result

echo "Step: to create flankingSeqs of NCL events" 

$NCLscan_hybrid_bin/FlankingSeq.sh \
   -nclscan $NCLscan \
   -fl 100 \
   -gtf $gtf \
   -g $genome_prefix.fa \
   -bedtools $bedtools_link  \
   -mps $NCLscan_hybrid_bin/merge_paired_sequences.py \
   -o $out/$out\_100bp

echo "Step: to align long reads against flankingSeqs"
# $minimap2_link -d $out/$out\_100bp_flanking_merged.mmi $out/$out\_100bp_flanking_merged.fa   
$minimap2_link -t 10 -x map-$long_type $out/$out\_100bp_flanking_merged.fa $longread -c --secondary=no > $out/tmp/$out\_to_FlankingRead.paf

cat $out/tmp/$out\_to_FlankingRead.paf | sort -k6,6 -k1,1 -k3,3n -k4,4n > $out/tmp/$out\_to_FlankingRead.sorted.paf
cat $out/tmp/$out\_to_FlankingRead.sorted.paf | awk -F'\t' '{print $6"\t"$0}' > $out/tmp/$out\_to_FlankingRead.sorted.paf.with_id.tmp
join -t$'\t' $out/tmp/$out\_to_FlankingRead.sorted.paf.with_id.tmp $out/$out\_100bp_NCL_flanking.length | awk -F'\t' '($10 - $(NF-1) > 10) && ($(NF-1) - $9 > 10)' | awk 'BEGIN{FS=OFS="\t"}{NF-=1; print $0}' | cut -f '2-' > $out/tmp/$out\_to_FlankingRead_overhang10.paf

cat $out/tmp/$out\_to_FlankingRead_overhang10.paf | awk -F'\t' '$11/$7 > 0.8 && $11/$10 < 2 && $12 == 60' > $out/tmp/$out\_to_FlankingRead_80.paf
cat $out/tmp/$out\_to_FlankingRead_80.paf | cut -f '6' | sort | uniq  > $out/tmp/$out\_to_FlankingRead_80.list

cat $out/tmp/$out\_to_FlankingRead_80.paf | awk -F'\t' '{print $0"\t"$(NF)-$8}' | awk -F'\t' '{print $(NF-2)"\t"$NF}' | sed 's/^cg:Z://g' | $NCLscan_hybrid_bin/get_qpos.py -f - > $out/tmp/$out\_to_FlankingRead_80.paf.qpos
paste $out/tmp/$out\_to_FlankingRead_80.paf $out/tmp/$out\_to_FlankingRead_80.paf.qpos | awk -F'\t' '($5=="+"){print $0"\t"$3+$NF}($5=="-"){print $0"\t"$4-$NF}' > $out/tmp/$out\_to_FlankingRead_80.junction.paf


# ----- Out of Circle: start ----- #

cat $out/tmp/$out\_to_FlankingRead_80.junction.paf | awk -F'\t' '{print $6"\t"$0}' > $out/tmp/$out\_to_FlankingRead_80.paf.with_id.tmp

echo -n > $out/tmp/All.list
echo -n > $out/tmp/split.bed
cat $out/tmp/$out\_to_FlankingRead_80.list | while read one
do 
   join -t$'\t' $out/tmp/$out\_to_FlankingRead_80.paf.with_id.tmp <(echo $one) | cut -f '2-' > $out/tmp/$one.paf
   cat $out/tmp/$one.paf | awk '{print $1}' | sort | uniq > $out/tmp/$one.list 

   cat $out/tmp/$one.paf | awk -F'\t' '{print $1"\t"0"\t"$NF"\n"$1"\t"$NF"\t"$2}' > $out/tmp/$one.split.bed

   cat $out/tmp/$one.list >> $out/tmp/All.list
   cat $out/tmp/$one.split.bed >> $out/tmp/split.bed
done

$seqtk_link subseq $longread $out/tmp/All.list > $out/tmp/All.fa
$seqtk_link subseq $out/tmp/All.fa $out/tmp/split.bed > $out/tmp/All_split.fa

echo "Step: to align splited reads against whole genome"
$minimap2_link -d $genome_prefix.mmi $genome_prefix.fa
$minimap2_link -t 10 -ax splice $genome_prefix.mmi $out/tmp/All_split.fa | $samtools_link view -bS - > $out/tmp/All.bam

$bedtools_link bamtobed -bed12 -i $out/tmp/All.bam | awk '$1~/^chr[0-9XY]/'| awk '$5>0'| sort -k4,4 > $out/tmp/All.bed12

cat $out/tmp/All.bed12 | awk -F'\t' '{print $4"\t"$0}' | sort -k1,1 > $out/tmp/All.bed12.sorted
cat $out/tmp/All.bed12.sorted | cut -f'1' | uniq -c | awk '$1==1 {print $2}'  > $out/tmp/All.bed12.uniq.list
join -t$'\t' $out/tmp/All.bed12.sorted $out/tmp/All.bed12.uniq.list | cut -f'2-' > $out/tmp/All.uniq.bed12


echo "Step: to diagnose split reads mapped to donor and acceptor: uniquely mapping and same chromsome"

cat $out/tmp/All.uniq.bed12 | awk -F'\t' '{print $4"\t"$0}' | awk 'BEGIN{FS="\t";OFS="\t"}{sub(/:.*/, "", $1); print $0}' | sort -k1,1 -k5,5 > $out/tmp/All.uniq.read_ID.bed12

cat $out/tmp/$out\_to_FlankingRead_80.list | while read one
do 
   join -t$'\t' $out/tmp/All.uniq.read_ID.bed12 $out/tmp/$one.list | cut -f '2-' > $out/tmp/$one.bed12

   chr_donor=$(echo $one | sed 's/:/\t/g' | awk '{print $1}')  
   chr_acceptor=$(echo $one | sed 's/:/\t/g' | awk '{print $4}')

   donor_count=$(cat $out/tmp/$one.bed12 | awk '$1==donor' donor=$chr_donor | wc -l)  
   acceptor_count=$(cat $out/tmp/$one.bed12 | awk '$1==acceptor' acceptor=$chr_acceptor | wc -l)  

   if [[ "$donor_count" -gt 0 ]] && [[ "$acceptor_count" -gt 0 ]]
      then mv $out/tmp/$one.bed12  $out/pass1/$one.bed12
      else mv $out/tmp/$one.bed12  $out/fail1/$one.bed12
   fi
done

ls $out/pass1 | sed 's/\./\t/g' | awk '{print $1}' | sort -k1,1 > $out/tmp/pass1.list

cat $NCLscan | awk '{print $1":"$2":"$3":"$4":"$5":"$6 "\t" $0}' | sort -k1,1 > $out/tmp/result.tmp
join $out/tmp/$out\_to_FlankingRead_80.list $out/tmp/result.tmp | tr ' ' \\t > $out/tmp/result.tmp2
cat $out/tmp/result.tmp2 | awk '$10==1 {print $1}' | sort > $out/tmp/intra.list
cat $out/tmp/result.tmp2 | awk '$10==0 {print $1}' | sort > $out/tmp/inter.list
join $out/tmp/intra.list $out/tmp/pass1.list  > $out/tmp/pass1_intra.list
join $out/tmp/inter.list $out/tmp/pass1.list  > $out/tmp/pass1_inter.list


echo "Step: at least one supported long read"
## intra ##
cat $out/tmp/pass1_intra.list | while read one
do
   chr=$(echo $one | sed 's/:/\t/g' | awk '{print $1}')  
   cp $out/pass1/$one.bed12  $out/pass2_intra/$one.bed12.tmp
   
   cat $out/pass2_intra/$one.bed12.tmp | awk '$5>=60' | cut -f '4' | sed 's/[:-]/\t/g' > $out/pieces.tmp
   cat $out/pieces.tmp | awk -F'\t' '$2==1{print $1":"$3+1"\t"$2"-"$3}' | sort -k1,1 > $out/pieces1.tmp1
   cat $out/pieces.tmp | awk -F'\t' '$2!=1{print $1":"$2"\t"$2"-"$3}' | sort -k1,1 > $out/pieces2.tmp1
   join -t$'\t' $out/pieces1.tmp1 $out/pieces2.tmp1 > $out/pieces12.tmp1

   cat $out/pieces12.tmp1  | sed 's/:/\t/g' | awk '{print $1":"$3}' > $out/pieces1.tmp2
   cat $out/pieces12.tmp1  | sed 's/:/\t/g' | awk '{print $1":"$4}' > $out/pieces2.tmp2
   cat $out/pieces1.tmp2 $out/pieces2.tmp2 | sort > $out/long.tmp1
   join -t$'\t' $out/pass2_intra/$one.bed12.tmp $out/long.tmp1 -1 4 -2 1 | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}' > $out/long.tmp2
   
   count_long=$(cat $out/long.tmp2 | wc -l)

   if [[ "$count_long" -gt 0 ]]; then
  
      cat $out/long.tmp2 | awk '$5==60' > $out/pass2_intra/$one.bed12
   fi

   rm -r -f $out/pass2_intra/$one.bed12.tmp
done

echo "Step: to check_intra_OutOfCircle"
## Out of circle: long read extends more than 100 bp on one side (upstream or downstream)  
$NCLscan_hybrid_bin/OutOfCircle.sh -input_folder $out/pass2_intra -o $out/$out -L 100

# ----- Out of Circle: end ----- #


echo  "Step: to check_intra_WithinCircle"
## Within Circle : obtain more than one pseudo-reference (200 bp NCL sequence) 
$NCLscan_hybrid_bin/WithinCircle.sh -input_folder $out/tmp/$out\_to_FlankingRead_80.junction.paf -c $config -o $out
##Within circle _view##
$NCLscan_hybrid_bin/BrowserView.sh -input_folder $out/WithinCircle_events




##intra_view##
$NCLscan_hybrid_bin/BrowserView.sh -input_folder $out/pass2_intra


## inter ##
cat $out/tmp/pass1_inter.list | while read one
do
   cp $out/pass1/$one.bed12 $out/pass2_inter/$one.bed12.tmp

   cat $out/pass2_inter/$one.bed12.tmp | awk '$5==60' | cut -f '4' | sed 's/[:-]/\t/g' > $out/pieces.tmp
   cat $out/pieces.tmp | awk -F'\t' '$2==1{print $1":"$3+1"\t"$2"-"$3}' | sort -k1,1 > $out/pieces1.tmp1
   cat $out/pieces.tmp | awk -F'\t' '$2!=1{print $1":"$2"\t"$2"-"$3}' | sort -k1,1 > $out/pieces2.tmp1
   join -t$'\t' $out/pieces1.tmp1 $out/pieces2.tmp1 > $out/pieces12.tmp1

   cat $out/pieces12.tmp1  | sed 's/:/\t/g' | awk '{print $1":"$3}' > $out/pieces1.tmp2
   cat $out/pieces12.tmp1  | sed 's/:/\t/g' | awk '{print $1":"$4}' > $out/pieces2.tmp2
   cat $out/pieces1.tmp2 $out/pieces2.tmp2 | sort > $out/long.tmp1
   join -t$'\t' $out/pass2_inter/$one.bed12.tmp $out/long.tmp1 -1 4 -2 1 | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}' > $out/long.tmp2

   count_long=$(cat $out/long.tmp2 | wc -l)

   if [[ "$count_long" -gt 0 ]]; then
  
      cat $out/long.tmp2 | awk '$5>=60' > $out/pass2_intra/$one.bed12
   fi

   rm -r -f $out/pass2_intra/$one.bed12.tmp

done

##inter_view##
$NCLscan_hybrid_bin/BrowserView.sh -input_folder $out/pass2_inter

rm -rf \
   $out/pieces.tmp \
   $out/piece1.tmp1 \
   $out/piece2.tmp1 \
   $out/piece1.tmp2 \
   $out/piece2.tmp2 \
   $out/piece12.tmp1 \
   $out/long.tmp1 \
   $out/long.tmp2


echo "Step: output intra result"

echo -n > $out/tmp/NCL_long_intra.tmp1

ls $out/pass2_intra > $out/tmp/pass2_intra.list
cat $out/tmp/pass2_intra.list  | while read one
do  
   cat $out/pass2_intra/$one | awk '{print $4}' | sed 's/:/\t/g' | awk 'NF{NF-=1};1' | sort | uniq | wc -l | awk '{print event "\t" $1}' event=$one | sed 's/\.bed12//g' >> $out/tmp/NCL_long_intra.tmp1
done 

cat $out/WithinCircle_tmp/FlankingRead_80_twice_intra.list | awk '{print $3}' | sort | uniq -c | awk '{print $2 "\t" $1}' |  sort -k1 > $out/WithinCircle_tmp/FlankingRead_80_twice_intra.tmp1
cat $out/tmp/NCL_long_intra.tmp1 | sort -k1 > $out/tmp/NCL_long_intra.tmp2
cat $out/tmp/NCL_long_intra.tmp2 $out/WithinCircle_tmp/FlankingRead_80_twice_intra.tmp1  | awk '{print $1}' | sort -k1 | uniq | sort -k1 > $out/tmp/NCL_long_intra.tmp3

join -o 1.1 2.2 $out/tmp/NCL_long_intra.tmp3 $out/tmp/NCL_long_intra.tmp2  -a1 -e0 > $out/tmp/NCL_long_intra.tmp4
join -o 1.1 1.2 2.2 $out/tmp/NCL_long_intra.tmp4 $out/WithinCircle_tmp/FlankingRead_80_twice_intra.tmp1 -a1 -e0 > $out/tmp/NCL_long_intra.tmp5
cat $out/tmp/NCL_long_intra.tmp5 | awk '{print $1 "\t" $2+$3}' > $out/tmp/NCL_long_intra.tmp6
join -o 1.1 1.2 2.2 2.3 $out/tmp/NCL_long_intra.tmp6 $out/$out\_OC.result -a1 -e0 > $out/tmp/NCL_long_intra.tmp7
join -o 1.1 1.2 1.3 1.4 2.2 2.3 $out/tmp/NCL_long_intra.tmp7 $out/$out\_WithinCircle.result -a1 -e0 > $out/tmp/NCL_long_intra.OC.WithinCircle
join $out/tmp/NCL_long_intra.OC.WithinCircle $out/tmp/result.tmp | tr ' ' \\t | sort | uniq > $out/$out\_long_intra.result

echo "Step: output inter result"

echo -n > $out/tmp/NCL_long_inter.tmp1

ls $out/pass2_inter > $out/tmp/pass2_inter.list  
cat $out/tmp/pass2_inter.list | while read one
do  
   cat $out/pass2_inter/$one | awk '{print $4}' | sed 's/:/\t/g' | awk 'NF{NF-=1};1' | sort | uniq | wc -l | awk '{print event "\t" $1}' event=$one | sed 's/\.bed12//g' >> $out/tmp/NCL_long_inter.tmp1
done 
cat $out/tmp/NCL_long_inter.tmp1 | sort -k1 > $out/tmp/NCL_long_inter.tmp2
join $out/tmp/NCL_long_inter.tmp2 $out/tmp/result.tmp | tr ' ' \\t | sort | uniq  > $out/$out\_long_inter.result

