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
     -t | --threads)
     threads=$2
     shift
     ;;
     *)

esac
shift
done
      
if [[ -z "$longread" ]]; then
   echo ""
   echo "Usage:"
   echo "./NCLscan-hybrid.sh -long [input long read fasta file] -long_type [pb or ont] -nclscan [NCLscan result file] -c [configure file link] -o [output_prefix_name] -t [num_threads]"
   echo ""
   exit
fi

if [[ -z "$threads" ]]; then
   threads=1
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

$NCLscan_hybrid_bin/FlankingSeq_exonic.sh \
   -nclscan $NCLscan \
   -fl 100 \
   -gtf $gtf \
   -g $genome_prefix.fa \
   -bedtools $bedtools_link  \
   -mps $NCLscan_hybrid_bin/merge_paired_sequences.py \
   -o $out/$out\_100bp

echo "Step: to align long reads against flankingSeqs"
$minimap2_link -d $out/$out\_100bp_flanking_merged.mmi $out/$out\_100bp_flanking_merged.fa   
$minimap2_link -t $threads -x map-$long_type $out/$out\_100bp_flanking_merged.fa $longread -c --secondary=no > $out/tmp/$out\_to_FlankingRead.paf

cat $out/tmp/$out\_to_FlankingRead.paf | sort -k6,6 -k1,1 -k3,3n -k4,4n > $out/tmp/$out\_to_FlankingRead.sorted.paf
cat $out/tmp/$out\_to_FlankingRead.sorted.paf | awk -F'\t' '{print $6"\t"$0}' > $out/tmp/$out\_to_FlankingRead.sorted.paf.with_id.tmp
join -t$'\t' $out/tmp/$out\_to_FlankingRead.sorted.paf.with_id.tmp $out/$out\_100bp_NCL_flanking.length | awk -F'\t' '($10 - $(NF-1) > 10) && ($(NF-1) - $9 > 10)' | awk 'BEGIN{FS=OFS="\t"}{NF-=1; print $0}' | cut -f '2-' > $out/tmp/$out\_to_FlankingRead_overhang10.paf

cat $out/tmp/$out\_to_FlankingRead_overhang10.paf | awk -F'\t' '$11/$7 > 0.8 && $11/$10 < 2 && $12 == 60' > $out/tmp/$out\_to_FlankingRead_80.paf
cat $out/tmp/$out\_to_FlankingRead_80.paf | cut -f '6' | sort | uniq  > $out/tmp/$out\_to_FlankingRead_80.list

cat $out/tmp/$out\_to_FlankingRead_80.paf | awk -F'\t' '{print $0"\t"$(NF)-$8}' | awk -F'\t' '{print $(NF-2)"\t"$NF}' | sed 's/^cg:Z://g' | $NCLscan_hybrid_bin/get_qpos.py -f - > $out/tmp/$out\_to_FlankingRead_80.paf.qpos
paste $out/tmp/$out\_to_FlankingRead_80.paf $out/tmp/$out\_to_FlankingRead_80.paf.qpos | awk -F'\t' '($5=="+"){print $0"\t"$3+$NF}($5=="-"){print $0"\t"$4-$NF}' > $out/tmp/$out\_to_FlankingRead_80.junction.paf


cat $out/tmp/$out\_to_FlankingRead_80.junction.paf | awk -F'\t' '{print $6"\t"$0}' > $out/tmp/$out\_to_FlankingRead_80.paf.with_id.tmp

cat $out/tmp/$out\_to_FlankingRead_80.junction.paf | cut -f '1' | sort | uniq > $out/tmp/All.list
cat $out/tmp/$out\_to_FlankingRead_80.junction.paf | awk -F'\t' '{print $1"\t"0"\t"$NF"\n"$1"\t"$NF"\t"$2}' > $out/tmp/split.bed

cat $out/tmp/$out\_to_FlankingRead_80.junction.paf | cut -f '1,6' | sort -k2,2 -k1,1 | uniq > $out/tmp/All.list.with_NCL_id

$seqtk_link subseq $longread $out/tmp/All.list > $out/tmp/All.fa
$seqtk_link subseq $out/tmp/All.fa $out/tmp/split.bed > $out/tmp/All_split.fa

echo "Step: to align splited reads against whole genome"

if [ ! -e "$genome_prefix.mmi" ]; then
   $minimap2_link -d $genome_prefix.mmi $genome_prefix.fa
fi

$minimap2_link -t $threads -ax splice $genome_prefix.mmi $out/tmp/All_split.fa | $samtools_link view -bS - > $out/tmp/All.bam

$bedtools_link bamtobed -bed12 -i $out/tmp/All.bam | awk '$1~/^chr[0-9XY]/'| awk '$5>0'| sort -k4,4 > $out/tmp/All.bed12

cat $out/tmp/All.bed12 | awk -F'\t' '{print $4"\t"$0}' | sort -k1,1 > $out/tmp/All.bed12.sorted
cat $out/tmp/All.bed12.sorted | cut -f'1' | uniq -c | awk '$1==1 {print $2}'  > $out/tmp/All.bed12.uniq.list
join -t$'\t' $out/tmp/All.bed12.sorted $out/tmp/All.bed12.uniq.list | cut -f'2-' > $out/tmp/All.uniq.bed12


echo "Step: to diagnose split reads mapped to donor and acceptor: uniquely mapping and same chromsome"

cat $out/tmp/All.uniq.bed12 | awk -F'\t' '{print $4"\t"$0}' | awk 'BEGIN{FS="\t";OFS="\t"}{sub(/:.*/, "", $1); print $0}' | sort -k1,1 -k5,5 > $out/tmp/All.uniq.bed12.with_read_id

touch $(cat $out/tmp/$out\_to_FlankingRead_80.list | awk '{print dir"/"$1".bed12"}' dir=$out/tmp)
join -t$'\t' <(cat $out/tmp/All.list.with_NCL_id | sort -k1,1 -k2,2) $out/tmp/All.uniq.bed12.with_read_id | sort -k2,2 -k1,1 -k6,6 | cut -f '2-' > $out/tmp/All.uniq.bed12.with_NCL_id
$NCLscan_hybrid_bin/split_file_by_first_column.py $out/tmp/All.uniq.bed12.with_NCL_id -o $out/tmp/ -s ".bed12"

$NCLscan_hybrid_bin/check_pass1.py $out/tmp/All.uniq.bed12.with_NCL_id > $out/tmp/pass1.list
mv $(cat $out/tmp/pass1.list | awk '{print dir"/"$1".bed12"}' dir=$out/tmp) $out/pass1/
mv $(join -t$'\t' $out/tmp/$out\_to_FlankingRead_80.list $out/tmp/pass1.list -v 1 | awk '{print dir"/"$1".bed12"}' dir=$out/tmp) $out/fail1/

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
   cat $out/pass1/$one.bed12 | sort -k4,4 > $out/tmp/$one.bed12.tmp
   #one of the end position of split long reads within upstream/downstream 10 bp of donor site
   #the other of the end position of split long reads within upstream/downstream 10 bp of acceptor site
   one_name=$(echo $one | sed 's/\.bed12//g')
   intra_chr=$(echo $one | sed 's/:/\t/g' | awk '{print $1}')
   intra_donor=$(echo $one | sed 's/:/\t/g' | awk '{print $2}')
   intra_acceptor=$(echo $one | sed 's/:/\t/g' | awk '{print $5}')
   upstream_10b_bed=$(echo $intra_donor $intra_acceptor | tr ' ' \\n | sort | head -n 1 | awk '{print chr "\t" $1-10 "\t" $1+10 "\t" "upstream"}' chr=$intra_chr)
   downstream_10b_bed=$(echo $intra_donor $intra_acceptor | tr ' ' \\n | sort | tail -n 1 | awk '{print chr "\t" $1-10 "\t" $1+10 "\t" "downstream"}' chr=$intra_chr)
   echo $upstream_10b_bed"\n"$downstream_10b_bed | tr  ' ' \\t | sort -k1,1 -k2,2n > $out/tmp/$one\_10bp.bed
   
   cat $out/tmp/$one.bed12.tmp | awk '{print $1 "\t" $2-1 "\t" $2 "\t " $4}' > $out/tmp/$one.bed12.tmp.OneSide
   cat $out/tmp/$one.bed12.tmp | awk '{print $1 "\t" $3-1 "\t" $3 "\t " $4}' > $out/tmp/$one.bed12.tmp.TheOtherSide
   cat $out/tmp/$one.bed12.tmp.OneSide $out/tmp/$one.bed12.tmp.TheOtherSide | sort -k1,1 -k2,2n > $out/tmp/$one.bed12.tmp.points
   $bedtools_link intersect -a $out/tmp/$one.bed12.tmp.points -b $out/tmp/$one\_10bp.bed -f 1 -wa -wb | awk '{print $4 "\t" $8}' | sort -k1,1 | uniq > $out/tmp/$one.bed12.tmp.list
   join $out/tmp/$one.bed12.tmp.list $out/tmp/$one.bed12.tmp  -1 1 -2 4 | tr '' \\t  > $out/tmp/$one.bed12.tmp2
   
   cat $out/tmp/$one.bed12.tmp2 | awk '$6>=60' | awk '{print $1"\t"$2 }' | sed -r 's/^(.+):([0-9]+)-([0-9]+)\t([a-z]+)$/\1\t\2\t\3\t\4/g' > $out/pieces.tmp
   cat $out/pieces.tmp | awk -F'\t' '$2==1{print $1":"$3+1"\t"$2"-"$3 "\t" $4}' | sort -k1,1 > $out/pieces1.tmp1
   cat $out/pieces.tmp | awk -F'\t' '$2!=1{print $1":"$2"\t"$2"-"$3 "\t" $4}' | sort -k1,1 > $out/pieces2.tmp1
   #one upstream and the other downstream
   join $out/pieces1.tmp1 $out/pieces2.tmp1 | awk '$3!=$5 {print $1 "\t" $2 "\t" $4}' > $out/pieces12.tmp1

   cat $out/pieces12.tmp1  | sed 's/:/\t/g' | awk '{print $1":"$3}' > $out/pieces1.tmp2
   cat $out/pieces12.tmp1  | sed 's/:/\t/g' | awk '{print $1":"$4}' > $out/pieces2.tmp2
   cat $out/pieces1.tmp2 $out/pieces2.tmp2 | sort > $out/long.tmp1
   join $out/tmp/$one.bed12.tmp $out/long.tmp1 -1 4 -2 1 | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}' > $out/long.tmp2
   
   cat $out/long.tmp2 | grep ':1-' > $out/long.tmp2.1 
   cat $out/long.tmp2 | grep -v ':1-' > $out/long.tmp2.2 
   $bedtools_link bed12tobed6 -i $out/long.tmp2.1 > $out/long.tmp2.1.bed
   $bedtools_link bed12tobed6 -i $out/long.tmp2.2 > $out/long.tmp2.2.bed
   #two split fragements overlapped at least 50bp
   $bedtools_link intersect -a $out/long.tmp2.1.bed -b $out/long.tmp2.2.bed  -s -wo | awk '$13 >=50' | awk -F'[\t:]' '$4==$11{print $0}' | awk '{print $4 "\n" $10}' | sort -k1,1 | uniq > $out/long.tmp2.overlap50.list
   join -t$'\t' <(cat $out/long.tmp2 | sort -k4) $out/long.tmp2.overlap50.list -1 4 -2 1 | sort | uniq | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}' > $out/long.tmp3

   count_long=$(cat $out/long.tmp3 | wc -l)

   if [ "$count_long" -gt 0 ]; then
      
      cat $out/long.tmp3 | awk '$5==60' > $out/pass2_intra/$one.bed12
   fi
    
   #rm -r -f $out/tmp/$one.bed12.tmp
   #rm -r -f $out/tmp/$one.bed12.tmp2

done

echo "Step: to check_intra_OutOfCircle"
## Out of circle: long read extends more than 100 bp on one side (upstream or downstream)  
$NCLscan_hybrid_bin/OutOfCircle.sh -input_folder $out/pass2_intra -o $out -L 100

##intra_view##
$NCLscan_hybrid_bin/BrowserView.sh -input_folder $out/OC_events

echo  "Step: to check_intra_WithinCircle"
join -t$'\t' $out/tmp/$out\_to_FlankingRead_80.paf.with_id.tmp $out/tmp/intra.list | cut -f '2-' > $out/tmp/$out\_to_FlankingRead_80.junction.intra.paf

## Within Circle : obtain more than one pseudo-reference (200 bp NCL sequence) 
$NCLscan_hybrid_bin/WithinCircle.sh -input_folder $out/tmp/$out\_to_FlankingRead_80.junction.intra.paf -c $config -o $out
##Within circle _view##
$NCLscan_hybrid_bin/BrowserView.sh -input_folder $out/WithinCircle_events




## inter ##
cat $out/tmp/pass1_inter.list | while read one
do
   cat $out/pass1/$one.bed12 | sort -k4,4 > $out/pass2_inter/$one.bed12.tmp

   cat $out/pass2_inter/$one.bed12.tmp | awk '$5==60' | cut -f '4' | sed -r 's/^(.+):([0-9]+)-([0-9]+)$/\1\t\2\t\3/g' > $out/pieces.tmp
   cat $out/pieces.tmp | awk -F'\t' '$2==1{print $1":"$3+1"\t"$2"-"$3}' | sort -k1,1 > $out/pieces1.tmp1
   cat $out/pieces.tmp | awk -F'\t' '$2!=1{print $1":"$2"\t"$2"-"$3}' | sort -k1,1 > $out/pieces2.tmp1
   join -t$'\t' $out/pieces1.tmp1 $out/pieces2.tmp1 > $out/pieces12.tmp1

   cat $out/pieces12.tmp1  | sed 's/:/\t/g' | awk '{print $1":"$3}' > $out/pieces1.tmp2
   cat $out/pieces12.tmp1  | sed 's/:/\t/g' | awk '{print $1":"$4}' > $out/pieces2.tmp2
   cat $out/pieces1.tmp2 $out/pieces2.tmp2 | sort > $out/long.tmp1
   join -t$'\t' $out/pass2_inter/$one.bed12.tmp $out/long.tmp1 -1 4 -2 1 | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}' > $out/long.tmp2
   
   cat $out/long.tmp2 | grep ':1-' > $out/long.tmp2.1 
   cat $out/long.tmp2 | grep -v ':1-' > $out/long.tmp2.2 
   $bedtools_link bed12tobed6 -i $out/long.tmp2.1 > $out/long.tmp2.1.bed
   $bedtools_link bed12tobed6 -i $out/long.tmp2.2 > $out/long.tmp2.2.bed
   #two split fragements overlapped at least 50bp
   $bedtools_link intersect -a $out/long.tmp2.1.bed -b $out/long.tmp2.2.bed  -s -wo | awk '$13 >=50' | awk -F'[\t:]' '$4==$11{print $0}' | awk '{print $4 "\n" $10}' | sort > $out/long.tmp2.overlap50.list
   join -t$'\t' <(cat $out/long.tmp2 | sort -k4) $out/long.tmp2.overlap50.list -1 4 -2 1  | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}' > $out/long.tmp3

   count_long=$(cat $out/long.tmp3 | wc -l)

   if [[ "$count_long" -gt 0 ]]; then
  
      cat $out/long.tmp3 | awk '$5>=60' > $out/pass2_inter/$one.bed12
   fi

   rm -r -f $out/pass2_inter/$one.bed12.tmp

done

##inter_view##
#$NCLscan_hybrid_bin/BrowserView.sh -input_folder $out/pass2_inter

rm -rf \
   $out/pieces.tmp \
   $out/pieces1.tmp1 \
   $out/pieces2.tmp1 \
   $out/pieces1.tmp2 \
   $out/pieces2.tmp2 \
   $out/pieces12.tmp1 \
   $out/long.tmp1 \
   $out/long.tmp2 \
   $out/long.tmp2.1 \
   $out/long.tmp2.2 \
   $out/long.tmp2.1.bed \
   $out/long.tmp2.2.bed \
   $out/long.tmp2.overlap50.list \
   $out/long.tmp3


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

