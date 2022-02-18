#!/usr/bin/env zsh

while [[ $# -gt 1 ]]
do 
key="$1"

case $key in 
     -long | --longread)
     longread=$(readlink -f $2)
     shift
     ;;
     -long_type | --longread_type)
     long_type=$2
     shift
     ;;
     -nclscan | --NCLscan)
     NCLscan=$(readlink -f $2)
     shift
     ;;
     -c | --config)
     config=$(readlink -f $2)
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
   echo "./NCLsl.sh -long [input long read fasta file] -long_type [pb or ont] -nclscan [NCLscan result file] -c [configure file link] -o [output_prefix]"
   echo ""
   exit
fi


source $config
BASEDIR=$(pwd)


rm -r -f $out
rm -r -f $out/tmp
rm -r -f $out/pass1
rm -r -f $out/fail1
rm -r -f $out/pass2_intra
rm -r -f $out/pass2_inter

mkdir $out
mkdir $out/tmp
mkdir $out/pass1
mkdir $out/fail1
mkdir $out/pass2_intra
mkdir $out/pass2_inter

echo -n > $out/$out\_long_intra.result
echo -n > $out/$out\_long_inter.result

echo "Step: to create flankingSeqs of NCL events" 

$NCLsl_bin/FlankingSeq.sh -nclscan $NCLscan -fl 100 -gtf $gtf -g $genome_prefix.fa -bedtools $bedtools_link  -mps $NCLsl_bin/merge_paired_sequences.py -o $out/$out\_100bp

echo "Step: to align long reads against flankingSeqs"
$minimap2_link -d $out/$out\_100bp_flanking_merged.mmi $out/$out\_100bp_flanking_merged.fa   
$minimap2_link -t 10 -x map-$long_type $out/$out\_100bp_flanking_merged.fa $longread > $out/tmp/$out\_to_FlankingRead.paf

echo -n > $out/tmp/$out\_to_FlankingRead_overhang10.paf
cat $out/$out\_100bp_NCL_flanking.length  | while read one
do 
   oneEvent=$(echo $one | awk '{print $1}')
   OneSideLength=$(echo $one | awk '{print $2}')
   cat $out/tmp/$out\_to_FlankingRead.paf | grep $oneEvent | awk '($9 -len > 10) && (len - $8 >10)' len=$OneSideLength >> $out/tmp/$out\_to_FlankingRead_overhang10.paf
done 

cat $out/tmp/$out\_to_FlankingRead_overhang10.paf | awk '$11/$7 > 0.8 && $11/$10 < 2 && $12 == 60 {print $6}'| sort | uniq  > $out/tmp/$out\_to_FlankingRead_80.list
cat $out/tmp/$out\_to_FlankingRead_overhang10.paf | awk '$11/$7 > 0.8 && $11/$10 < 2 && $12 == 60 {print $0}'| sort | uniq  > $out/tmp/$out\_to_FlankingRead_80.paf


echo -n > $out/tmp/All.list
echo -n > $out/tmp/split.bed
cat $out/tmp/$out\_to_FlankingRead_80.list | while read one
do 
   cat $out/tmp/$out\_to_FlankingRead_80.paf | grep $one | awk '{print $1}' | sort | uniq > $out/tmp/$one.list 
   cat $out/tmp/$out\_to_FlankingRead_80.paf | grep $one | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' | awk '{print $0 "\t" int(($4-$3+1)/2)}' | awk '{print $0 "\t" $3+$5}' | awk '{print $1 "\t" 0 "\t" $6 "\n" $1 "\t" $6 "\t" $2}' > $out/tmp/$one.split.bed
   cat $out/tmp/$one.list >> $out/tmp/All.list
   cat $out/tmp/$one.split.bed >> $out/tmp/split.bed
done

$seqtk_link subseq $longread $out/tmp/All.list > $out/tmp/All.fa
$seqtk_link subseq $out/tmp/All.fa $out/tmp/split.bed > $out/tmp/All_split.fa

echo "Step: to align splited reads against whole genome"
#$minimap2_link -d $genome_prefix.mmi $genome_prefix.fa
$minimap2_link -t 10 -ax splice $genome_prefix.fa $out/tmp/All_split.fa | $samtools_link view -bS - > $out/tmp/All.bam
$bedtools_link bamtobed -bed12 -i $out/tmp/All.bam | awk '$1~/^chr[0-9XY]/'| awk '$5>0'| sort -k4 > $out/tmp/All.bed12
cat $out/tmp/All.bed12 | awk '{print $4}' | uniq -c | awk '$1==1 {print $2}'  > $out/tmp/All.bed12.uniq.list
cat $out/tmp/All.bed12 | grep -f $out/tmp/All.bed12.uniq.list > $out/tmp/All.uniq.bed12 

echo "Step: to diagnose split reads mapped to donor and acceptor: uniquely mapping and same chromsome"
cat $out/tmp/$out\_to_FlankingRead_80.list | while read one
do 
   cat $out/tmp/All.uniq.bed12 | grep -f $out/tmp/$one.list > $out/tmp/$one\.bed12
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
   cat $out/pass2_intra/$one.bed12.tmp | awk '$5>=60' | awk '$1==chr' chr=$chr | awk '{print $4}' | sed 's/:/\t/g' | grep '1-' | sort | uniq | awk '{split($2,a,"-"); print $1":"a[2]+1 "\t" $2}' | sort -k1 > piece1.tmp1
   cat $out/pass2_intra/$one.bed12.tmp | awk '$5>=60' | awk '$1==chr' chr=$chr | awk '{print $4}' | sed 's/:/\t/g' | grep -v '1-' | sort | uniq | awk '{split($2,a,"-"); print $1":"a[1] "\t" $2}' | sort -k1 > piece2.tmp1
   join piece1.tmp1 piece2.tmp1 | tr ' ' \\t  > piece12.tmp1
   cat piece12.tmp1  | sed 's/:/\t/g' | awk '{print $1":"$3}' > piece1.tmp2 
   cat piece12.tmp1  | sed 's/:/\t/g' | awk '{print $1":"$4}' > piece2.tmp2
   cat piece1.tmp2 piece2.tmp2 | sort -k1 > long.tmp1
   join $out/pass2_intra/$one.bed12.tmp long.tmp1 -1 4 -2 1  | tr ' ' \\t | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}' > long.tmp2
   
   count_long=$(cat long.tmp2 | wc -l)

   if [[ "$count_long" -gt 0 ]]; then
  
      cat long.tmp2 | awk '$5==60' > $out/pass2_intra/$one.bed12
   fi

   rm -r -f $out/pass2_intra/$one.bed12.tmp
done

echo "Step: to check_intra_OutOfCircle"
## Out of circle: long read extends more than 100 bp on one side (upstream or downstream)  
$NCLsl_bin/OutOfCircle_100bp.sh -input_folder $out/pass2_intra -o $out/$out

echo  "Step: to check_intra_WithinCircle"
## Within Circle : obtain more than one pseudo-reference (200 bp NCL sequence) 
$NCLsl_bin/WithinCircle.sh -input_folder $out/tmp/$out\_to_FlankingRead_80.paf -c $config -o $out
##Within circle _view##
$NCLsl_bin/BrowserView.sh -input_folder $out/WithinCircle_events




##intra_view##
$NCLsl_bin/BrowserView.sh -input_folder $out/pass2_intra


## inter ##
cat $out/tmp/pass1_inter.list | while read one
do
   cp $out/pass1/$one.bed12 $out/pass2_inter/$one.bed12.tmp
   cat $out/pass2_inter/$one.bed12.tmp | awk '$5==60' | awk '{print $4}' | sed 's/:/\t/g' | grep '1-' | sort | uniq | awk '{split($2,a,"-"); print $1":"a[2]+1 "\t" $2}' | sort -k1 > piece1.tmp1
   cat $out/pass2_inter/$one.bed12.tmp | awk '$5==60' | awk '{print $4}' | sed 's/:/\t/g' | grep -v '1-' | sort | uniq | awk '{split($2,a,"-"); print $1":"a[1] "\t" $2}' | sort -k1 > piece2.tmp1
   join piece1.tmp1 piece2.tmp1 | tr ' ' \\t  > piece12.tmp1
   cat piece12.tmp1  | sed 's/:/\t/g' | awk '{print $1":"$3}' > piece1.tmp2 
   cat piece12.tmp1  | sed 's/:/\t/g' | awk '{print $1":"$4}' > piece2.tmp2
   cat piece1.tmp2 piece2.tmp2 | sort -k1 > long.tmp1
   join $out/pass2_inter/$one.bed12.tmp long.tmp1 -1 4 -2 1  | tr ' ' \\t | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}' > long.tmp2
   count_long=$(cat long.tmp2 | wc -l)

   if [[ "$count_long" -gt 0 ]]; then
   
      cat long.tmp2 | awk '$5>=60' > $out/pass2_inter/$one.bed12
   fi

   rm -r -f $out/pass2_inter/$one.bed12.tmp

done

##inter_view##
$NCLsl_bin/BrowserView.sh -input_folder $out/pass2_inter

rm -r -f piece1.tmp1
rm -r -f piece2.tmp1
rm -r -f piece1.tmp2
rm -r -f piece2.tmp2
rm -r -f piece12.tmp1
rm -r -f long.tmp1
rm -r -f long.tmp2

echo "Step: output intra result"
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
ls $out/pass2_inter > $out/tmp/pass2_inter.list  
cat $out/tmp/pass2_inter.list | while read one
do  
   cat $out/pass2_inter/$one | awk '{print $4}' | sed 's/:/\t/g' | awk 'NF{NF-=1};1' | sort | uniq | wc -l | awk '{print event "\t" $1}' event=$one | sed 's/\.bed12//g' >> $out/tmp/NCL_long_inter.tmp1
done 
cat $out/tmp/NCL_long_inter.tmp1 | sort -k1 > $out/tmp/NCL_long_inter.tmp2
join $out/tmp/NCL_long_inter.tmp2 $out/tmp/result.tmp | tr ' ' \\t | sort | uniq  > $out/$out\_long_inter.result

