#!/usr/bin/env zsh
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
     -nclscan) 
     NCLresult=$(readlink -f $2)
     shift
     ;; 
     -fl)
     FlankLen=$2
     shift
     ;;
     -g|--genome)
     genome=$(readlink -f $2)
     shift
     ;;
     -gtf)
     gtf=$(readlink -f $2)
     shift
     ;;
     -bedtools)
     bedtoolslink=$(readlink -f $2)
     shift
     ;;
     -mps)
     mpslink=$(readlink -f $2)
     shift
     ;;
     -o) 
     output=$2
     shift
     ;; 
     *)

esac
shift
done

if [[ -z "$NCLresult" ]]; then
   echo ""
   echo "Usage:"
   echo "./FlankingSeq.sh -nclscan [input NCL file] -fl [flanking length] -gtf [annotation.gtf] -g [genome.fa] -bedtools [bedtools v2.25.0 link] -mps [merge_paired_sequences.py link] -o [output_prefix]"
   echo ""
   exit
fi

if [[ -z "$FlankLen" ]]; then
   echo ""
   echo "Usage:"
   echo "./FlankingSeq.sh -nclscan [input NCL file] -fl [flanking length] -gtf [annotation.gtf] -g [genome.fa] -bedtools [bedtools v2.25.0 link] -mps [merge_paired_sequences.py link] -o [output_prefix]"
   echo ""
   exit
fi

if [[ -z "$output" ]]; then
   echo ""
   echo "Usage:"
   echo "./FlankingSeq.sh -nclscan [input NCL file] -fl [flanking length] -gtf [annotation.gtf] -g [genome.fa] -bedtools [bedtools v2.25.0 link] -mps [merge_paired_sequences.py link] -o [output_prefix]"
   echo ""
   exit
fi

echo "flanking length: "$FlankLen
flanklen=$(echo $FlankLen | awk '{print $1}')


## exonic regions ##
cat $gtf | awk '$3=="exon"' | awk '{print $1 "\t" $4-1 "\t" $5 "\t" $3 "\t" 0 "\t" $7}' | sort -k1,1 -k2,2n | $bedtoolslink merge -s | awk '{print $1 "\t" $2 "\t" $3 "\t" "exonic" "\t" 0 "\t" $4}'> exonic_ranges.bed
cat $NCLresult > NCLinput


cat NCLinput | awk '$3=="+" {print $1 "\t" $2-FL "\t" $2 "\t" $1":"$2":"$3":"$4":"$5":"$6".1" "\t" 0 "\t" $3}' FL=$flanklen | sort -k1,1 -k2,2n > NCL_Dplus.bed
cat NCLinput |awk '$3=="-" {print $1 "\t" $2-1 "\t" $2+FL-1 "\t" $1":"$2":"$3":"$4":"$5":"$6".1" "\t" 0 "\t" $3}' FL=$flanklen | sort -k1,1 -k2,2n > NCL_Dminus.bed
cat NCLinput | awk '$6=="+" {print $4 "\t" $5-1 "\t" $5+FL-1 "\t" $1":"$2":"$3":"$4":"$5":"$6".2" "\t" 0 "\t" $6}' FL=$flanklen | sort -k1,1 -k2,2n > NCL_Aplus.bed
cat NCLinput |awk '$6=="-" {print $4 "\t" $5-FL "\t" $5 "\t" $1":"$2":"$3":"$4":"$5":"$6".2" "\t" 0 "\t" $6}' FL=$flanklen | sort -k1,1 -k2,2n > NCL_Aminus.bed

cat NCL_Dplus.bed NCL_Dminus.bed | sort -k1,1 -k2,2n > NCL_flank.1.bed
cat NCL_Aplus.bed NCL_Aminus.bed | sort -k1,1 -k2,2n > NCL_flank.2.bed
$bedtoolslink intersect -a exonic_ranges.bed -b NCL_flank.1.bed -s -wb | awk '{print $1 "\t" $2 "\t" $3 "\t" $10 "\t" $5 "\t" $6}' | sort -k1,1 -k2,2n | awk -F"\t" '!_[$4]++' > NCL_flank.1_exonic.bed 
$bedtoolslink intersect -a exonic_ranges.bed -b NCL_flank.2.bed -s -wb | awk '{print $1 "\t" $2 "\t" $3 "\t" $10 "\t" $5 "\t" $6}' | sort -k1,1 -k2,2n | awk -F"\t" '!_[$4]++' > NCL_flank.2_exonic.bed 

cat NCL_flank.1_exonic.bed NCL_flank.2_exonic.bed | sort -k1,1 -k2,2n > NCL_flank_exonic.bed

$bedtoolslink getfasta -fi $genome -bed NCL_flank_exonic.bed -s -name -fo NCL_flank_exnic.fa
$mpslink NCL_flank_exnic.fa  $output\_flanking_merged.fa

cat NCL_flank.1_exonic.bed | awk '{print $4 "\t" $3-$2}' | sed 's/\.1//g' | tr ' ' \\t | sort -k1,1 > NCL_flank.1_exonic.length
cat NCL_flank.2_exonic.bed | awk '{print $4 "\t" $3-$2}' | sed 's/\.2//g' | tr ' ' \\t | sort -k1,1 > NCL_flank.2_exonic.length
join NCL_flank.1_exonic.length NCL_flank.2_exonic.length > $output\_NCL_flanking.length


echo "Output:" $output"_flanking_merged.fa;" $output"_NCL_flanking.length"



rm -r -f exonic_ranges.bed
rm -r -f NCLinput
rm -r -f NCL_Dplus.bed
rm -r -f NCL_Dminus.bed
rm -r -f NCL_Aplus.bed
rm -r -f NCL_Aminus.bed
rm -r -f NCL_flank.1.bed
rm -r -f NCL_flank.2.bed
rm -r -f NCL_flank.1_exonic.bed
rm -r -f NCL_flank.2_exonic.bed
rm -r -f NCL_flank_exonic.bed
rm -r -f NCL_flank_exnic.fa


