 #!/usr/bin/env zsh

#genome data link
gtf="/media/user/UCLA_ASD_analyses/longread_project_tmp/NCLsl/hg38/gencode.v28.annotation.gtf"
genome_prefix="/media/user/UCLA_ASD_analyses/longread_project_tmp/NCLsl/hg38/GRCh38.p12.genome"

#external tools link
bedtools_link="/media/user/UCLA_ASD_analyses/longread_project_tmp/NCLsl/bedtools2/bin/bedtools"
samtools_link="samtools"
seqtk_link="/media/user/UCLA_ASD_analyses/longread_project_tmp/NCLsl/seqtk/seqtk"
minimap2_link="/media/user/UCLA_ASD_analyses/longread_project_tmp/NCLsl/minimap2/minimap2"



out="ont_to_cirbase_HEK293_Rep1"
CircOne="chr11:33287511:+:chr11:33286413:+"

  
   #echo -n > $out/WithinCircle_tmp/$CircOne.split.bed
   #cat $out/WithinCircle_tmp/CircOne_readID.list | while read ReadOne
   #do
   #   echo $ReadOne
   #   cat $out/WithinCircle_tmp/$CircOne.split.tmp1 | grep -P $ReadOne'\t' | sort -nk6 >  $out/WithinCircle_tmp/ReadOne.split.tmp1
   #   echo $ReadOne | awk '{print SRR "\t" 0}' SRR=$ReadOne > $out/WithinCircle_tmp/ReadOne.split_Start.tmp1
   #   ReadOneEnd=$(cat $out/WithinCircle_tmp/$CircOne.split.tmp1 | grep $ReadOne  | head -n 1 | awk '{print $2}')
   #   cat $out/WithinCircle_tmp/$CircOne.split.tmp1 | grep -P $ReadOne'\t' | awk '{print SRR "\t" $6}' SRR=$ReadOne >> $out/WithinCircle_tmp/ReadOne.split_Start.tmp1
   #  
   #   ReadOneSplitNum=$(cat $out/WithinCircle_tmp/ReadOne.split_Start.tmp1 | wc -l)
   #  
   #   cat $out/WithinCircle_tmp/$CircOne.split.tmp1 | grep -P $ReadOne'\t' | awk '{print $6}' > $out/WithinCircle_tmp/ReadOne.split_End.tmp1
   #   echo $ReadOneEnd >> $out/WithinCircle_tmp/ReadOne.split_End.tmp1
   #   paste $out/WithinCircle_tmp/ReadOne.split_Start.tmp1 $out/WithinCircle_tmp/ReadOne.split_End.tmp1  >> $out/WithinCircle_tmp/$CircOne.split.bed
   #done

   #$seqtk_link subseq $out/tmp/All.fa $out/WithinCircle_tmp/$CircOne.split.bed > $out/WithinCircle_tmp/$CircOne.split.fa
   #echo "Step: to align split sequences of long reads against whole genome"
   #$minimap2_link -ax splice $genome_prefix.fa $out/WithinCircle_tmp/$CircOne.split.fa | $samtools_link view -bS - > $out/WithinCircle_tmp/$CircOne.bam
   #$bedtools_link bamtobed -bed12 -i $out/WithinCircle_tmp/$CircOne.bam | awk '$1~/^chr[0-9XY]/'| awk '$5==60'| sort -k4  > $out/WithinCircle_tmp/$CircOne.tmp.bed12
   
   
   #echo -n > $out/WithinCircle_tmp/WithinCircle_reads_tmp/$CircOne.bed12
   cat $out/WithinCircle_tmp/$CircOne.tmp.bed12 | awk '{print $4}' | sed 's/:/\t/g' | awk '{print $1}' | sort  | uniq > $out/WithinCircle_tmp/CircOne_readID.list2
   cat $out/WithinCircle_tmp/CircOne_readID.list2 | while read ReadOne
   do
      #ReadOne="SRR10612050.1017435"
      echo $ReadOne
      cat $out/WithinCircle_tmp/$CircOne.tmp.bed12 | grep -P $ReadOne':' > $out/WithinCircle_tmp/ReadOne.tmp.bed12
      
      fragmentNum=$(cat $out/WithinCircle_tmp/$CircOne.split.bed | grep -P $ReadOne'\t' | wc -l)
      bedNum=$(cat $out/WithinCircle_tmp/ReadOne.tmp.bed12 | awk '{print $4}' | wc -l)
      chrNum=$(cat $out/WithinCircle_tmp/ReadOne.tmp.bed12 | awk '{print $1}' | uniq | wc -l) 
      chrSplit=$(cat $out/WithinCircle_tmp/ReadOne.tmp.bed12 | awk '{print $1}' | uniq | head -n 1)  
      startSplit=$(cat $out/WithinCircle_tmp/ReadOne.tmp.bed12 | awk '{print $2 "\n" $3}' | sort | head -n 1)
      endSplit=$(cat $out/WithinCircle_tmp/ReadOne.tmp.bed12 | awk '{print $2 "\n" $3}' | sort -r | head -n 1)
    
      echo $fragmentNum
      echo $bedNum
      echo $chrNum
      echo $chrSplit
      echo $startSplit
      echo $endSplit
      
      #echo "Step: to report a circRNA long read"
       if [ "$bedNum" -eq $(($fragmentNum+1)) ] && [ "$chrNum" -eq 1 ] && [ $(($intraMin-10)) -le "$startSplit" ] && [ "$endSplit" -le $(($intraMax+10)) ] && [ "$intraChr"=="$chrSplit" ]
       then 
             cat $out/WithinCircle_tmp/ReadOne.tmp.bed12 >> $out/WithinCircle_tmp/WithinCircle_reads_tmp/$CircOne.bed12
       fi 
   done
