#!/usr/bin/env zsh

#run "ncl_exon_number_extractor.py" to extract exon number of EB and belonging to which trpt ID
cat test_dat/circ_EB_col9_test.txt | /media/treeslab/MCF7_TS_2022/NCLscan-hybrid/bin/ncl_exon_number_extractor.py gencode.v28.annotation.gtf > MCF7_pacbio_MQ0_test/tmp/ncl_exon_number.txt 
#pip intall bedparse (require "bedparse" software)
bedparse gtf2bed gencode.v28.annotation.gtf > MCF7_pacbio_MQ0_test/tmp/annotation.gtf.bed
cat MCF7_pacbio_MQ0_test/tmp/annotation.gtf.bed | grep 'ENST00000262215.7' > MCF7_pacbio_MQ0_test/tmp/OneTrpt.bed12
