## Manual of NCLscan-hybrid
####   Version: 1.1.0
#### NCLscan-hybrid, a tool using long-read sequencing (Pabio/Nanopore) to validate non-col-linear (NCL) transcripts (fusion, trans-splicing, and circular RNA) 
-----------------------------


### Requirements

- Python
- bedtools==v2.25.0
- samtools
- minimap2
- seqtk

We recommand to use conda to install the dependencies.


### Installation

```
git clone https://github.com/TreesLab/NCLscan-hybrid.git
```


### Usage

```
./NCLscan-hybrid.sh \
    -long [input long read fasta/fastq file] \
    -long_type [pb or ont] \
    -nclscan [NCLscan result file] \
    -c [configure file] \
    -o [out_prefix_name] \
    -t [number of threads]
```

#### Parameters

Parameter | Description
:-------- | :-----------------
 -long FILE | Long reads dataset.(FASTA or FASTQ)
 -long_type TYPE | The type of the long reads dataset. ('pb' or 'ont')
 -nclscan FILE | The results file from NCLscan.
 -c CONFIG_FILE | Config file.
 -o PREFIX | Prefix for output files.
 -t INT | Number of threads.


#### The format of NCLscan results


 \#  | Column 
 :-: | :----- 
  1  | chr (donor) 
  2  | pos (donor)  
  3  | strand (donor) 
  4  | chr (acceptor) 
  5  | pos (acceptor) 
  6  | strand (acceptor) 
  7  | gene_symbol (donor) 
  8  | gene_symbol (acceptor) 
  9  | is_intragenic


The remaining columns generated by NCLscan are optional for NCLscan-hybrid.



### Outputs

- PREFIX.long_intra.result
- PREFIX.long_inter.result 


##### PREFIX.long_intra.result

 \#  | Column | Description 
 :-: | :----- | :---------- 
  1  | NCL_event_id |  
  2  | #supporting_reads | 
  3  | has_reads_out_of_circle | 
  4  | #reads_out_of_circle | 
  5  | has_reads_rolling_circle | 
  6  | #reads_rolling_circle | 
 7 ~ N | | The remaining columns are from the original input file. 


##### PREFIX.long_inter.result

 \#  | Column | Description 
 :-: | :----- | :---------- 
  1  | NCL_event_id |  
  2  | #supporting_reads | 
 3 ~ N | | The remaining columns are from the original input file. 




### Visualization

To visualize the alignments of supporting reads of an supported NCL event, upload the BED files in the following directories to the UCSC genome browser.

- pass2_intra_BrowserView/
- WithinCircle_events_BrowserView/
- pass2_inter_BrowserView/


