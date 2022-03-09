## Manual of NCLscan-hybrid
####   Version: 1.0.0
#### NCLscan-hybird, a tool using long-read sequencing (Pabio/Nanopore) to validate non-col-linear (NCL) transcripts (fusion, trans-splicing, and circular RNA) 
-----------------------------


### Requirements

- Python
- bedtools==v2.25.0
- samtools
- minimap2
- seqtk

recommand to use conda to install the dependencies.


### Usage

```bash
NCLscan-hybrid.sh \
    -long [input long read fasta/fastq file] \
    -long_type [pb or ont] \
    -nclscan [NCLscan result file] \
    -c [configure file] \
    -o [out_prefix_name]
```


### Input files

- long-reads dataset (fasta / fastq)
- the result file from NCLscan
- the config file

### Outputs

- long_intra.result
- long_inter.result 

#### Visualization

- pass2_intra_BrowserView/
- WithinCircle_events_BrowserView/
- pass2_inter_BrowserView/
