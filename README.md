# TP53_SGE



## Instructions for use


### Generating counts from fastqs

The generation of mutant sequence counts from fastqs obtained via targeted sequencing
involved several steps.


#### Strategy to obtain sequence counts from multiplexed paired-end libraries

1. Demultiplexing the libraries

   Sequenced libraries typically contained multiplexed samples or replicates.
   Based on a set of sample-specific barcodes that comprised the fist n nucleotides
   of the R1 and R2 read pairs, the initial step involved demultiplexing the libraries
   and trimming the adapter barcodes.
   
2. Merging paired-end reads

   After demultiplexing, paired-end reads were merged into a single sequence.
   Since the targeted sequencing resulted in overlapping paired-end reads, we
   used this fact to eliminate read pairs that showed divergence in the overlapping
   region. As a result, these pairs could not be merged and were discarded.

3. Counting the sequences

   As our mutant library contained every possible SNV, the merged sequences would potantially
   contain sequences that would differ by only a single nucleotide. Since this precluded
   error-tolerant alignment strategies, counting sequences is reduced to a simple
   exact-matching problem. Moreover, since our mutant library sequences all contained constant flanking regions,
   we trimmed the merged reads so that only the mutatable sequence remained and then simply counted the occurence
   of each synthetic sequence directly.


#### The run script

   The [main.py](main.py) python script contains the complete run script used to create these counts.
   
   It uses pypipegraph2 as a job-scheduling framework in combination with a python environment 
   on a NixOS based server infrastructure.


### Obtaining Relative Fitness Scores (RFS)

The different lengths of the targeted exons resulted in a bias that precluded a direct
comparison of the raw counts or the relative frequencies of mutant sequences between the
exons. In order to obtain a normalized and comparable score across all exons, performed the following normalization:

To calculate the relative frequencies (variant abundances), the read count was divided by the total 
number of matched reads. From this ratio, we obtained Enrichment scores ($ES$) as the negative log2 fold change 
of the variant abundance in treated versus control conditions. 

However, this $ES$ is dependent on the relative amount of wild-type-like and loss-of-function variants in a cell population, which varies between different libraries. To obtain a score that is comparable across different libraries and
screens, the $ES$ was further nromalized into a relative fitness score ($RFS$) by the follwing formula:

$RFS_{ex}(ES) = (\frac{ES - \tilde x_{ex}^{non}}{\tilde x_{ex}^{non} - \tilde x_{ex}^{syn}}) * 2 + 1$

with $\tilde x_{ex}^{non}$ denoting the median of the scores
for all nonsense mutations in a specific exon $ex$ and $\tilde x_{ex}^{syn}$ denoting the median of all synonymous mutations
in this exon.

$RFS$ scores were calculated for each replicate, then, as our total score, we obtained the median ($RFS_{median}$) over
all replicates. 


The calculation can be found in the following jupyter notebook:

[RFS.ipynb](RFS.ipynb)
   
The same transformation was performed on the Enrich2 scores, the corresponding calculation can be found in:

[RFS_enrich.ipynb](RFS_enrich.ipynb)

### Performing Enrich2 analyses

The Enrich2 package was used to perform an additional analysis on our variant counts.
More precisely, we used Enrich2 in count-mode by supplying our variant counts as 
"Identifiers Only" SeqLib. These can be found in the folder "enrich2_files".

The analysis was run via the script [run_enrich.py](run_enrich.py), calling the command line option enrich2_cmd in a
dedicated conda environment as instructed.

For our analysis, we used the configuration file [Exon5678_identifier](enrich2_files/Exon5678_identifier), which is also found in the "enrich2_files" folder.

We used Enrich2 log ratios as scoring method, with DMSO treated samples as T0 and N3a treated samples as T1.
As normalization method we used the library size (full) option.


### Calculating p-values

We used the RFS-transformed scores and standard errors of the enrich2 analysis as calculated above (in [RFS_enrich.ipynb](RFS_enrich.ipynb)) to perform one-sided z-tests for each non-synonymous variant under the null hypothesis that the variant’s score is equal or lower from the weighted mean of the synonymous (WT-like) variants. The alternative hypothesis being that the scores obtained from non-synonymous variants are significantly higher. The z-tests are then adjusted for multiple hypothesis testing using Benjamini-Hochberg correction.

The calculation can be found in this jupyter notebook: [z_statistics_enrich.ipynb](z_statistics_enrich.ipynb)

### How to use it 

The listed python scripts should be run in a python environment with the 
installed dependencies listed below.
Adjust the variables in the section marked as "initialize variables"
to the specifics of your folder structure, if needed.

Jupyter notebooks can be run as-is inside any python environment with pandas 
(version >= 2.0) and jupyter installed.


### Demo

the file [demo.py](demo.py) contains a demo on how to run the counting analysis.

It can be run in a python environment with the appropriate dependencies installed by calling

```bash
python demo.py
```

on the command line.

The neccessary input files are deposited in the demo folder:

[NovaSeq_incoming_sequences_Exon5demo.xlsx](demo/NovaSeq_incoming_sequences_Exon5demo.xlsx) contains
an excel sheet with example barcodes akin to those used in our analysis for the demultiplexing step.

Another excel sheet in the same file contains a couple of predefined sequences that should be found in
the incoming fastq files.

Expected output is a folder with counted sequences and a summary table, that
will list each of the Exon 5 input sequences with a count of 1.

Execution time should be 10 seconds on an average computer.

## System requirements

### Python

All source code was run on python, version >= 3.8.


### Dependencies

The analysis depends on the following python packages that were installed on our python environment:

- pandas (version = 2.1)
- numpy (version = 1.26)
- mbf: https://github.com/imTMarburg/mbf (version = 0.1)
- pypipegraph2: https://github.com/TyberiusPrime/pypipegraph2/ (version = 3.0)
- mmdemultiplex: https://github.com/MarcoMernberger/mmdemultiplex (version = 0.5.0)
- counting_sequences: https://github.com/MarcoMernberger/counting_sequences (version = 0.5)
- NGMerge: https://github.com/jsh58/NGmerge (version = 0.1)
- cutadapt: https://github.com/marcelm/cutadapt (version = 4.1)

The version numbers indicate the versions our scripts have been tested on.

In addition, the packages above may have further dependencies.

### Operating System

The analyses were run on NixOS servers with the anysnake2 package manager. 

Jupyter notebooks were created and run on jupyter-server version >= 1.23.



### Installation

Install pandas, numpy, mbf and pypipegraph2 from pypi via:

```bash
pip install <package>
```

The other packages can be ontained via github using:

```bash
git clone <url>
```

For cutadapt and NGmerge follow package instructions on how to build them.

Typical install time on an everage computer should be less than 5 minutes.
