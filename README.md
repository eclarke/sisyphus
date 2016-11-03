# sisyphus: enSemble In-silico aSsembly of PHage Und Similar

Yes, this is the worst acronym I could come up with.

Sisyphus is an ensemble assembly method adapted from Deng et al, _Nucleic Acids Research_, 2014. 
It drops Metavelvet in favor of IDBA_UD and a few other tweaks. 
Importantly, the method uses Snakemake for cluster deployability and resumable operation.
There are also a few additions for blasting and annotating the resulting contigs.

## Usage
Install the dependencies (bioconda is probably your best bet):

- vsearch
- pyfaidx
- pear
- idba_ud
- soapdenovo
- abyss
- cap3
- rust-bio-tools

If you're using the binaries distributed with Deng's EnsembleAssembler, you can skip soapdenovo, cap3 and abyss.

Download the included Snakefile and set up the configuration, notably:
- `fastqs`: should point to the folder containing all your reads (paired). Use absolute paths.
- `filename_fmt`: how your samples are named. Use brackets for the text region with the sample name {sample} and read pair indicator {rp}. For example, MySample5_28_R1.fastq would be `{sample}_{junk}_{rp}.fastq`.
- `bin`: the folder containing the binaries from EnsembleAssembler.


## Method

Sisyphus assembles reads separately with ABySS, SOAPdenovo, and IDBA-UD. 
It also partitions the reads and assembles each partition using ABySS, then collects the resulting contigs.
All the individual resulting contigs from the previous four parts are assembled using the OLC assembler CAP3.


## Anecdotal results

This generally returns longer, more correct contigs than other methods, though sensitivity may be questionable. 
From a spike of Pseudomonas phage phi6, we were able to recover most of the tripartite genome in three contigs.
