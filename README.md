# SnakeMake Teaching Example

![Static Badge](https://img.shields.io/badge/maintainer-Garth_Kong-green)

In bioinformatics, it is a very common task to analyze published data using a different approach. However, the choice of (or a lack of) a pipeline manager often dictates how reproducible an analysis can be. In this simple (but realistic) scenario, suppose I have a list of NCBI short read archive (SRA) IDs that I want to analyze. Each SRA ID is related to a sample like the following: 

```
sra          sample
SRR8383512   K562_1_H3K4me1
SRR8383513   K562_2_H3K4me1
SRR8383514   K562_1_H3K4me2
SRR8383515   K562_2_H3K4me2
SRR8383516   K562_1_H3K4me3
SRR8383517   K562_2_H3K4me3
SRR8435051   K562_1_IgG
SRR8435052   K562_2_IgG
```

Given the above list, I want to do the following steps:
1. Download CUT&Tag datasets from the NCBI SRA database
2. Rename the files from SRA ID to eligible sample name
3. Take the first 250K reads
4. Align the reads to a reference genome

The current SnakeMake pipeline will streamline the above steps. Please read the following sections on how to invoke the pipeline. Lastly, please see the accompanying powerpoint slide.

# Quick Start

```
# snakemake version 7
snakemake -j 8 --use-conda

# snakemake version 8
snakemake -j 8 --software-deployment-method conda
```