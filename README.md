# EukRep #
Classification of Eukaryotic and Prokaryotic sequences from metagenomic datasets

Two important notes: 
* EukRep is designed to miss as little eukaryotic sequence as possible and makes predictions on every seqeunce, even if it has low confidence. Because of this, it has a false positive rate around 2% and you should not expect the entire output to be eukaryotic sequence. This also means if your metagenome doesn't have a eukaryotic genome present (which is more common than not in our experience) you will still get output from EukRep.
* EukRep is intended to be used as one part of a larger pipeline. For obtaining high quality gene predictions and binning identified eukaryotic contigs as described in "Genome-reconstruction for eukaryotes from complex natural microbial communities" (West et al. in review), see methods section https://doi.org/10.1101/171355

-or-

See a provided example workflow (work in progess)
https://github.com/patrickwest/EukRep_Pipeline

# Installation #
* Requires Python3
* Easy install with pip
```
$ pip install EukRep
```
# Example Usage #
* Identify and output sequences predicted to be of eukaryotic origin from a fasta file:
```
$ EukRep -i <Sequences in Fasta format> -o <Eukaryote sequence output file>
```
* Identify and output both sequences of eukaryotic and prokaryotic origin from a fasta file:
```
$ EukRep -i <Sequences in Fasta format> -o <Eukaryote sequence output file> --prokarya <Prokaryote sequence output file>
```
