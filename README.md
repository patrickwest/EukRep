# EukRep #
Classification of Eukaryotic and Prokaryotic sequences from metagenomic datasets

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

# Obtaining Eukaryotic Bins #

EukRep is intended to be used as one part of a larger pipeline. For obtaining high quality gene predictions and binning identified eukaryotic contigs as described in "Genome-reconstruction for eukaryotes from complex natural microbial communities" (West et al. in review), see methods section https://doi.org/10.1101/171355

-or-

See a provided example workflow (work in progess)
https://github.com/patrickwest/EukRep_Pipeline

# Adjusting Stringency #

The stringency of identifying eukaryotic contigs can be adjusted with -m. The false positive rate (FPR) and false negative rate (FNR) for the strict, balanced, and lenient modes are shown below. Default is balanced. Prior to version 0.6.5, lenient was the default.

20kb

<img src="https://github.com/patrickwest/EukRep/blob/master/images/20kb_fpr.png" width="400">

5kb

<img src="https://github.com/patrickwest/EukRep/blob/master/images/5kb_fpr.png" width="400">

Data was obtained by running EukRep on 20kb and 5kb fragmented scaffolds from genomes from mock novel phyla.

# Important Caveat #

In our experience, most metagenomes do not have a eukaryotic genome present; however, EukRep has a false positive rate and you will still receive output in these cases. 
