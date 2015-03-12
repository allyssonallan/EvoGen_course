# EvoGen_course
Material for Evolutionary Genomics course - data analysis module

## Requirements

Please note that this tutorial has been designed and tested on unix and mac.

Mac users need the developer tools, as [here](https://developer.apple.com/xcode/) or installed using [homebrew](http://brew.sh/).

If you are using Windows please consider using a virtual machine.
Some quick guidelines can be found [here](http://www.howtogeek.com/170870/5-ways-to-run-linux-software-on-windows).
You can use VMware Player and VirtualBox, and then you can install a light linux distribution like Xubuntu or Lubuntu.
Alternatively, you can use 'cygwin' to recreate a command line environment like in unix/mac.

You should have installed at least SAMtools, FastQC, ANGSD, ngsTools (see details below) as these will be the main programs used.
No knowledge of command line scripting is required.

## Software

You need the following software to run most of these examples:

 - [SAMtools](http://samtools.sourceforge.net/) and [VCFtools](http://vcftools.sourceforge.net/)

 - [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

 - [ANGSD](http://popgen.dk/wiki/index.php/ANGSD)

 - [ngsTools](https://github.com/mfumagalli/ngsTools)

 - [ngsDist](https://github.com/fgvieira/ngsDist)

To run some additional optional exercises you will also need [Picard](http://picard.sourceforge.net/), [Sickle](https://github.com/najoshi/sickle) and [Scythe](https://github.com/vsbuffalo/scythe).

For parsing and plotting you will also need [R](http://www.r-project.org/) and [Perl](http://www.perl.org/) installed. 
You will also need some R packages to run all examples, and these are: 'ggplot2', 'VennDiagram', 'optparse'.

Please read individual documentation for proper installation of each software. As a general reference, we provide below instructions on how to download and install main programs used.


## Data

We will use some example datasets:

\begin{itemize}

\item Raw reads data of chimpmunks are from \href{http://www.ncbi.nlm.nih.gov/pubmed/24118668}{his} paper, and are kindly provided by Tyler Linderoth and Ke Bi.
In this paper, Authors sequenced around 4Mbp from early 20th century samples. We will use a very small subset of the original dataset.

\item BAM files of butterflies are a small sample subset of the original study and we will analyse around 1 Mbp. More information about the original aligned files can be obtained at the original paper \href{http://www.ncbi.nlm.nih.gov/pubmed/22759293}{here}. In this study, authors analysed a butterfly Next-Generation Sequencing dataset sequenced using Illumina GAII technology. The dataset consists of a total of 381 samples of Lycaeides idas, Lycaeides melissa and an hybrid collection (Jackson Hole). DNA resequencing was conducted on custom reduced genomic complexity libraries (RAD-seq). We will analysed only a subset of 20 samples.

\item Other BAM files that will be used are a subset of human sequencing data from the 1000 Genomes Project. More information on this project can be found [[http://www.1000genomes.org/][here]]. This dataset comprises 33 individuals of European descent.

\item Genotype likelihoods from inbred samples will be generated on-the-run.

\end{itemize}





