
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
You will also need some R packages to run all examples, and these are: 'ggplot2', 'VennDiagram', 'optparse', 'fields'.

Please read individual documentation for proper installation of each software. As a general reference, we provide below instructions on how to download and install main programs used.


## Data

We will use some example data sets:

 - Raw reads data of **chipmunks** are from [this](http://www.ncbi.nlm.nih.gov/pubmed/24118668) paper, and are kindly provided by Tyler Linderoth and Ke Bi.
In this paper, Authors sequenced around 4Mbp from early 20th century samples. We will use a very small subset of the original dataset.

 - BAM files of **butterflies** are a small sample subset of the original study and we will analyse around 1 Mbp.
More information about the original aligned files can be obtained at the original paper [here](http://www.ncbi.nlm.nih.gov/pubmed/22759293).
In this study, authors analysed a butterfly Next-Generation Sequencing dataset sequenced using Illumina GAII technology.
The dataset consists of a total of 381 samples of Lycaeides idas, Lycaeides melissa and an hybrid collection (Jackson Hole).
DNA resequencing was conducted on custom reduced genomic complexity libraries (RAD-seq).
We will analysed only a subset of 20 samples.

 - Other BAM files that will be used are a subset of **human** sequencing data from the 1000 Genomes Project.
More information on this project can be found [here](http://www.1000genomes.org/).
This dataset comprises 33 individuals of European descent.

 - Genotype likelihoods from inbred samples will be generated on-the-run.

All these files can be downloaded [[http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/input.tar.gz][here]]. Untar these files by typing:
<src>
tar -xvf input.tar.gz
rm input.tar.gz
</src>
This file is rather big (~550M)
You can use a smaller version (~100M) if you prefer, but you will not be able to run examples on human dataset.
If this is the case, you can download it [[http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/input_light.tar.gz][here]].
<src>
tar -xvf input_light.tar.gz
rm input_light.tar.gz
mv input_light/ input/ # if you run this you will not be able to use human examples
</src>

I added already generated output files, as a back up:

As a back up, you can find my output files [[http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/myoutput.tar.gz][here]].
<src>
tar -xvf myoutput.tar.gz
rm myoutput.tar.gz
</src>


## Preparation

First, you can create a folder where we will put all material (and results) for this day.
```
mkdir workshop
cd workshop
mkdir output
```

----------

To install the latest version of ANGSD and ngsTools the best way is to first install git, for instance as described [here](http://git-scm.com/downloads).

```
git clone --recursive https://github.com/mfumagalli/ngsTools.git
cd ngsTools
make clean
make
```
You may see a lot of warnings, but as long as they are not errors (and compilation fails), please ignore them.

To view a list of possible options, and check that everything worked for instance for ANGSD, simply type:
```
./angsd/angsd

```

------------

SAMtools can be downloaded [[http://sourceforge.net/projects/samtools/files/][here]]. Please use version 0.1.19. Then type:
<src>
tar xvjf samtools-0.1.19.tar.bz2
rm samtools-0.1.19.tar.bz2
cd samtools-0.1.19
make
make razip # this is optional
cd ..
</src>
To check that everything went fine, type the following commands and check whether help messages are printed:
<src>
samtools-0.1.19/samtools
samtools-0.1.19/bcftools/bcftools
</src>

VCFtools can be obtained [[http://sourceforge.net/projects/vcftools/files/][here]]. Version 0.1.11 is fine. Then type:
<src>
tar xzvf vcftools_0.1.11.tar.gz
rm vcftools_0.1.11.tar.gz
cd vcftools_0.1.11
make
cd..
</src>

-------

For data filtering, we will use:

- [[https://github.com/vsbuffalo/scythe][Scythe]]

<src>
git clone https://github.com/vsbuffalo/scythe.git
cd scythe
make all
cd ..
</src>

- [[https://github.com/najoshi/sickle][Sickle]]

<src>
git clone https://github.com/najoshi/sickle.git
cd sickle
make
cd ..
</src>

- [[http://picard.sourceforge.net/][Picard]]

Download the latest zipped version [[http://sourceforge.net/projects/picard/files/][here]] and unzip it.
You need [[http://www.java.com/en/][Java]] to run Picard tool.

You may want to install FastQC as well, available [[http://www.bioinformatics.babraham.ac.uk/projects/fastqc/][here]].

For the solution of an exercise we will use [[https://code.google.com/p/cutadapt/][cutadapt]].
You can use the git repository to install it from [[https://github.com/marcelm/cutadapt][here]].
You need Python to run it.

------

Lastly, we will also assume that you have perl (for parsing) and R (with Rscript, for plotting) installed in your /usr/bin directory. To test if this is the case, type:
<src>
perl --help
R --help
Rscript --help
</src>
and check that this prints out help messages. Otherwise install them following instructions reported [[http://www.perl.org/get.html][here]] and [[http://www.r-project.org/][here]]. You will also need some R packages to run all examples, and these are: 'ggplot2', 'VennDiagram', 'optparse'. To install a package in R, open R and type:
<src>
install.packages("ggplot2")
install.packages("VennDiagram")
install.packages("optparse")
</src>
You will be asked to select a mirror, and eventually to first install some missing dependencies. Close R by typing 'q()' and type 'n' if you are asked whether to save the current workspace.

Also [[http://en.wikipedia.org/wiki/AWK][awk]] will be sometimes used to parse files.

Files will be un/compressed using gunzip, bzip2 and bunzip2, some details can be found [[http://osxdaily.com/2012/05/29/create-extract-bz2-mac-os-x/][here]] and [[http://linux.about.com/library/cmd/blcmdl1_bunzip2.htm][here]].

A Perl script will use Statistics::Distributions package, available for download [[http://search.cpan.org/~mikek/Statistics-Distributions-1.02/Distributions.pm][here]]. 
The easiest way to install it is to run:
<src>
sudo cpan Statistics::Distributions
</src>
If it fails, download the .tar.gz file, then unzipped it and install it:
<src>
tar -xvzf Statistics-Distributions-1.02.tar.gz
cd Statistics-Distributions-1.02
perl Makefile.PL 
make test
cd..
</src>
You need to make sure that the package has been correctly added your Perl directory.

Likewise you may need to install this additional package [[http://search.cpan.org/~pmqs/IO-Compress-2.064/lib/IO/Compress/Bzip2.pm][IO-Compress]]:
<src>
sudo cpan IO::Compress::Bzip2
# if it does not work try: sudo cpan force install IO::Compress::Bzip2
</src>
or manually:
<src>
tar -xvzf IO-Compress-2.064.tar.gz
cd IO-Compress-2.064
perl Makefile.PL
cd ..
</src>
and Getopt which should be already installed, otherwise see [[http://search.cpan.org/~jhi/perl-5.8.1/lib/Getopt/Std.pm][here]].
Before being worried that something failed, read below for a quick test to check that this perl script indeed works.

------------

You will also need some scripts available [[http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/scripts.tar.gz][here]]. Untar these files by typing:

<src>
tar -xvf scripts.tar.gz
rm scripts.tar.gz
</src>

To test that you have all required packages in Perl type:
<src>
perl ../scripts/SNPcleaner.pl --help
</src>
and see if it prints out something or errors.

Optionally, you can download [[http://www.broadinstitute.org/gatk/][GATK]] and [[https://github.com/ekg/freebayes][FreeBayes]], as they will be briefly mentioned and discussed as supplementary information.

*** Misc

In case no internet connection will be available, you can download some files we will generate [[http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/data.tar.gz][here]].
I will also include PDFs of some of the cited references.
<src>
tar -xvf data.tar.gz
rm data.tar.gz
</src>

As a back up, you can find my output files [[http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/myoutput.tar.gz][here]].
<src>
tar -xvf myoutput.tar.gz
rm myoutput.tar.gz
</src>
















