
**WORKFLOW**:
LAB WORK > SEQUENCING > LOW-LEVEL DATA (reads)

Sequencing machines produce millions of DNA sequences.
Some basic quality control checks are essential to avoid that your biological conclusions are not biased by some unusual features of the data.

## Filtering reads

**WORKFLOW**:
LAB WORK > SEQUENCING > LOW-LEVEL DATA (reads) > FILTERING READS

First step is to filter our and/or trim reads in order to remove low quality data.
Cleaning reads is an important step to decrease the likelihood that alignment and sequencing errors are mistaken as SNPs.

There are several ways and tools one can use to achieve this goal.
Here we will use some custom scripts (mostly borrowed from F.G. Vieira and T. Linderoth) and several external programs.
We will also use FastQC, available [here](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

------------

Sequencing machines usually output Quality Control (QC) files which can then be interpreted to analyze and solve which problems might affect your data.
It is convenient to perform a QC analysis before and after your filtering procedure, to check that all problems have been solved.

**WORKFLOW**:
LAB WORK > SEQUENCING > LOW-LEVEL DATA (reads) > FILTERING READS (FastQ)

In a **FastQ** file quality scores are associated at each called base.
These scores are usually in -10log10(e) where e is the error rate.
This is also called a *Phred quality score*.
Therefore, a Q (quality score) of 10 implies an error every 10 bases, so base call accuracy of 90%; Q20 is an accuracy of 99% and Q30 of 99.9%.
Thus, each base pair has a raw quality score associated.
Calibrated quality scores (using known sequences for instance) are generally more accurate.

Let`s have a look at an example:

```
   less input/chipmunk/IND01_R1.fq
   @HWI-ST745_0097:1:1101:1001:1000#0/1
   ATATGAACATGGAATAGGAGGGTGGTTGAAGTCCCAATGGTATGACTTACTGAGGCTCCATCCGGGATCGAGAGCTTCTTACGGACACGTATTTGCGACC
   +
   CCCCCCCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJIJJHHJIJJJJJJJEIJJJJJJHHHHFFFFDEEEDECC>>>ACCAC
```

Lines are:

order | value
----- | -----
first | sequence identified
second | raw sequence
third | additional description
fourth | quality values

More details can be found [here](http://maq.sourceforge.net/fastq.shtml).
Usually, these files have file extension .fq or .fastq.

Quality values are ordered ASCII characters

  `!""#$%&()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_''abcdefghijklmnopqrstuvwxyz{|}~`

Therefore, you can see that the highest quality is in the middle of the read.

-----------------

FastQ files are useful to retrieve some basic statistics on the quality of your data.
There are several reasons why you may have low quality reads, like bad quality of DNA template, problems with reagents, and adapter contamination.

With the program **FastQC**, we can quickly inspect:

Statistic | Description
--------- | -----------
per base sequence quality | overview of the range of quality values across <br>
 | all bases at each position (quality usually degrades during the run, thus towards the end of the read);<br>
per sequence quality scores | check whether a subset of sequences have <br>
 | global low quality scores;<br>
per base sequence content | check whether there is difference in base compositions <br>
 | along the read (e.g. due to contamination or degradation);<br>
per base GC content | (see above);<br>
per sequence GC content | the expectation should be a Normal distribution;<br>
per base N content | check the amount of not-called bases;<br>
sequence length distribution | check lengths of reads;<br>
duplicate sequences | high levels of (exact) duplication may suggest enrichment bias;<br>
over-represented sequences | over-representation may indicate contamination;<br>
over-represented K-mers | check the relative enrichment along the read.
<br>

Open a .fq file in the folder `input/chimpmunk/` with FastQC (for instance `IND02_R1.fq`).
Inspect the statistics, especially *per base sequence quality* and *sequence duplication levels*.
Do we need to filter out some reads?

As an example, [here](http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/IND02_R1.fq_fastqc.zip) you can find the zipped summary report for `IND02_R1.fq`.

Here we can see some of these QC plots: 

- [per base quality](http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/per_base_quality.pdf)

- [per sequence quality](http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/per_sequence_quality.pdf)

-------------------

There are numerous steps to clean your data at this stage.
Here we do not attempt to provide an exhaustive way to achieve this goal, but rather suggesting some simple steps to take.
Now we perform a very basic filtering reads:

Reads filtering | Description | Notes
--------------- | ----------- | -----
removing adapters | remove adapters from your reads | particularly useful for de novo assembly; you need a list of adapter sequences <br>
trimming low quality bases | remove bases below a certain threshold of quality and/or up to a minimum length | done especially towards the end of the reads; you can remove such bases (hard clipping) or ignore them (soft clipping) <br>
removing contamination from human and bacteria | in case you suspect contamination | 
removing low complexity regions | remove long streches of As,Ts,Cs, and Gs | 
removing duplicate reads | remove reads which arised by PCR duplication | due to over-amplification of libraries during PCR <br>
merging overlapping paired-ends reads | in case of low quality degraded ancient DNA, or small insert size | 

We will use and cite several tools, like `scythe`, `sickle`, and `picard`.

------------------------

First, we want to **remove adapters** in our reads (if any).
Inclusion of adapters in raw reads may seriously affect your assembly.
We will use `scythe` for this purpose, and fasta files of adapters sequences.

Since our example data is paired-end, we need to perform this process for both files.
You can specify the minimum match length argument (-n) and the minimum length of sequence (-M) to keep after trimming:
```
./scythe/scythe -a input/chipmunk/TruSeq2-PE.fa -M 20 -n 5 -p 0.1 input/chipmunk/IND02_R1.fq > output/chipmunk.R1.adapt.fq

	#Adapter Trimming Complete
	#contaminated: 309, uncontaminated: 1622, total: 1931
	#contamination rate: 0.160021

./scythe/scythe -a input/chipmunk/TruSeq2-PE.fa -M 20 -n 5 -p 0.1 -m output/scythe_matches.txt input/chipmunk/IND02_R2.fq > output/chipmunk.R2.adapt.fq

	#Adapter Trimming Complete
	#contaminated: 382, uncontaminated: 1549, total: 1931
	#contamination rate: 0.197825

```
We can even have a look at all matches found by the program
```
less output/scythe_matches.txt

 p(c|s): 1.000000; p(!c|s): 0.000000; adapter: P7_index1
 HWI-ST745_0097:2:1101:1009:1000#0/2
 CAAGCAGAAGACGGCATACGAGATcgtgatGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
 ||||||||||||||||||||||||      ||||||||||||||||||||||||||||||||||
 CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
 @@@@@@B?DFFFEDHHHFJGIIJJJIHIJDICEGIJFHHIJGIGIJHIIIIIIJFGJIJCHGII
 [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]

```

As an additional note, in case of large files, it is recommended to gzip your fastq files, e.g. using `pigz` or `gzip`.

----------------------

Second, let us **trim both ends of the reads** up to a certain threshold of quality and length.
We can use `sickle` for this purpose.
Sickle uses sliding windows along with quality and length thresholds to determine when quality is sufficiently low to trim both ends of reads.
A list of options (for paired-ends reads) can be obtained typing:
```
./sickle/sickle pe
...
-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.
-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.
-n, --discard-n, Discard sequences with any Ns in them.
...
```
More details can be found [here](https://github.com/najoshi/sickle). 

```
./sickle/sickle pe -t sanger -f output/chipmunk.R1.adapt.fq -r output/chipmunk.R2.adapt.fq -o output/chipmunk.R1.trim.fq -p output/chipmunk.R2.trim.fq -s output/chipmunk.RU.trim.fq -q 10 -l 35

 FastQ paired records discarded: 8 (4 pairs)
 FastQ single records discarded: 112 (from PE1: 55, from PE2: 57)
```

Here I used a very relaxed threshold on minimum quality (you can use 20 if you want to be more conservative) and roughly one third of the original reads length.
We can check how many more reads we discard with a stricter filtering.

```
./sickle/sickle pe -t sanger -f output/chipmunk.R1.adapt.fq -r output/chipmunk.R2.adapt.fq -o output/chipmunk.R1.trim.fq -p output/chipmunk.R2.trim.fq -s output/chipmunk.RU.trim.fq -q 20 -l 50

 FastQ paired records discarded: 12 (6 pairs)
 FastQ single records discarded: 131 (from PE1: 63, from PE2: 68)
```

Outputs are both single files and a combined file.

Minimum quality should be set to 10 and at least 50 bases should be retain to obtain a reliable mapping.


###EXERCISE
There are many programs to perform these tasks and the idea it to try different tools with different options and choose the one you feel more comfortable with.
For instance [cutadapt](https://code.google.com/p/cutadapt/), [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), and [FastX-toolkit](http://hannonlab.cshl.edu/fastx_toolkit/) are valid alternatives to remove adapters, in a more extensive way, and trim low quality bases.
Perform reads filtering using one these programs on your dataset or on provided fastq files.

Note that trimming options are a tradeoff between increased accuracy and reads removal.
In case of high-depth data it is reasonable to use a stringent threshold on quality.
In case of low quality data, you may want to be more relaxed in order not to remove too much data (e.g. and thus having very short reads more difficult to align).

As a further check of sanity, once you have called variants in your dataset, you can check their distribution along each read.
If your trimming worked, this distribution should be even across each read.

As a reference on various methods for read trimming see [this](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0085024) paper.

Now we can check the new QC report using FastQC.
Open file: `output/chipmunk.R2.trim.fq`.
Can you see any improvement?

You can see some plots here:

* [per base quality](http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/per_base_quality_trim.pdf)

* [per sequence quality](http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/per_sequence_quality_trim.pdf)

Which problems still persist?

* [K-mer distribution](http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/kmer_profiles.pdf)


--------

**WORKFLOW**:
LAB WORK > SEQUENCING > LOW-LEVEL DATA (reads) > FILTERING READS (FastQ , SAM) > REMOVE READS

Other additional checks and analyses can be performed at this stage to further refine our data.
These will depend on the overall quality of your dataset and your specific goals.

For instance, FastQC results show that there is a non negligible fraction of **duplicate reads** we may need to remove such, as these may be PCR duplicates. 
We can use Picard to do this, after you have aligned your reads:
```
# (do not run this command) 
#java -jar ./picard-tools-1.108/picard-tools-1.108/MarkDuplicates.jar INPUT=myfile.sam OUTPUT=myfile.filtered.sam METRICS_FILE=duplicates.txt REMOVE_DUPLICATES=true
```
although you can use your own scripts (and processing the raw reads directly) and other programs as well (e.g. SAMtools).

Another important step you may want to consider is to remove reads from potential **contaminant** sources. 
You can identify and remove reads that align to the human genome (if you are not analyzing human data...) or other contaminats (e.g. bacteria).
Several scripts can be found for this purpose, for instance the ones provided by Sonal Singhal [here](https://sites.google.com/site/mvzseq/original-scripts-and-pipelines).

If you have paired-end data, low quality data and small insert sizes, some of the reads will overlap each other. 
You need to be aware of this otherwise you may have biases in SNP calls and expression counts.
In general, **merging paired-end reads** improves sequence quality. 
You can use [FLASH](http://ccb.jhu.edu/software/FLASH/) for this purpose:
```
# (do not run this command)
# from FastQ files
#./flash file_R1.fq file_R2.fq --min-overlap=NUM --max-overlap=NUM --max-mismatch-density=NUM ...
```
and you can set the minimum and maximum overlap length between two reads, as well as the mismatch ratio in the overlap.
Another tool to merge paired-end reads is [PEAR](http://sco.h-its.org/exelixis/web/software/pear/).

You can use GATK to perform a **local realignment** around indels, in particular when information on known indels is available.
Likewise, when external information is available, you can recalibrate of quality scores (those provided by the machines are known to be not accurate) using known variants (or polymorphisms in non-variable regions like the human Y chromosome).

An additional program is [FreeBayes](https://github.com/ekg/freebayes), a variant caller software, which performs an indel realignment step internally.
It also automatically detect haplotypes, so in theory you do not need a base quality recalibration.  

When you are analyzing ancient samples, **damage** patterns affect DNA sequences generated by high-throughput sequencing machines.
It is essential to look at the damaging pattern to avoid overcalling SNPs, and remove substitutions that are likely to originate from post-mortem degradation.
A program to estimate parameters of ancient DNA damage pattern is [MapDamage](http://ginolhac.github.io/mapDamage/][MapDamage).

There is a vast amount of programs for NGS data processing, most of them collected at [OmicTools](http://omictools.com/][OmicTools).
For your reference, a comprehensive script that performs all above cited filtering steps can be found in the downloaded `scripts/` folder (named `scrubReads.pl`), kindly provided by Ke Bi at UC Berkeley.

*** Reads alignment

**WORKFLOW**:
... FILTERED READS > ASSEMBLY

In this tutorial we will not practise how to perform a de-novo assembly or reads alignment due to time constraints.

Assemblers merge short reads into contigs, then into scaffolds, and finally into chromosomes.
At least 30X of depth is required for a de novo alignment.

There are several tools that can be used to align reads:

* [BWA](http://bio-bwa.sourceforge.net/)
* [Bowties2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
* [NovoAlign](http://www.novocraft.com/main/index.php)
* [SOAP](http://soap.genomics.org.cn/soapaligner.html)
* [Stampy](http://www.well.ox.ac.uk/project-stampy)
* [SNAP](http://snap.cs.berkeley.edu/)
* [Abyss](http://evomics.org/learning/assembly-and-alignment/abyss/)
* ...

BWA, Botwie2 and NovoAlign have been shown to perform the best.
It is important to test with different programs and different options.
If you do not have a good reference, it is recommended to use a slow but accurate aligner like NovoAlign.

Each alignment has a quality score (MAPQ) associated, similar to a Phred score.
These scores can be used for further filtering of low quality mapping of reads.


*** Filtering sites and individuals

**WORKFLOW**:
... FILTERED READS > ASSEMBLY > MAPPED READS

Now we have aligned reads, usually in SAM/BAM format. 
The great majority of read aligners output alignments in SAM format.

**BAM** file is SAM compressed in the BGZF format.
This enables good compression and efficient random access.
The latter can be achieved through BAM **indexing**, and allows to quickly retrieve positions without first scanning the whole alignment.
BAM files must be sorted before being indexed.
As a recall, other two file formats particularly common in genomics are BED and GFF for annotation purposes.
More details can be found [here](http://samtools.sourceforge.net/SAMv1.pdf).

**SAMtools** offers a good set of programs to visualize and handle SAM/BAM files.
Another program to visualize reads alignments is [Tablet](http://bioinf.scri.ac.uk/tablet/index.shtml).
We will later see how SAMtools can be useful to assign SNPs and genotypes.

Let us have a look at all the options in SAMtools to convert from SAM to BAM:
```
./samtools-0.1.19/samtools view
```
Therefore, if we want to convert a file from SAM to BAM we need to run something like:
```
#./samtools view -bS in.sam > out.bam
#./samtools view out.bam > in.sam # the opposite, from BAM to SAM 
```

Let us use the butterly dataset to illustrate these example.
```
# convert BAM to SAM
./samtools-0.1.19/samtools view input/lyca/bam/lan_10_09f.bam > output/lan_10_09f.sam
less -S output/lan_10_09f.sam
```

The other commands we will use are `sort` and `index`, to quickly retrieve positions without the need to scan the entire data first.
In case we use a reference sequence (and visualize the alignment), we need to index its sequence as well.
`sort` is a prerequiste for `index`.

```
./samtools-0.1.19/samtools sort input/lyca/bam/lan_10_09f.bam output/lan_10_09f.sorted
./samtools-0.1.19/samtools index output/lan_10_09f.sorted.bam
./samtools-0.1.19/samtools faidx input/lyca/referenceseq.fasta
```

We can now have a look at the reads alignment, using `tview` utility in SAMtools.
It assumes that indexed files are in the same directory of the original BAM and FASTA files.
```
./samtools-0.1.19/samtools tview output/lan_10_09f.sorted.bam input/lyca/referenceseq.fasta
```
Type `q` or `ctrl-c` to close it.
This can be useful to inspect specific regions like indels or inversions.
A more advanced tool for viewing alignements is [Tablet](http://bioinf.scri.ac.uk/tablet/index.shtml).

-------------

**WORKFLOW**:
... FILTERED READS > ASSEMBLY > MAPPED READS > FILTERING INDIVIDUALS

Now that we have our aligned reads data and had a look at the alignment, we can proceed with data filtering.
Specifically, now we will use some custom scripts and SAMtools to filter out sites and samples with unusual features.

First thing we may want to do is to **remove abnormal samples**, typically with an unusual sequencing depth compared to the average.
Please note that we need to remove samples with extremely both low and high depth compared to the rest.
You should also consider how many samples you have and you are willing to lose at this stage.

For instance, in the original complete dataset on butterflies, we had 381 samples.
Let us check the distribution of per-sample average depth [here](http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/depth_distrib.pdf).

Which samples should we remove?
To get this distribution you can run the pipeline described below a first time (or better a quicker and more relaxed one) to get such distribution.
Then you remove those aberrant samples (simply remove them from the list of BAM files) and rerun everything.

---------------------

**WORKFLOW**:
... > MAPPED READS > FILTERING SITES

We can now **filter out sites** with low quality or unusual depth.
The latter can be achieved by running the pipeline described below once and check the empirical distribution of per-site depth.
There are many steps one can take in this regard.
Here we will highlight some of them.

We will use a combination of custom scripts and SAMtools for this purpose.
Additionally once can simply use available tools like **[VCFtools](http://vcftools.sourceforge.net/docs.html)** and the `vcftools.pl` script included in [BCFtools](http://samtools.github.io/bcftools/). 

**WORKFLOW**:
... > MAPPED READS > FILTERING SITES (mpileup)

As a first step we will use SAMtools to apply mapped reads quality filters.
We will use **mpileup** utilities implemented in SAMttols for this purpose.

First, let us make sure that our BAM files are sorted and indexed:
```
ls input/lyca/bam/*bam > input/lyca/bams.list
for i in input/lyca/bam/*.bam; do samtools-0.1.19/samtools index $i; done
./samtools-0.1.19/samtools faidx input/lyca/referenceseq.fasta
```

Now we are ready to run the mpileup command on a subset of the data:
```
./samtools-0.1.19/samtools mpileup -AI -q 0 -Q 20 -C 50 -f input/lyca/referenceseq.fasta -b input/lyca/bams.list -r ref_contig:1-10000 > output/lyca.mpileup 2> /dev/null
```

Options `-q -Q -C` are related to reads quality filtering:
* -Q sets a cutoff on the Base Alignment Quality (BAQ), a metric used to remove false SNPs caused by nearby INDELs.<br>
* -C reduces the effect of reads with excessive mismatches.<br>
* -q removes reads with low mapping quality.<br>

Let us have a look at the resulting file.
```
less -S output/lyca.mpileup

 ...
 ref_contig      125     T       18      ..................      GHCHHG@FHHGH>HHHH3      16      ................        ?HBHG2GHHHHGHHHE        30      .....
 ref_contig      126     T       19      ...................     GHGIIHC<IIIGI>IIIHG     17      .................       B0I4HH9GFHDIDIIIE       30      .....
 ref_contig      127     T       19      ...................     GH8IIH@>III<I<IIIGE     17      .................       >:IEHH6GHDHIGHIIG       30      .....
 ref_contig      128     C       19      ...................     IHGIIG@;IIHGI>GIIFE     17      .................       /II?EH:IEBBI>HIGD       29      .....
 ...
```
This file is in (multiple)-**pileup** format.
It represents bases information at each position (on rows).

Rows are defined as:

nr | value
-- | -----
1 | chromsome <br>
2 | position <br>
3 | reference base <br>
4 | total sequencing depth <br>
5 | read bases (compared to the reference) <br>
6 | base qualities <br>
7 | mapping qualities <br>
.  | (after the ^ character for the end read) <br>

Columns 4-5-6 are repeated for each individual.

Bases are encoded as matches (`.` and `,` for forward and reverse strand) and mismatches (`ACGTN` or `acgtn` for forward and reverse strand), while `+` and `-` indicate indels.
`$` represents the end of a read, `^` the start of a read.
More details can be found [here](http://samtools.sourceforge.net/pileup.shtml) and [here](http://samtools.sourceforge.net/samtools.shtml).

-------

**WORKFLOW**:
... > MAPPED READS > FILTERING SITES (mpileup , VCF/BCF)

For convenience, we convert this format (which anyway is useful for visual inspection) to a different one, called BCF (binary variant caller format).
We will later describe this format when practising SNP calling with SAMtools.
So far, just keep in mind that it reports genotypes for each individual at each site along with several metrics of quality.
Therefore, it is extremely suitable for implementing some data filtering at this stage.

VCF (the non-binary version of BCF) is becoming the standard format for storing variable sites from sequencing data. It works with SNPs, indels and structural variations.
VCF files have a header and a data section. The header contains a line starting with `#` with the name of each field, plus several lines starting with `##` containing additional information. The data section is TAB delimited and each line has at the least the first 8 of these fields:

column | value
------ | -----
CHROM | Chromosome name
POS | Position (1-based)
ID | Identifier
REF | Reference (base)
ALT | Alternative (base)
QUAL | Probability of all samples being homozygous reference (and thus not variable)
FILTER | List of filter that this variant failed to pass
INFO | List of misc information
FORMAT | Format of the individual genotypes
... | Individual genotypes

--------------

We can rerun mpileup command and pipe it output dircetly to **bcftools**, a program that analyzes BCF files.
We need to specify `-g` option to print out BCF format (and not pileup format).
Moreover, `-D` and `-S` options will output depth and strand bias p-values, metrics that we will use for data filtering. 
As you can see bcftools is primarly used to call SNPs and genotypes.
Let us ignore this now and print out all sites, whether they might be variable or not in our sample. 

```
./samtools-0.1.19/samtools mpileup -u -AIDS -q 0 -Q 20 -C 50 -f input/lyca/referenceseq.fasta -b input/lyca/bams.list -r ref_contig:1-10000 2> /dev/null | ./samtools-0.1.19/bcftools/bcftools view -gI - > output/lyca.raw.vcf 2> /dev/null
```
Again, for large datasets it is recommended to gzip output files (e.g. by piping out `pbzip2 > out.vcf.bz2`). 

You can quickly inspect this file:
```
tail -n 5000 output/lyca.raw.vcf | head -n 5

 ref_contig	3141	.	G	.	45.2	.	DP=178;AF1=0;AC1=0;DP4=165,0,0,0;MQ=47;FQ=-45.1	PL:DP:SP	0:7:0	0:6:0	0:9:0	0:9:0:11:0	0:12:0	0:2:0	0:17:0	0:9:0	0:4:0	0:10:0	0:6:0	0:12:0	0:10:0	0:3:0	0:4:0	0:11:0	0:11:0	0:7:0	0:5:0
 ref_contig	3142	.	G	.	45.2	.	DP=178;AF1=0;AC1=0;DP4=165,0,0,0;MQ=47;FQ=-45.1	PL:DP:SP	0:7:0	0:6:0	0:9:0	0:9:0:11:0	0:12:0	0:2:0	0:17:0	0:9:0	0:4:0	0:10:0	0:6:0	0:12:0	0:10:0	0:3:0	0:4:0	0:11:0	0:11:0	0:7:0	0:5:0
 ref_contig	3143	.	T	.	45	.	DP=178;RPB=1.709737e+00;AF1=0;AC1=0;DP4=154,0,1,0;MQ=47;FQ=-44.9;PV4=1,0.41,0.42,0.47	PL:DP:SP	0:7:0	0:6:0	0:8:0	0:8:0	0:10:0	0:12:0	0:2:0	0:17:0	0:7:0	0:4:0	0:10:0	0:6:0	0:11:0	0:8:0	0:3:0	0:4:0	0:9:0	0:11:0	0:7:0:5:0
 ref_contig	3144	.	A	.	45.2	.	DP=178;AF1=0;AC1=0;DP4=166,0,0,0;MQ=47;FQ=-45.1	PL:DP:SP	0:7:0	0:6:0	0:10:0	0:9:0:11:0	0:12:0	0:2:0	0:17:0	0:9:0	0:4:0	0:10:0	0:6:0	0:12:0	0:10:0	0:3:0	0:4:0	0:11:0	0:11:0	0:7:0	0:5:0
 ref_contig	3145	.	T	.	44.6	.	DP=178;AF1=0;AC1=0;DP4=162,0,0,0;MQ=47;FQ=-44.4	PL:DP:SP	0:7:0	0:6:0	0:10:0	0:8:0:11:0	0:12:0	0:2:0	0:17:0	0:9:0	0:3:0	0:10:0	0:6:0	0:12:0	0:9:0	0:3:0	0:4:0	0:11:0	0:11:0	0:7:0	0:4:0

```

You can see at each positions several metrics are reported, like the total depth (DP), mapping quality (MP), strand bias (SP), and so on.
We can make full use of this file to filter our entries that do not meet our quality requirements.

For instance, here we will use a custom Perl script, although there are many available tools for that.
The scope here is to show the main steps for data filtering of sites, rather than giving a fixed pipeline.
Our aim here is also to illustrate how filtering indeed changes the data and how important it is to keep track of such changes.


-----------------------

**WORKFLOW**:
... > MAPPED READS > FILTERING SITES (mpileup , VCF/BCF)

Main quality control filters can be grouped as:

1. Depth
* Minimum site read depth <br>
* Maximum site read depth <br>
* Even depth across individuals <br>
* Minimum number of alternate alleles <br>

2. Bias and other quality aspects
* Minimum RMS mapping quality for SNPs <br>
* Mapping quality bias <br>
* Strand bias <br>
* Allele bias in potential heterozygotes <br>
* Base quality bias <br>
* Distance from end of read bias <br>

3. Hardy-Weinberg Equilibrium
* Exact test of HWE <br>
* Excess of heterozygotes <br>
* Excess of homozygotes <br>

4. Biallelic site filter

5. Mutation type removal filters


Read depth cutoffs can be chosen by first investigating the empirical distribution.
Again, we may want to remove sites with excessively high or low depth compared to the distribution.

Strand bias measures the imbalance for the number of reads at each strand. The end read bias test looks for unusually high/low number of alternate alleles at the of the reads, which could be due to bad trimming.
Allele bias in potential heterozygotes tests for potential barcode/allele swapping.

Let us see the options for this script (implemented by T. Linderoth and F.G. Vieira and others).
This script reads a VCF file and filter SNPs based on a set of rules. 
It is similar to SAMtools `varfilter.pl` script but it extends some of the filters.
Please note that this script is only for SNP filtering and will ignore INDELs (which we have already filtered out using SAMtools).

```
perl ../scripts/SNPcleaner.pl --help

Usage: 
SNPcleaner.pl [OPTIONS] <infile.vcf>
or
cat <infile.vcf> | SNPcleaner.pl [OPTIONS]

############ OPTIONS ############

input files:

--pop_info|-P	FILE	input file with population information (each line should list `sample_name	pop_ID`) for HWE filtering
--anc|-A	FILE	ancestral-state fasta file (with FAI in same directory)
--pileup|-G	FILE	pileup format file (all sites and individuals in input vcf must be in pileup)
--exons|-X	FILE 	BED format file of exonic regions (sorted from lowest to highest numbered contig)

coverage filters:

--minDepth|-d	INT	minimum site read depth [2]
--maxDepth|-D	INT	maximum site read depth [1000000]
--minIndiv|-k	INT	minimum number of individuals with at least [-u INT]X coverage (requires SNPcleaner -u and mpileup -D) [1] 
--minIndiv_cov|-u	INT	minimum individual coverage threshold used for -k (requires SNPcleaner -k and mpileup -D) [0]
--minalt|-a	INT	minimum number of alternate alleles per site [0]

bias and other quality-aspect filters:

--minRMSmap|-Q	INT	minimum RMS mapping quality for SNPs [10]
--mapqual|-f	FLOAT	min p-value for map quality bias [0]
--strand_ind|-S	FLOAT	min p-value for strand bias from combining p-values across individuals [0.0001]
--strand_site|-s	FLOAT	min p-value for strand bias determined from read counts summed over all individuals at the site [0.0001]
--allele_bias|-T	FLOAT	min p-value for allele bias in potential heterozygotes [0.001]
--hetero_llh|-R	FLOAT	skip called homozygotes with heterozygote likelihood less than -R FLOAT for allele bias filter [0.05] 	
--basequal|-b	FLOAT	min p-value for base quality bias [1e-100]
--endbias|-e	FLOAT	min p-value for biased distance of alternate bases from ends of reads (indication of misalignment) [0.0001]

Hardy-Weinberg equilibrium filters:

--hwe|-h		FLOAT	min p-value for exact test of HWE (two-tailed) [0]
--hetero_excess|-H	FLOAT	min p-value for exact test of excess heterozygotes [0]
--hetero_deficit|-L	FLOAT	min p-value for exact test of deficient heterozygotes [0]
--inbreed_coef|-F	FLOAT	inbreeding coefficient value [0]
--rmv_nonHWEexons|-g	filter-out exons containing at least one SNP out of HWE (requires -X)

mutation type filters:

--mutation_rmv|-M	STR	mutation type(s) to remove (ex: `-M CT_GA` means remove C<=>T and G<=>A)
--one_dir|-w	remove mutations defined by -M in one direction (ex: `-M CT_GA` means remove C=>T and G=>A) (requires -A)
--alt_excess|-J	FLOAT	min p-value for excess substitutions defined by -M in sites called as nonvariable (requires -A if -w) [1e-06] 
--error|-E	FLOAT	sequencing error rate (required by -J for filtering mutation types) [0.01]

general filters:

--nonvar|-v	process nonvariant sites (in addition to variants)
--keep_nonbinary|-2	keep non-biallelic sites
--exclude_contigs|-r	FILE	list of contigs/chromosomes to exclude (each line should list a contig name exactly as it appears in the input vcf file)
--rmv_nonexonic|-t	filter-out non-exonic sites (requires -X)

output:

--bed|-B	FILE	name of dumped BED format file for sites that pass all filters
--failed_sites|-p	FILE	name of dumped file containing sites that failed at least one filter (bziped)
--vcfout|-o	FILE	name of dumped vcf file containing sites that passed filters
--ind_depth|-I	FILE	dumped file with individual mean depth

##########

Notes:

Some of the filters rely on annotations generated by SAMtools/BCFtools.
To use the eveness-of-coverage filters (options -k and -u), -D must be used with satmools mpileup.
If option -s or -S is set to 0, that particular strand bias filter is not performed.
It is recomended to use mpileup -I to ignore indels.
Characters in front of filtered sites (dumped with option -p) indicate filters that the site failed to pass.

```

So there are a lot of options.
Input file is a VCF.
Some of these options require the mpileup file, which we have already produced.
Moreover, we may need to give an ancestral sequence (indexed, `input/lyca/referenceseq.fasta`) for some mutation-type filters (e.g. for ancient DNA).


-----

**WORKFLOW**:
... > FILTERING SITES > DEPTH

First we want to remove sites with extremely low or high depth.
Let us get the empirical distribution.

```
cat output/lyca.raw.vcf | perl ../scripts/SNPcleaner.pl -v -a 0 -k 10 -u 2 -d 1 -S 0 -s 0 -G output/lyca.mpileup -M NN_NN -o output/lyca.filt.vcf
# it should print how many sites were processed and the individual mean depths
```

Let us get the distribution of per-site depth values (script by F.G Vieira):
```
cat output/lyca.filt.vcf | tr ";" "\t" | awk 'BEGIN{print "chr_pos\tdepth"} !/#/{sub("DP=","",$8); if(rand()<=1) print $1"_"$2"\t"$8}' > output/lyca.sdepth
Rscript scripts/plotDepth.R output/lyca.sdepth output/lyca.sdepth.pdf
open output/lyca.sdepth.pdf # on mac
```

You can investigate the distribution of per-site depth:
```
     0%    2.5%      5%    7.5%     10%   12.5%     15%   17.5%     20%   22.5% 
  42.00   43.00   53.00   54.00   61.00   63.00   63.00   72.00   77.00   84.00 
    25%   27.5%     30%   32.5%     35%   37.5%     40%   42.5%     45%   47.5% 
 116.00  151.00  154.00  157.00  163.00  176.00  178.00  180.00  183.00  183.00 
    50%   52.5%     55%   57.5%     60%   62.5%     65%   67.5%     70%   72.5% 
 191.00  196.00  236.65  250.00  274.00  303.00  313.00  327.00  334.00  371.00 
    75%   77.5%     80%   82.5%     85%   87.5%     90%   92.5%     95%   97.5% 
 397.00  399.00  571.00  655.00  667.00  716.25  748.00  873.00 1067.00 1131.00 
   100% 
1576.00
```
and decide which cutoffs are more suitable for your scope.


---------------

**WORKFLOW**:
... > FILTERING SITES > BIAS

We are now ready to perform the first round of filtering.
Here some parameters:

Parameter | Value | Description
--------- | ----- | -----------
-G | output/lyca.mpileup | mpileup file <br> 
-P | input/lyca/clst.txt | this gives the correspondence sample-pop for HWE filtering <br>
-p | output/lyca.bad.bz2 | list of filtered sites <br>
-d | 50 | minimum site depth <br>
-D | 900 | maximu site depth <br>
-k,-u | 10,2 | at least 10 samples must have 2 reads <br>
-a | 0 | minimum number of alternate alleles per site <br>
-Q | 20 | minimum RMS mapping quality for SNPs <br>
-f | 0.001 | min p-value for map quality bias <br>
-S | 0.001 | min p-value for strand bias from combining p-values across individuals <br>
-T | 0.001 | min p-value for allele bias in potential heterozygotes <br>
-b | 0.001 | min p-value for base quality bias <br>
-e | 0.001 | min p-value for biased distance of alternate bases from ends of reads <br>
-h | 0.001 | min p-value for exact test of HWE <br>
-v | process nonvariant sites <br>

Let us write down the command:
```
OPT1='-a 0 -k 10 -u 2 -d 50 -D 900 -Q 20'
OPT2='-f 0.001 -S 0.001 -s 0 -T 0.001 -b 0.001 -e 0.001 -h 0.001 -v'
cat output/lyca.raw.vcf | perl scripts/SNPcleaner.pl $OPT1 $OPT2 -G output/lyca.mpileup -P input/lyca/clst.txt -M NN_NN -o output/lyca.filt.vcf -p output/lyca.bad.bz2
```

How many sites did we filter? How many sites were retained?
We can inspect the removed sites to better tune our filtering parameters.
```
bunzip2 output/lyca.bad.bz2
less -S output/lyca.bad 
```

Here some examples of filtered sites:

* low depth:
```
kd	ref_contig	50160	.	A	.	23.7	.	DP=1;;AC1=25;FQ=-23.7	PL:DP:SP	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0
```

* N in the reference
```
Nkd	ref_contig	81199	.	N	.	5.42	.	DP=1;AF1=1;AC1=40;DP4=0,0,1,0;MQ=50;FQ=-23.7	PL:DP:SP	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	26:1:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0	0:0:0
```

* HWE test:
```
kdh(p=0.000108018128040375;0)	ref_contig	51494	.	C	A	8.93	.	\
DP=14;VDB=5.960000e-02;RPB=8.745357e-01;AF1=0.6656;AC1=27;DP4=1,0,2,0;MQ=37;FQ=3.43;PV4=1,0.14,0.11,1	GT:PL:DP:SP:GQ	0/1:0,0,0:0:0:3	0/1:0,0,0:0:0:3	0/1:0,0,0:0:0:3	0/1:0,0,0:0:0:3	0/1:20,3,0:1:0:4	0/1:26,3,0:1:0:4	0/0:0,3,34:1:0:4	0/1:0,0,0:0:0:3	0/1:0,0,0:0:0:3	0/1:0,0,0:0:0:3	0/1:0,0,0:0:0:3	0/1:0,0,0:0:0:3	0/1:0,0,0:0:0:3	0/1:0,0,0:0:0:3	0/1:0,0,0:0:0:3	0/1:0,0,0:0:0:3	0/1:0,0,0:0:0:3	0/1:0,0,0:0:0:3	0/1:0,0,0:0:0:3	0/1:0,0,0:0:0:3
```
Indeed here all individuals are heterozygotes.

As you can see we have many sites that failed `S` test (strand bias) so maybe we have been too strict.
For other tests, it seems that all sites passed so maybe we have been too relaxed.
Let us check the proportion of sites filetered out by each test.
```
Rscript ../scripts/plotFilters.R output/lyca.bad output/lyca.bad.pdf 7000
open output/lyca.bad.pdf # on mac
```

We can rerun it by adjusting our parameters from the considerations.
```
OPT1='-a 0 -k 7 -u 1 -d 20 -D 900 -Q 20'
OPT2='-f 0.001 -S 1e-100 -s 0 -T 0.05 -b 0.05 -e 0.05 -h 0.05 -v'
cat output/lyca.raw.vcf | perl scripts/SNPcleaner.pl $OPT1 $OPT2 -B output/lyca.good.bed -G output/lyca.mpileup -P input/lyca/clst.txt -M NN_NN -o output/lyca.filt.vcf -p output/lyca.bad2.bz2
bunzip2 output/lyca.bad2.bz2
```

Let us check sites that failed due to quality bias:

* base quality bias
```
b	ref_contig	50209	.	C	.	24.8	.	DP=39;RPB=1.633743e+00;AF1=0.02792;AC1=1;DP4=34,0,1,0;MQ=44;FQ=-24.7;PV4=1,0.024,0.031,1	PL:DP:SP	0:2:0	0:1:0	0:5:0	0:0:0	0:2:0	0:4:0	0:1:0	0:1:0	0:2:0	0:1:0	0:2:0	0:1:0	0:2:0	15:4:0	0:1:0	0:0:0	0:1:0	0:0:0	0:1:0	0:4:0
grep 50209 output/lyca.mpileup 
ref_contig	50209	C	2	..	GH	1	.	H	5	.....	HGIIE	0	*	*	2	..	GI	4	....	IHHI	1	..	HG	1	.	I	2	..	6G	1	.	I	2	..	GG	4	...T	G5D<	1	.	D	0	*	*	....	EHHI
```

* map quality bias
```
f	ref_contig	50222	.	T	.	8.03	.	DP=39;RPB=1.528942e+00;AF1=0.04621;AC1=1;DP4=33,0,1,0;MQ=44;FQ=-7.86;PV4=1,1,3.9e-05,1	PL:DP:SP	0:2:0	0:1:0	0:5:0	0:0:0	0:2:0	0:5:0	0:1:0	0:0:0	0:2:0	0:1:0	0:2:0	0:1:0	0:1:0	0:4:0	0:1:0	0:0:0	36:1:0	0:0:0	0:1:0	0:4:0
grep 50222 output/lyca.mpileup 
ref_contig	50222	T	2	..	GI	1	.	G	5	.....	HIIEH	0	*	*	2	..	>I	5	.....	IEHGI	1	..	G:	1	.	G	2	..	HG	1	.	F	1	.	I	4	....	>H@E	1	.	G	0	*	*	....	EDFI
```

* strand bias
```
DS	ref_contig	55250	.	C	.	82.8	.	DP=567;AF1=0;AC1=0;DP4=551,0,0,0;MQ=46;FQ=-82.7	PL:DP:SP	0:31:0	0:36:0	0:34:0	0:22:0	0:42:0	0:27:0	0:33:0	0:33:0	0:15:0	0:18:0	0:28:0	0:22:0	0:27:0	0:14:0	0:18:0	0:21:0	0:33:0	0:25:0	0:37:0	0:35:0
mac:try matteo$ grep 55250 output/lyca.mpileup
ref_contig	55250	C	31	...............................	IGI7EIHHDI@E?GBIGHIIIIIIIDIHHGI	36	....................................	HGGIEHEAHIEIGHGIBH>GIHFGGBH=HHGDIHHH	34	..................................	IIIHIHHIID=GHIIHI@GICI;EID;HAIIHHE	22	......................	BHDGIDIIHBDIGHGHGGHH?C	42	..........................................	@=IHHIHIHIIDHIIIDDHHIGIGDIHHIHGII7IGDEI@II	27	...........................	IEIHDII5BHIIHIIDHHHGIAHIIIH	33	.................................	<IGIIIIEHGGHIE>IHIHGGGGIDIGIGIHHD	33	.................................	HBIHIIGIIIIIIHIHGH?III?>BIH@D=IIH	15	...............	GHGHIIGIIGIHBGG	18	..................	>HHIGGDEHHIIHIHGII	28	............................	HGHGEHIHIDHIHIGIIIIH?III=HII	22	......................	IIGIGDBGIHBIIHIIBIIGFI	27	...........................	8II<?DHIGDIHHGIHFIHIIHADHIH	14	..............	IFGHDGDGHHHHDH	18	..................	IIDBIDGBGCIH>IDGII21	.....................	IGFHIH<IHIEIIGIDIHGIH	33	.................................	HGIBIIGFHIGHHIHDIDHHDIHHIHDHIHHH;	25	.........................	DIBHIIBDEDDH:BADIGIBDH:II	37	.....................................	IHIIGHHIGGIGHIAHHHGGIIHIHIBIH7IGHHFDI	35	...................................	?FIIBIEIIIHG?IIIHHIIGIAEIGHDDHIHHIH
```

-------------

**WORKFLOW**:
... > FILTERING SITES > CHECK

Now we check the resulting Site Frequency Spectrum (SFS).
We will see later how to accurately estimate it from sequencing reads.
So far, we will quickly compute it from the resulting VCF file.

```
perl -ne 'print "$1\n" if /AF1=([^,;]+)/' output/lyca.filt.vcf > output/lyca.vcf.freqs
Rscript scripts/plotSFS.R output/lyca.vcf.freqs output/lyca.vcf.sfs.pdf
open output/lyca.vcf.sfs.pdf # on mac
```
You can also get it [here](http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/lyca.vcf.afs.pdf).

The output for this filtering step is a BED file (along with the filtered VCF file) with the positions that passed the quality filters. 
```
head output/lyca.good.bed
```
This file can be used for all downstream analyses, along with the final VCF with only the valid sites.



