
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





