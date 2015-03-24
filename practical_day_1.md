
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

- [per base quality](http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/per_base_quality.pdf][per base quality)

- [per sequence quality](http://palin.popgen.dk/mfumagalli/Workshop/ANU/web/per_sequence_quality.pdf][per sequence quality)

-------------------







