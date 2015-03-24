
**WORKFLOW**:
LAB WORK > SEQUENCING > LOW-LEVEL DATA (reads) <br>

Sequencing machines produce millions of DNA sequences.
Some basic quality control checks are essential to avoid that your biological conclusions are not biased by some unusual features of the data.

*** Filtering reads

**WORKFLOW**:
LAB WORK > SEQUENCING > LOW-LEVEL DATA (reads) > FILTERING READS

First step is to filter our and/or trim reads in order to remove low quality data.
Cleaning reads is an important step to decrease the likelihood that alignment and sequencing errors are mistaken as SNPs.

There are several ways and tools one can use to achieve this goal.
Here we will use some custom scripts (mostly borrowed from F.G. Vieira and T. Linderoth) and several external programs.
We will also use FastQC (available [[http://www.bioinformatics.babraham.ac.uk/projects/fastqc/][here]]).

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



