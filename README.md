
# EvoGen_course

Evolutionary Genomics course - Data Analysis module

13th-17th April 2015

## Material

Instructions to download all material necessary for the practical sessions are given [here](https://github.com/mfumagalli/EvoGen_course/blob/master/install.md).

## Agenda

### Day 1

 *	Lecture: From raw NGS data to genotypes ([pdf](https://github.com/mfumagalli/EvoGen_course/blob/master/slides_day_1.pdf))

 *	Paper discussion: slides ([pdf](https://github.com/mfumagalli/EvoGen_course/blob/master/slides_day_1_paper.pdf)) and paper ([pdf](http://cteg.berkeley.edu/~nielsen/wordpress/wp-content/uploads/2013/01/Nielsen-R.-et-al.-2011.pdf))

 *	Practical session: 
	+	[data filtering](https://github.com/mfumagalli/EvoGen_course/blob/master/filtering.md)
	+       [genotype calling](https://github.com/mfumagalli/EvoGen_course/blob/master/genocall.md)

 *	Additional material:

	+       genotype calling using:
		-	[SAMtools](https://github.com/mfumagalli/EvoGen_course/blob/master/genocall_samtools.md)
		-	[BEAGLE](https://github.com/mfumagalli/EvoGen_course/blob/master/imputation.md)

### Day 2

 *	Lecture: SNP calling and advanced methods for evolutionary inferences from NGS data ([pdf](https://github.com/mfumagalli/EvoGen_course/blob/master/slides_day_2.pdf))

 *	Paper discussion: slides ([pdf](https://github.com/mfumagalli/EvoGen_course/blob/master/slides_day_2_paper.pdf)) and papers ([link](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0079667)) ([link](http://journal.frontiersin.org/article/10.3389/fgene.2012.00066/abstract))

 *	Practical session: 

	+       [SNP calling](https://github.com/mfumagalli/EvoGen_course/blob/master/snpcall.md)
	+	advanced methods to estimate [SFS](https://github.com/mfumagalli/EvoGen_course/blob/master/sfs.md)

 *	Additional material: 

	+       SNP calling using [SAMtools](https://github.com/mfumagalli/EvoGen_course/blob/master/snpcall_samtools.md)
	+	estimation of [inbreeding](https://github.com/mfumagalli/EvoGen_course/blob/master/inbreeding.md) coefficients
	+	advanced methods to calculate summary statistics using:
		-	[ANGSD](https://github.com/mfumagalli/EvoGen_course/blob/master/lowcov.md)
		-	[ngsTools](https://github.com/mfumagalli/EvoGen_course/blob/master/lowcov_ngstools.md)
		-	[ngsTools/ngsDist](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md)



### Day 3

 *	Lecture: Population structure and demographic inferences ([pdf](https://github.com/mfumagalli/EvoGen_course/blob/master/slides_day_3.pdf))

 *	Journal discussion: slides ([pdf](https://github.com/mfumagalli/EvoGen_course/blob/master/slides_day_3_paper.pdf)) and paper ([pdf](https://github.com/mfumagalli/EvoGen_course/blob/master/Moltke_AJHG_2015.pdf))

 *	Practical exercises:

	+	estimating demographic parameters through simulations ([R script](https://github.com/mfumagalli/EvoGen_course/blob/master/practise_day_3.R))
	+	inference of population splits and mixture events using TreeMix ([bash script](https://github.com/mfumagalli/EvoGen_course/blob/master/practise_day_3_extra.txt))

### Day 4

 *	Lecture: Detecting natural selection ([pdf](https://github.com/mfumagalli/EvoGen_course/blob/master/slides_day_4.pdf))

 *	Journal discussion: slides ([pdf](https://github.com/mfumagalli/EvoGen_course/blob/master/slides_day_4_paper.pdf)) and paper ([pdf](https://github.com/mfumagalli/EvoGen_course/blob/master/Liu_Cell_2014.pdf))

 *	Practical exercises:

	+	selection scan in the human genome and assessment of significance ([R script](https://github.com/mfumagalli/EvoGen_course/blob/master/practise_day_4.R))
	+	summary statistics using [ngsTools/ngsDist](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md)

	


## Credits

Some materials have been borrowed (and then adapted) from [Thorfinn Korneliussen](http://scholar.google.co.uk/citations?user=-YNWF4AAAAAJ&hl=en), [Anders Albrechtsen](http://popgen.dk/albrecht/web/WelcomePage.html), [Tyler Linderoth](http://scholar.google.com/citations?user=dTuxmzkAAAAJ&hl=en), [Filipe G. Vieira](http://scholar.google.com/citations?user=gvZmPNQAAAAJ&hl=en).
Sequencing data on butterfiles has been kindly provided by [Zach Gompert](https://gompertlab.wordpress.com/).
ANGSD has been developed by Thorfinn Korneliussen and Anders Albrechtsen. 
Filipe G. Vieira implemented the inbreeding estimation. 
Tyler Linderoth wrote most of data filtering scripts together with [Sonal Singhai](https://systemsbiology.columbia.edu/people/sonal-singhal), [Ke Bi](http://scholar.google.ca/citations?user=ymcwERQAAAAJ), and Filipe G. Vieira.








