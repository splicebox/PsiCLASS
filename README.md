PsiCLASS
=======

Described in: 

Song, L., Sabunciyan, S., and Florea, L. (2018). A multi-sample approach increases the accuracy of transcript assembly, *Submitted*. 

	Copyright (C) 2018- and GNU GPL by Li Song, Liliana Florea

Includes portions copyright from: 

	samtools - Copyright (C) 2008-, Genome Research Ltd, Heng Li
	
Commands, scripts and supporting data for the paper can be found [here](https://github.com/splicebox/PsiCLASS_paper/).

### What is PsiCLASS?

PsiCLASS is a reference-based transcriptome assembler for single or multiple RNA-seq samples. Unlike conventional methods that analyze each sample separately and then merge the outcomes to create a unified set of meta-annotations, PsiCLASS takes a multi-sample approach, simultaneously analyzing all RNA-seq data sets in an experiment. PsiCLASS is both a transcript assembler and a meta-assembler, producing  separate transcript sets for the individual samples and a unified set of meta-annotations. The algorithmic underpinnings of PsiCLASS include using a global subexon splice graph, statistical cross-sample feature (intron, subexon) selection methods, and an efficient dynamic programming algorithm to select a subset of transcripts from among those encoded in the graph, based on the read support in each sample. Lastly, the set of meta-annotations is selected from among the transcripts generated for individual samples by voting. While PsiCLASS is highly accurate and efficient for medium-to-large collections of RNA-seq data, its accuracy is equally high for small RNA-seq data sets (2-10 samples) and even for individual samples, and therefore can be effectively used both as a multi-sample and as a conventional single-sample assembler. 

### Install

1. Clone the [GitHub repo](https://github.com/splicebox/psiclass), e.g. with `git clone https://github.com/splicebox/psiclass.git`
2. Run `make` in the repo directory

You will find the executable files in the downloaded directory. If you want to run PsiCLASS without specifying the directory, you can either add the directory of PsiCLASS to the environment variable PATH or create a soft link ("ln -s") of the file "psiclass" to a directory in PATH.

PsiCLASS depends on [pthreads](http://en.wikipedia.org/wiki/POSIX_Threads) and samtools depends on [zlib](http://en.wikipedia.org/wiki/Zlib).


### Usage

	Usage: ./psiclass [OPTIONS]
		Required:
                	-b STRING: paths to the alignment BAM files; use comma to separate multiple BAM files
                        	or
               		--lb STRING: path to the file listing the alignment BAM files
       		Optional:
               		-s STRING: path to the trusted splice file (default: not used)
          		-o STRING: prefix of output files (default: ./psiclass)
          		-p INT: number of threads (default: 1)
			-c FLOAT: only use the subexons with classifier score <= than the given number. (default: 0.05)
			-vn INT |-vf FLOAT : minimum number | fraction of samples a transcript must appear in to be reported 
				(defaults: 3, 0.25)
              		--stage INT:  (default: 0)
                     		0-start from the beginning - building the splice site file for each sample
                     		1-start from building the subexon file for each samples
                     		2-start from combining the subexon files across samples
                     		3-start from assembling the transcripts for each sample
                     		4-start from voting the consensus transcripts across samples
	
### Practical notes

*Alignment compatibility.* PsiCLASS has been tuned to run on alignments generated with the tools [HISAT](https://ccb.jhu.edu/software/hisat/index.shtml) and [STAR](https://github.com/alexdobin/STAR). 

When running PsiCLASS with STAR alignments, run STAR with the option `--outSAMstrandField intronMotif`, which will include the XS field indicating the strand in the BAM alignments. Further, when including alignments with *non-canonical splice sites*, use the provided `addXS` executable to add the XS field:

	samtools view -h in.bam | ./addXS reference_genome.fa | samtools view -bS - > out.bam

*Trusted introns from other sources.* By default, PsiCLASS determines a set of trusted introns from the input spliced alignments, to use in building the global subexon graph. Alternatively, the user can supply an external set of trusted introns, for instance extracted from the GENCODE gene annotations or judiciously selected from the input data using a tool like [JULIP](https://github.com/Guangyu-Yang/JULiP). This file must contain three columns:

	chr_name start_site end_site
	
*Voting optimization.* The default parameters for voting may not be optimal for all types of data, for instance a lower voting cutoff may be more appropriate for sparse rRNA depleted total RNA samples. To determine a better cutoff value, one can run the voting stage (see [Usage](#usage) above) with different cutoffs, and assess the performance against a reference set of gene annotations, such as [GENCODE](https://www.gencodegenes.org). Note that the per sample sets of transcripts will remain unchanged.        

*Add gene name.* For many experiments, it would be convenient to have the canonical gene name to associate with each transcript. PsiCLASS provides the program "add-genename" for such purposes. "add-genename" takes the GTF of gene annotation and a file listing the raw GTF files as input. It will generate a new GTF file for each input raw GTF file by appending the canonical gene names. If a gene is not found in the annotation, "add-genename" will use "novel_INT" to represent its gene name. 

### Input/Output

The primary input to PsiCLASS is a set of BAM alignment files, one for each RNA-seq sample in the analysis. The program calculates a set of subexon files and a set of splice (intron) files, for the individual samples. (Optionally, one may specify a path to an external file of trusted introns as explained [above](#practical-notes).) The output consists of one GTF file of transcripts for each sample, and the GTF file of meta-annotations produced by voting, stored in the output directory:

	Sample-wise GTF files: (psiclass)_sample_{0,1,...,n-1}.gtf
	Meta-assembly GTF file: (psiclass)_vote.gtf

where indices 0,1,...,n-1 match the order of the input BAM files.

Subexon and splice (intron) files, and other auxiliary files, are in the subdirectories:

	Intron files: splice/*
	Subexon graph files: subexon/*
	Log file: (psiclass)_classes.log

### Example

The directory './example' in this distribution contains two BAM files, along with an example of a BAM list file. Run PsiCLASS with:

	./psiclass -b example/s1.bam,example/s2.bam

or

	./psiclass --lb example/slist

The run will generate the files 'psiclass_sample_0.gtf' for 's1.bam', 'psiclass_sample_1.gtf' for 's2.bam', and the file 'psiclass_vote.gtf' containing the meta-assemblies.

### Terms of use

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received (LICENSE.txt) a copy of the GNU General
Public License along with this program; if not, you can obtain one from
http://www.gnu.org/licenses/gpl.txt or by writing to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
### Support

Create a [GitHub issue](https://github.com/splicebox/PsiCLASS/issues).
