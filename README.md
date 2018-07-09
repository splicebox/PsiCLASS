PsiCLASS
=======

Described in: 

Song, L., Sabunciyan, S., and Florea, L. (2018). A multi-sample approach increases the accuracy of transcript assembly, *Submitted*.

Copyright (C) 2018- and GNU GPL by Li Song, Liliana Florea

Includes portions copyright from:
[`samtools`](https://github.com/samtools/samtools) - Copyright (C) 2008-, Genome Research Ltd, Heng Li

### What is PsiCLASS?

PsiCLASS is a reference-based transcriptome assembler for single or multiple RNA-seq samples. Unlike conventional methods that analyze each sample separately and then merge the outcomes to create a unified set of meta-annotations, PsiCLASS takes a multi-sample approach, simultaneously analyzing all RNA-seq data sets in an experiment. PsiCLASS is both a transcript assembler and a meta-assembler, producing  separate transcript sets for the individual samples and a unified set of meta-annotations. The algorithmic underpinnings of PsiCLASS include using a global subexon splice graph, statistical cross-sample feature (intron, subexon) selection methods, and an efficient dynamic programming algorithm to select a subset of transcripts from among those encoded in the graph, based on the read support in each sample. Lastly, the set of meta-annotations is selected from among the transcripts generated for individual samples by voting. While PsiCLASS is highly accurate and efficient for medium-to-large collections of RNA-seq data, its accuracy is equally high for small RNA-seq data sets (2-10 samples) and even for individual samples, and therefore can be effectively used both as a multi-sample and as a conventional single-sample assembler. 

### Install

1. Clone the [GitHub repo](https://github.com/splicebox/psiclass), e.g. with `git clone https://github.com/splicebox/psiclass.git`
2. Run `make` in the repo directory

PsiCLASS depends on [pthreads](http://en.wikipedia.org/wiki/POSIX_Threads) and samtools depends on [zlib](http://en.wikipedia.org/wiki/Zlib).


### Usage

#### Basic usage

	Usage: ./psiclass [OPTIONS]
		Required:
                	-b STRING: paths to the alignment BAM files; use comma to separate multiple BAM files
                        	or
               		--lb STRING: path to the file listing the alignment BAM files
       		Optional:
               		-s STRING: path to the trusted splice sites file (default: not used)
          		-o STRING: prefix of output files (default: ./psiclass)
          		-t INT: number of threads (default: 1)
              		--stage INT:  (default: 0)
                    		0-start from the beginning - building the splice site file for each sample
                   		1-start from building the subexon file for each samples
                     		2-start from combining the subexon files across samples
                     		3-start from assembling the transcripts for each sample
                     		4-start from voting the consensus transcripts across samples
	
#### Advanced usage

Alternatively, one can run the different components of the program in succession, for instance:

### Practical notes

*Alignment compatibility.* PsiCLASS has been tuned to run on alignments generated with the tools [HISAT](https://ccb.jhu.edu/software/hisat/index.shtml) and [STAR](https://github.com/alexdobin/STAR). At this point, its behavior on alignment data generated with other tools is unknown. 

*Voting optimization.* The default parameters for voting may not be optimal for all types of data, for instance a lower voting cutoff may be more appropriate for sparse rRNA depleted total RNA samples. To determine a better cutoff value, one can run the voting tool (see [Advanced usage](#advanced-usage) above) with different cutoffs, and assess the performance against a reference set of gene annotations, such as [GENCODE](https://www.gencodegenes.org), using the included tool 'grader'. Note that the per sample sets of transcripts will remain unchanged.        

### Input/Output

The primary input to PsiCLASS is a set of BAM alignment files, one for each RNA-seq sample in the analysis. The program calculates a set of subexons files and a set of splice (intron) files, for the individual samples. (Optionally, one may specify a path to an external file of trusted introns.) The output consists of one GTF file of transcripts for each sample, and the GTF file of meta-annotations produced by voting, stored in the output directory:

	Sample-wise GTF files: (psiclass)_sample_{0,1,...,n-1}.gtf
	Meta-assembly GTF file: (psiclass)_vote.gtf

where indices 0,1,...,n-1 match the order of the input BAM files.

Subexon and splice (intron) files, and other auxiliary files, are in the subdirectories:

	Intron files: splice/*
	Subexon graph files: subexon/*
	Log file: ./(psiclass)_classes.log

### Example

The directory './example' in this distribution contains two BAM files, along with an example of a BAM list file. You can run PsiCLASS with:

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
