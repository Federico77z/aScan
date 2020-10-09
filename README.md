# 1. INSTALLING aScan

A C++11 compiler (e.g., g++) is needed to build aScan.  

aScan depends on the following libraries that should be installed in your system before compiling aScan:

1) GSL - GNU Scientific Library https://www.gnu.org/software/gsl/
2) zlib - https://zlib.net/
3) boost C++ libraries - https://www.boost.org/ 
4) bamtools https://github.com/pezmaster31/bamtools

For convenience, the bamtools library source is included in the aScan package, build it as follows.

If you prefer, you can skip all the installation steps and use instead the aScan Docker image (see 1.3).

# 1.1 BUILDING bamtools

The cmake utility (https://cmake.org/) is required to compile the bamtools library. 

From the aScan folder type:

cd bamtools

mkdir build

cd build

cmake ..

make

For further support on installing the bamtools refer to https://github.com/pezmaster31/bamtools/wiki/Building-and-installing.

# 1.2 BUILDING aScan

To compile aScan (using g++) run the following command from within the aScan folder:

g++ aScan_dev.cpp -o aScan -I bamtools/src/ -lbamtools -L bamtools/build/src/api -lgsl -lgslcblas -lz -std=c++11 -Wall -pthread -O2

That assuming you compiled the bamtools library as described in 1.1. If you prefer to point to another bamtools installation, just modify the command line to match your bamtools path.

You can now move the resulting aScan executable to your local bin PATH to remove the need to specify the path to the executable each time you use it.

# 1.3 aScan Docker image

A docker image for aScan is available at https://quay.io/repository/pmandreoli/ascan. See the aScanDocker repo at https://github.com/pmandreoli/aScanDocker for further info.

# 2. USING aScan

aScan is a user-friendly bioinformatics tool devised expressly to reliably detect allele-specific expression from genomic and transcriptomic data by the same individual. 
As such, the primary input consists of VCF and BAM files containing matched genomic variation and RNA-Seq data. A reference GTF transcripts annotation file is also required and, finally, an optional transcript-ID/gene-name association file can be supplied. 

# 2.1 SYNOPSIS

aScan --rna rnaseq_bamfile --vcf vcf_file --gtf gtf_file [--nametable transcript/gene_correspondence_file] [-p thread_num] [--filter filter_key_word]

--rna rnaseq_bamfile 
BAM file containing RNA-Seq reads mapped to a reference genome. The BAM file does not need to be sorted. Please notice that aScan does not currently support reads mapped to the reference transcriptome. If your BAM file has been mapped to the transcriptome, you can use the rsem-tbam2gbam utility (https://deweylab.github.io/RSEM/) to convert it to genome coordinates.

--vcf vcf_file 
VCF file with genomic variants from the same individual as the RNA-Seq data. Please notice that aScan currently does not support multi-VCF files. aScan considers only single nucleotide substitutions. aScan can work both with phased or unphased VCF files.

--gtf gtf_file 
GTF file with reference transcripts annotation. Please notice that aScan will consider only "exon" named features. The attribute column of the GTF file must include a transcript_id value for each exon.

-nt, --nametable transcript/gene_correspondence_file
A tab-separated tabular file associating transcript IDs (first column) to gene names or gene IDs (second column).

-p threadnum
The number of threads to use (default 1). While aScan supports multithreading, the execution speed bottleneck is usually represented by the drive read bandwidth, meaning that normally the benefits of using more than 2 or 3 threads are minimal. 

--filter filter_key_word
Using this option, all the variants flagged with filter_key_word in the FILTER field of the VCF file will be ignored. Use this option to discard low-quality variant calls.

# 2.2 EXAMPLE AND TEST DATASET

A small example dataset is included in the test_dataset folder of the aScan package. The dataset can also be used as a test to check that your aScan build works correctly.

The example can be run as follows:

cd test_dataset
aScan --rna heart_S12_chr19.bam --vcf S12_def_chr19.vcf --gtf hg19_chr19.gtf -nt transcript-gene_table.txt

The produced output file heart_S12_chr19.bam_aScan.txt should be mostly equal to the heart_S12_chr19.bam_aScan_reference.txt included for reference. Tiny differences in the FDR values can be expected on different systems. Also, the rank position of genes with identical p-values (see OUTPUT) may slightly differ among different aScan runs.

# 3.0 OUTPUT

The output file consists of a tab-separated tabular file with nine fields.  

GENE_NAME	the name of gene or transcript.

HZ_POS_N	number of heterozygous positions overlapping the gene/transcript according to the provided VCF file

HZ_POS		genomic coordinates of the heterozygous positions overlapping the gene/transcript.

H_HAPLO		high expression haplotype as reconstructed by aScan (0 for reference allele, 1 for alternative allele) for unphased VCFs.  
		First parental haplotype (0 for reference allele, 1 for alternative allele) for phased VCFs.

L_HAPLO		low expression haplotype as reconstructed by aScan (0 for reference allele, 1 for alternative allele) for unphased VCFs.
		Second parental haplotype (0 for reference allele, 1 for alternative allele) for phased VCFs.

H_COV		Total number of reads supporting the high expression haplotype for unphased VCFs or the first parental haplotype for phased VCFs.

L_COV		Total number of reads supporting the low expression haplotype for unphased VCFs or the second parental haplotype for phased VCFs.

PV		allele-specific expression p-value as computed by aScan.

FDR		False Discovery Rate associated with the allele-specific expression p-value.

Genes are sorted by ascending p-value. Only genes/transcripts overlapping with at least one heterozygous substitution are reported in the output.

Genes with an associated FDR lower than an arbitrary threshold (usually 0.01 or 0.05) can be considered as showing a genuine allele-specific expression bias.
 
