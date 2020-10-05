# 1. INSTALLING aScan

aScan depends on the following libraries that should be installed in your system prior to compiling aScan:

1) GSL - GNU Scientific Library https://www.gnu.org/software/gsl/
2) zlib - https://zlib.net/
3) boost C++ libraries - https://www.boost.org/ 
4) bamtools https://github.com/pezmaster31/bamtools

For convenience, the bamtools library source is included in the aScan package, build it as follows.

# 1.1 BUILDIING bamtools

The cmake utility (https://cmake.org/) is required to compile the bamtools library. 

From the aScan folder type:

cd bamtools

mkdir build

cd build

cmake ..

make

For further support on installing the bamtools refer to https://github.com/pezmaster31/bamtools/wiki/Building-and-installing

# 1.2 BUILDING aScan

To compile aScan run the following command from within the aScan folder

g++ aScan_dev.cpp -o aScan -I bamtools/src/ -lbamtools -L bamtools/build/src/api -lgsl -lgslcblas -lz -std=c++11 -Wall -pthread -O2

This command assumes you compiled the bamtools library as described in 1.1. If you prefer to point to another bamtools installation, just modify the command to match your bamtools path.

You can now move the resulting aScan executable to your local bin PATH to remove the need to specify the path to the executable each time you use it.

# 2. USING aScan

aScan is a user-friendly bioinformatics tool devised expressly to reliably detect allele-specific expression from genomic and transcriptomic data by the same individual. 
As such, the primary input consists of VCF and BAM files containing matched genomic variation and RNA-Seq data. A reference GTF transcripts annotation file is also required and, finally, an optional transcript-ID/gene-name association file can be supplied. 

# 2.1 Synopsis

aScan --rna rnaseq_bamfile --vcf vcf_file --gtf gtf_file --only_first_mate [-nt transcript/gene_correspondence_file] [-p thread_num] [--filter filter_key_word]

--rna rnaseq_bamfile 
BAM file containing RNA-Seq reads mapped to a reference genome. The BAM file does not need to be sorted. Please notice that aScan does not currently support reads mapped to the reference transcriptome. If your BAM file has been mapped to the transcriptome, you can use the rsem-tbam2gbam utility (https://deweylab.github.io/RSEM/) and convert it to genome coordinates.

--vcf vcf_file: VCF file with genomic variants from the same individual as the RNA-Seq data. Please notice that aScan currently does not support multi-VCF files.

--gtf gtf_file: GTF file with reference transcripts annotation. Please notice that aScan will consider only "exon" named features. The 


