# CircParser
CircParser: a novel streamlined pipeline for circular RNA structure and host gene prediction in non-model organisms

  CircParser is simple and fast pipeline that uses the outputs from the most common circular RNAs in silico prediction tools (CIRI, CIRI2, CircExplorer2, find_circ, and circFinder) to annotate circular RNAs, assigning
presumable host genes from local or public databases such as National Center for Biotechnology Information (NCBI). Also this pipeline can discriminate circular RNAs
based on their structural components (exonic, intronic, exon-intronic or intergenic) using genome annotation file.


## Dependencies and requirements
Operating system: Linux only.

You need:
- Samtools: We recommend the newests versions of SAMtools (e.g. > 1.4.1)
- Bedtools
- Blast

## Installation
Install dependencies, you can skip this step if these packages are already installed on your system

- Samtools 

            1. wget https://github.com/samtools/samtools/releases/download/1.4.1/samtools-1.4.1.tar.bz2 -O samtools.tar.bz2
            2. tar -xjvf samtools.tar.bz2 
            3. cd samtools-1.4.1/
            4. ./configure
            5. make
            6. make install
 
- Bedtools 

            1. wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz
            2. tar -zxvf bedtools-2.28.0.tar.gz
            3. cd bedtools2
            4. ./configure
            5. make
            
- Blast 

Get the compiled executables from this URL:

```
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```

Decompress the archive. For example:

```
tar xvfz ncbi-blast-2.9.0+-x64-linux.tar.gz
```

Add the `bin` folder from the extracted archive to your path. For example, add
the following line to your `~/.bashrc` file:

```
export PATH="/PATH/TO/ncbi-blast-2.9.0+/bin":$PATH
```

And change the `/PATH/TO` part to the path where you have put the extracted
archive.

## Download CircParser

Download the program like this:
```
git clone https://github.com/SharkoTools/CircParser.git
```           

## Usage


## Command line
-  -b              circRNA input file (required)
-  -g,--genome     reference genome file (required)
-  -t,--tax        NCBI TaxID (optional)
-  -a              genome annotation file, gff/gff3 file (optional)
-  -c,--ciri       input circRNA from CIRI|CIRI2 in silico predictors, (default: input from CircExplorer2, find_circ, circFinder, and BED files)
-  --np            prohibition for coordinates merging (optional)
-  --threads       number of threads (CPUs) for BLAST search (default: 8)
-  -h,--help       show this help message and exit
-  --version,-v    current version
