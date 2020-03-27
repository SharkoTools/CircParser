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

### Schema
![Image alt](https://github.com/SharkoTools/CircParser/blob/master/Figure_1.png)


### How to run
```
perl CircParser.pl -b INPUT_FILE_BED –genome REF_GENOME

Enter the bin size for merge circleRNA: 0
Enter the max blast results for each circleRNA (in case when TaxID is empty): 1

```
The output is a file called `finish.table.sort.txt` with the following columns:

* Gene ID
* Gene coordinates start
* Gene coordinates end
* Host gene for circRNAs ID
* Host gene for circRNAs
* Number of circRNAs
* Minimum size, bp
* Maximum size, bp
* Structure


| Gene ID 	| Gene coordinates start 	| Gene coordinates end 	| Host gene for circRNAs ID 	| Host gene for circRNAs 	| Number of circRNAs 	| Minimum size, bp 	| Maximum size, bp 	| Structure 	|
|-------------	|------------------------	|----------------------	|---------------------------------------	|---------------------------------------------------------------------------------------------------------------------------	|--------------------	|------------------	|------------------	|--------------------------------------------------------------------------------------------------------------------------	|
| NC_031987.2 	| 3424075 	| 3424478 	| gi\|1434974748\|ref\|XM_025904114.1\| 	| PREDICTED: Oreochromis niloticus cyclin-T2 (LOC100698210), transcript variant X6, mRNA 	| 1 	| 432 	| 432 	| intron-exon-exon-intron- 	|
| NC_031987.2 	| 33880079 	| 33880385 	| gi\|1434972051\|ref\|XM_025903599.1\| 	| PREDICTED: Oreochromis niloticus titin (LOC100702396), transcript variant X22, mRNA 	| 1 	| 337 	| 337 	| intron-exon- 	|
| NC_031987.2 	| 27189185 	| 27201070 	| gi\|1434974608\|ref\|XM_005475351.3\| 	| PREDICTED: Oreochromis niloticus CD209 antigen-like protein A (LOC102078188), mRNA 	| 1 	| 11916 	| 11916 	| intron-exon-exon-intron- 	|
| NC_031987.2 	| 22316258 	| 22317046 	| gi\|1434972940\|ref\|XM_005475492.4\| 	| PREDICTED: Oreochromis niloticus ABI family member 3 binding protein (abi3bp), transcript variant X4, mRNA 	| 1 	| 819 	| 819 	| exon-exon-intron-exon-intron- 	|
 

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

## Cite this as

Nedoluzhko A, Sharko F, Rbbani MG, Teslyuk A, Konstantinidis I, Fernandes JMO. 2020. CircParser: a novel streamlined pipeline for circular RNA structure and host gene prediction in non-model organisms. PeerJ 8:e8757 https://doi.org/10.7717/peerj.8757
