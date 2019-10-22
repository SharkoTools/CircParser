# CircParser
CircParser: a novel streamlined pipeline for circular RNA structure and host gene prediction in non-model organisms

  CircParser is simple and fast pipeline that uses the outputs from the most common circular RNAs in silico prediction tools (CIRI, CIRI2, CircExplorer2, find_circ, and circFinder) to annotate circular RNAs, assigning
presumable host genes from local or public databases such as National Center for Biotechnology Information (NCBI). Also this pipeline can discriminate circular RNAs
based on their structural components (exonic, intronic, exon-intronic or intergenic) using genome annotation file.


# Requirements:
-
-
-

# Command line:
-  -b              circRNA input file (required)
-  -g,--genome     reference genome file (required)
-  -t,--tax        NCBI TaxID (optional)
-  -a              genome annotation file, gff/gff3 file (optional)
-  -c,--ciri       input circRNA from CIRI|CIRI2 in silico predictors, (default: input from CircExplorer2, find_circ, circFinder, and BED files)
-  --np            prohibition for coordinates merging (optional)
-  --threads       number of threads (CPUs) for BLAST search (default: 8)
-  -h,--help       show this help message and exit
-  --version,-v    current version
