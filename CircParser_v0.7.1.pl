#!/usr/bin/perl
# perl script.pl

use strict;


use Getopt::Long 'HelpMessage';



=head1 NAME

license - get license texts at the command line!

=head1 SYNOPSIS

  -b              circRNA input file (required)
  -g,--genome     reference genome file (required)
  -t,--tax        NCBI TaxID (optional)
  -a              genome annotation file, gff/gff3 file (optional)
  -c,--ciri       input circRNA from CIRI|CIRI2 in silico predictors, (default: input from CircExplorer2, find_circ, circFinder, and BED files)
  --np            prohibition for coordinates merging (optional)
  --threads       number of threads (CPUs) for BLAST search (default: 8)
  -h,--help       show this help message and exit
  --version,-v    current version
  

=head1 VERSION

0.71

=cut

my $version = '0.71';



sub usage {
   my $message = $_[0];
   if (defined $message && length $message) 
	 {
      	  $message .= "\n"
         unless $message =~ /\n$/;
   	  }

   my $command = $0;
   $command =~ s#^.*/##;

   print STDERR (
      $message,
      "usage: $command -b circRNA.bed -g genome.fa [-a gff] [-t TaxID] [<-c|-ciri>]\n\n",
      "Use $command --help to see all command-line options.\n");



   die("\n")
}




my $name_bed;
my $name_fasta;
my $name_gff;
my $TaxID = ".";
my $input_programm = 0;
my $no_merge = 0;

my $blast_threads = 8;


Getopt::Long::GetOptions(
   'help|h|?' => \&help,
   'b=s' => \$name_bed,
   'genome|g=s' => \$name_fasta,
   'tax|t=s' => \$TaxID,
   'a=s' => \$name_gff,
   'ciri|c' => \$input_programm,
   'threads' => \$blast_threads,
   'np' => \$no_merge,
   'version|v' =>   sub { Version () })

or usage ();
 
usage("The circRNA bed file must be specified.")
   unless defined $name_bed;

usage("The genome file must be specified.")
   unless defined $name_fasta;

sub help {
    HelpMessage(1);
}

sub Version {
    print $version."\n";
	exit (1);
}
 
 

 

print "Enter the bin size for merge circleRNA: ";
my $bin_size = <STDIN>;  
chomp $bin_size;  

print "Enter the max blast results for each circleRNA (in case when TaxID is empty): ";
my $blast_max = <STDIN>;  
chomp $blast_max;  

print "Enter the path to local blastDB [default: online]: ";
my $blast_path = <STDIN>;  
chomp $blast_path;  
  


#if ($blast_path ne  '')
#{
#print "Enter the Number of threads (CPUs) to use in the BLAST search [default: 8]: ";
#$blast_threads = <STDIN>;  
#chomp $blast_threads;  

#$ENV{BLASTDB}=$blast_path;
#$ENV{BLASTDB}="/mss1/export/galaxy/database/NR/Nt";
#system ('export BLASTDB=/mss1/export/galaxy/database/NR/');
#}


if ($no_merge)
{
if ($input_programm) 
		{
	
	 	system ('awk \'FNR>1 {print $2 "\t" $3 "\t" $4 "\t" $3 "\t" $4}\' '. $name_bed  .' | sort -k1,1 -k2,2n -k3,3n -  | uniq  > merge_bed.txt');
		}
		else
			{
			system ('awk \'{print $1 "\t" $2 "\t" $3 "\t" $2 "\t" $3}\' '. $name_bed  .' | sort -k1,1 -k2,2n -k3,3n -  | uniq  > merge_bed.txt');
			}
               
}
else
	{
	if ($input_programm) 
		{
	
	 	system ('awk \'FNR>1 {print $2 "\t" $3 "\t" $4}\' '. $name_bed  .' | sort -k1,1 -k2,2n -k3,3n - | bedtools merge -i -  -c 2,3 -o collapse -d '. $bin_size .' > merge_bed.txt');
		}
		else
			{
			system ('awk \'{print $1 "\t" $2 "\t" $3}\' '. $name_bed  .' | sort -k1,1 -k2,2n -k3,3n - | bedtools merge -i -  -c 2,3 -o collapse -d '. $bin_size .' > merge_bed.txt');
			}
       }

open (FF, "merge_bed.txt"); 

open(my $f_out, '>', "finish.table.txt");

#print $f_out "TaxId"."\t"."Gene ID"."\t"."Gene coordinates"."\t"."Host gene for circRNAs ID"."\t"."Host gene for circRNAs"."\t"."Number of circRNAs"."\t"."Minimum size, bp"."\t"."Maximum size, bp"."\n";
my $header = "Gene ID"."\t"."Gene coordinates start"."\t"."Gene coordinates end"."\t"."Host gene for circRNAs ID"."\t"."Host gene for circRNAs"."\t"."Number of circRNAs"."\t"."Minimum size, bp"."\t"."Maximum size, bp"."\t"."Structure"."\n";



my @fields;
 

mkdir './parse_result/';
mkdir './parse_result/fasta/';


while (my $line = <FF>) 
{
chomp $line;


my @fields = split(/\t/, $line);
my $ID = @fields[0];
my @fields_coord1 = split(/,/, @fields[3]);
my @fields_coord2 = split(/,/, @fields[4]);

my $size1 = @fields_coord1;

my $gene_coord = @fields_coord1[0]."\t".@fields_coord2[@fields_coord2-1];



	my $filename;
	for (my $i = 0; $i < $size1; $i++)
		{
		$filename = "parse_result/$ID"."_"."@fields_coord1[0]"."_"."@fields_coord2[0]".".bed";
		open(my $fh, '>>', $filename);
		print $fh $ID."\t".@fields_coord1[$i]."\t".@fields_coord2[$i]."\n";      
		close $fh;
		}
	
	 
	
	my $fasta_out = "parse_result/fasta/$ID"."_"."@fields_coord1[0]"."_"."@fields_coord2[0]".".fa";
	
	system('bedtools getfasta -fi '. $name_fasta .' -bed '. $filename .' -fo '. $fasta_out);
	
	unlink './blast_result.table';
	
	if ($blast_path ne  '') 
		{
		system('blastn -query '. $fasta_out .' -db '. $blast_path .' -perc_identity 90  -outfmt "6 qseqid sseqid pident evalue bitscore stitle staxids" -out blast_result.table -max_target_seqs 1000 -max_hsps 1 -num_threads '. $blast_threads .'');
		}
		else
			{
			 system('blastn -query '. $fasta_out .' -db nt -remote -perc_identity 90  -outfmt "6 qseqid sseqid pident evalue bitscore stitle staxids" -out blast_result.table -max_target_seqs 1000 -max_hsps 1');
			}
	
	system ('sort -nk5 -r  blast_result.table > blast_result.sort');

	my $value_seq_awk = qx(awk 'BEGIN{RS = ">" ; ORS = ""}NR==2{ min=length($2); max=length($2); next} {sum ++;} max < length($2) {max=length($2)} min > length($2) {min=length($2)} END {print min"\t" max "\t" sum}' $fasta_out);
	my @value_seq = split(/\t/, $value_seq_awk);

 	open(FB,  "blast_result.sort");

	my $line = <FB>;
	chomp $line;

	my $find_taxid;
	my @fields = split(/\t/, $line);

	my @fields_id = split(/:/, @fields[0]);

	close FB;

	
	open(FB,  "blast_result.sort");

	if ($TaxID == ".")
		{

		 while (my $line = <FB>)
		{
	        chomp $line;
		 my @fields = split(/\t/, $line);


 		 

	  		if (index(@fields[5], "contig") == -1 && 
			    index(@fields[5], "assembly") == -1 && 
			    index(@fields[5], "chromosome") == -1 &&
	                  index(@fields[5], "scaffold") == -1  &&
		           index(@fields[5], "clone") == -1  && 
			    index(@fields[5], "genome") == -1 && 
		           index(@fields[5], "linkage group") == -1 &&
			    index(@fields[5], "uncharacterized") == -1)
				{
				$find_taxid = $find_taxid + 1;

				if ($find_taxid == 1)
					{
					print $f_out @fields_id[0]."\t".$gene_coord."\t".@fields[1]."\t".@fields[5]."\t".@value_seq[2]."\t".@value_seq[0]."\t".@value_seq[1]."\t".$TaxID."\n";
					}


				elsif ($find_taxid <= $blast_max)
				{
				#print $f_out "."."\t"."."."\t"."."."\t".@fields[5]."\t"."."."\t"."."."\t"."."."\t".$TaxID."\n";
				print $f_out @fields_id[0]."\t".$gene_coord."\t".@fields[1]."\t".@fields[5]."\t".@value_seq[2]."\t".@value_seq[0]."\t".@value_seq[1]."\t".$TaxID."\n";
				}
				   else
					{	
					last;
					}
				}
		

 		}

		if (!$find_taxid)
		{
		#print $f_out @fields_id[0]."\t".$gene_coord."\t".@fields[1]."\t".@fields[5]."\t".@value_seq[2]."\t".@value_seq[0]."\t".@value_seq[1]."\t".$TaxID."\n";
		print $f_out @fields_id[0]."\t".$gene_coord."\t".@fields[1]."\t"."NOT ASSIGNED"."\t".@value_seq[2]."\t".@value_seq[0]."\t".@value_seq[1]."\t".$TaxID."\n";
	
		}


		}
		else
		 {
			 while (my $line = <FB>)
				{
	       			 chomp $line;
		
 				 my @fields = split(/\t/, $line);

	  				if ($TaxID == @fields[6])
						{
						print $f_out @fields_id[0]."\t".$gene_coord."\t".@fields[1]."\t".@fields[5]."\t".@value_seq[2]."\t".@value_seq[0]."\t".@value_seq[1]."\t".$TaxID."\n";
					 	$find_taxid = 1;
						last;
						}
				}
				if (!$find_taxid)
				{
				#print $f_out @fields_id[0]."\t".$gene_coord."\t".@fields[1]."\t".@fields[5]."\t".@value_seq[2]."\t".@value_seq[0]."\t".@value_seq[1]."\t"."."."\n";
				print $f_out @fields_id[0]."\t".$gene_coord."\t".@fields[1]."\t"."NOT ASSIGNED"."\t".@value_seq[2]."\t".@value_seq[0]."\t".@value_seq[1]."\t"."."."\n";
				}
			}

 
 
		close FB;


	

		

	}


 

#system('bedtools getfasta -fi '. $name_fasta .' -bed parse_result/trash.bed -fo parse_result/fasta/trash.fa');

close FF;
close $f_out;


if ($name_gff ne  '')
{

system ("samtools faidx  $name_fasta ");

system ('cut -f 1,2  '. $name_fasta .'.fai | sort -k1,1 > chrom.sizes');

system ("awk '\$1 ~ /^#/ {next}  {if (\$3 == \"region\") next;} {print \$1 \"\\t\" \$4 \"\\t\" \$5}' $name_gff  | sort -k1,1 -k2,2n -k3,3n | uniq > in_sorted.bed");
 
system ("awk '\$1 ~ /^#/ {next;} {if (\$3 == \"exon\") print \$1 \"\\t\" \$4 \"\\t\" \$5 \"\\t\" \"exon-\"}' $name_gff | sort -k1,1 -k2,2n -k3,3n | uniq  > exon_sorted.bed");

system ("bedtools complement -i in_sorted.bed -g chrom.sizes | awk '{if ((\$3-\$2)<2) next; print \$0 \"\\t\" \"intergenic-\";}'  > intergenic_sorted.bed");


system("cat exon_sorted.bed intergenic_sorted.bed | sort -k1,1 -k2,2n -k3,3n > exon_inter.bed");
system("bedtools complement -i exon_inter.bed -g chrom.sizes |  awk '{if ((\$3-\$2)<2) next; print \$0 \"\\t\" \"intron-\";}' > intron_sorted.bed");

system ('cat intron_sorted.bed intergenic_sorted.bed exon_sorted.bed | sort -k1,1 -k2,2n -k3,3n > feature.bed');

system ('bedtools intersect -wa -wb  -a finish.table.txt -b feature.bed | bedtools groupby -i - -g 1,2,3,4,5,6,7,8 -c 13  -o concat  > finish.txt');

system ('mv finish.txt finish.table.txt');

}

my $cmd=join(' ', 'sed', '-i', "'1i $header'" , 'finish.table.txt');
	
system $cmd;


system("(head -n 1 finish.table.txt && tail -n +2 finish.table.txt | sort -t'\t' -nk 6 -r ) > finish.table.sort.txt");
system ('rm -r ./parse_result/');
system ('rm finish.table.txt blast_result.sort blast_result.table merge_bed.txt');
system ('rm in_sorted.bed exon_sorted.bed intergenic_sorted.bed exon_inter.bed intron_sorted.bed feature.bed');
 

