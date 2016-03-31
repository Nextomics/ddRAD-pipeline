# ddRAD-pipeline
find_Restriction_Enzyme_cutting_site.pl -- Run to find Restriction Enzyme Cutting Site basing on Enzyme target recognition 
sequence.

DESCRIPTION

1.This script can find the enzyme target recognition sequence ,including Positive ,Reverse,Negative,Reverse-Negative;

  For example,when you input a enzyme recognition sequence like this:AAGGCT;
			The Positive (from left to right) 		:AAGGCT.
			The Reverse  (from right to left) 		:TCGGAA.
			The Negative (from left to right) 		:TTCCGA.
			The Reverse-Negative (from right to left) 	:AGCCTT.

2.The input of this script  at least includ one target recognition sequence and one genome file with fasta format,you can 
  also input more than one target recognition sequences ,but only aceppt one genome file,if you do not known fasta format,
  please visit http://en.wikipedia.org/wiki/Fasta_format for detail.Besides, you also need to input a lable to display 
  which cases you want to care,for detail please see below.

									 
3.The Result will output to screen,so you'd better redirect a output file.
  The Result form like this:	LOC_Os01g04960  1007    -       AAGGCT  4
			  	LOC_Os01g04960  3388    +       AAGGCT  1
			  	LOC_Os01g04960  3979    -       AAGGCT  3
			  	LOC_Os01g53980  205     -       AAGGCT  4
			  	LOC_Os01g53980  210     +       AAGGCT  2
			  	LOC_Os01g53980  861     -       AAGGCT  4
			  	LOC_Os01g05530  1245    +       AAGGCT  2

  The 1th column is the scaffold or chromosome.

  The 2th column is the start position of the Restriction Enzyme Cutting site ,please attention,the first base in a 
  chromosome(scaffold) is numbered 0.

  The 3th column is the strand,either '+' or '-'.

  The 4th colum is the Enzyme which can recognize this position.

  The 5th colum is use to show the label of recognition sequence,you sholuld input a lable to select which cases you 
  want to consider,no default.
	1 means The Positive Sequence.
	2 means The Reverse Sequence.
	3 means The Negative Sequence.
	4 means The Reverse-Negative Sequence.

AUTHOR

Jiang Hu : mooldhu\@gmail.com

Usage

perl find_Restriction_Enzyme_cutting_site.pl AAGGCT rice.fa 1234 > your.result.file;

perl find_Restriction_Enzyme_cutting_site.pl AAGGCT GCGGCA GGGCA rice.fa 13 > your.result.file;

perl find_Restriction_Enzyme_cutting_site.pl AAGGCT GCGGCA GGGCA ... rice.fa 123 > your.result.file;
