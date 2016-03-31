#!/usr/bin/perl -w
use strict;

my $Help= <<USAGE;
################################################USAGE########################################################################


NAME

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


#############################################################################################################################
USAGE
die "$Help\n" unless (@ARGV>=3);


my ($Cut_site,$Fasta,$Lable,%Hash);

($Cut_site,$Fasta,$Lable)=&read_para(@ARGV);
#print "$Lable\n";
die "Please input a lable like 'perl find_Restriction_Enzyme_cutting_site.pl AAGGCT rice.fa 123 > your.result.file' use 1234;\n" if ($Lable!~/[1234]/ || $Lable=~/[^1234]/);

&read_fasta($Fasta,\%Hash);

for my $Enzyme_temp (@{$Cut_site}){
	my ($Enzyme,$Enzyme_length,$posi_strand,$posi_rev_strand);
	my ($miuns_strand,$miuns_rev_strand);
#	my (%posi_coor,%posi_rev_coor,%miuns_coor,%miuns_rev_coor);
	$Enzyme=uc ($Enzyme_temp);
	$posi_strand=$Enzyme;
	$Enzyme_length=length $Enzyme;
	$posi_rev_strand=&posi_rev_strand($Enzyme);
	$miuns_strand=&miuns_strand($Enzyme);
	$miuns_rev_strand=&posi_rev_strand($miuns_strand);
#	print "$posi_strand\t$posi_rev_strand\t$miuns_strand\t$miuns_rev_strand\n";
	for my $key (keys %Hash){
		my (%posi_coor,%posi_rev_coor,%miuns_coor,%miuns_rev_coor)=((),(),(),());
		%posi_coor=&find_coordinate($posi_strand,$key);
		%posi_rev_coor=&find_coordinate($posi_rev_strand,$key);
		%miuns_coor=&find_coordinate($miuns_strand,$key);
		%miuns_rev_coor=&find_coordinate($miuns_rev_strand,$key);
		&print_coor($key,$Enzyme,$Enzyme_length,\%posi_coor,\%posi_rev_coor,\%miuns_coor,\%miuns_rev_coor);
#		print "$key\n";
	}
}



################Sub Route############
sub print_coor {
	my ($chr,$Enzyme,$len,$posi_coor,$posi_rev_coor,$miuns_coor,$miuns_rev_coor)=@_;
	my (@coor_temp,@temp);
	for my $i (@{$posi_coor->{$chr}}){
		push @coor_temp,[$i,'+',1];
	}
	for my $i (@{$posi_rev_coor->{$chr}}){
		push @coor_temp,[$i,'+',2];
	}
	for my $i (@{$miuns_coor->{$chr}}){
		push @coor_temp,[$i,'-',3];
	}
	for my $i (@{$miuns_rev_coor->{$chr}}){
		push @coor_temp,[$i,'-',4];
	}
	@temp = sort {$a->[0] <=> $b->[0]} @coor_temp;
	for (my $i=0;$i<=$#temp;$i++){
		print "$chr\t$temp[$i][0]\t$temp[$i][1]\t$Enzyme\t$temp[$i][2]\n" if ($Lable=~m/$temp[$i][2]/);
	}
}



sub miuns_strand {
	my $sqe=shift;
	my @split=split //,$sqe;
	my $miu_sqe;
	for (my $i=0;$i<=$#split;$i++){
		if($split[$i]=~/A/){
			$miu_sqe.='T';
		}elsif($split[$i]=~/T/){
			$miu_sqe.='A';
		}elsif($split[$i]=~/G/){
			$miu_sqe.='C';
		}elsif($split[$i]=~/C/){
			$miu_sqe.='G';
		}else{
			die "Error:Illegal base in Enzyme $sqe\n";
		}
	}
	return $miu_sqe;
}




sub posi_rev_strand {
	my $sqe=shift;
	my @split=split //,$sqe;
	my $rev_sqe;
	for (my $i=$#split;$i>=0;$i--){
		$rev_sqe.=$split[$i];
	}
	return $rev_sqe;
}


sub find_coordinate {
	my $Enzyme_site=shift;
	my $chr=shift;
	my %Temp_coor;
	my $index=0;
	while ($index!=-1){
#		print "$chr\t$Enzyme_site\n";
		$index=index ($Hash{$chr}[0],$Enzyme_site,$index+1);
#		print "$chr\t$index\n";
		push @{$Temp_coor{$chr}},$index unless ($index==-1);
	}
	return %Temp_coor;
}

sub read_fasta {
	my ($infile,$hash)=@_;
	open IN,$infile || die "Cannot open $infile:$!\n";
	$/='>';<IN>;$/="\n";
	while(<IN>){
		my ($id,$seq,$SEQ,$length);
		if (/^(\S+)/){
			$id=$1;
		}else{
			die "No access number found in header line of fasta file:$infile!\n";
		}
		if ( $id=~/\|/ ) {
			die "No '|' allowed in the access number of fasta file:$infile!\n";
		}
		$/=">";
		$seq=<IN>;
		chomp $seq;
		$seq=~s/\s//g;
		$SEQ=uc ($seq);
		$length=length $SEQ;
		$hash->{$id}=[$SEQ,$length];
		$/="\n";
	}
	close IN;
}


sub read_para {
	my ($a,$b,$c);
	for (my $i=0;$i<$#_-1;$i++){
		push @{$a},$_[$i];
	}
	$b=$_[-2];
	$c=$_[-1];
	return ($a,$b,$c);
}
##################END####################
