#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
die "Usage: perl $0   <grouping.map.list>   <all.sca.ab>   <outputDir>" unless(@ARGV==3);
my (%scaffold, %number);
open (LIST, "$ARGV[0]") || die $!;
#读取list
while (<LIST>) {
    chomp;
    my $file = $_;
    my $lg = (split /\./, basename($file))[0];
    open (IN, "$file") || die $!;
    my $num = 0;
    while (<IN>) {
        chomp;
        next if ($_ =~ /^;|^group|^\s*$/);
        $num ++;
        my @inf = split;
        my ($sca, $pos) = $inf[0] =~ /(\w+)_(\d+)/;
        $scaffold{$sca}{$lg}{$num} = join("\t", $pos, $inf[1]);
        $number{$sca}{$lg} ++;
        #print "$sca\t$lg\t$num\t$inf[1]\n";
    }
    close IN;
}
close LIST;

my (%position, %orientation);
foreach my $key (sort{$a cmp $b} keys %scaffold) {
    my (%LGpos, %LGorder);
    foreach my $key1 (sort{$a cmp $b} keys %{$scaffold{$key}}) {
        my @line = sort{$a <=> $b} keys %{$scaffold{$key}{$key1}};
        my (%posTemp, %disTemp, %mapTemp);
        my $phase = $line[0];
        for (my $i=0; $i<@line; $i++) {
            $phase = $line[$i] if ($i>0 && $line[$i] - $line[$i-1] != 1);
            $posTemp{$phase}++;
            my @tmp = split /\t/, $scaffold{$key}{$key1}{$line[$i]};
            push @{$disTemp{$phase}}, $tmp[0];
            push @{$mapTemp{$phase}}, $tmp[1];
        }
        my @order = sort {$posTemp{$b} <=> $posTemp{$a}} keys %posTemp; # achieve the phase code which have maxmum snp loci of a specific scaffold in a group
        $phase = $order[0];
        my @map = @{$mapTemp{$order[0]}};
        my @pos = @{$disTemp{$order[0]}};
        for (my $j=1; $j<@order; $j++) {
            if ($posTemp{$order[$j]} == $posTemp{$phase}) {
                my $phase2 = &phase($phase, \@map, $order[$j], \@{$mapTemp{$order[$j]}});
                if ($phase2 == $order[$j]) {
                    @map = @{$mapTemp{$order[$j]}};
                    @pos = @{$disTemp{$order[$j]}};
                    $phase = $order[$j];
                }
            }
        }
        @{$LGpos{$key1}} = @pos;
        $LGorder{$key1} = $phase;
    }
    
    my @LG = keys %LGpos;
    my $finalLG = $LG[0];
    if (scalar(@LG) > 1) {
        my $aa = $number{$key}{$LG[0]};
        for (my $i=1; $i<@LG; $i++) {
            if ($aa < $number{$key}{$LG[$i]}) {
                $aa = $number{$key}{$LG[$i]};
                $finalLG = $LG[$i];
            }
        }
    }
    $position{$finalLG}{$LGorder{$finalLG}} = $key;
    
    my ($forward, $backward) = (0, 0);
    my @snp = @{$LGpos{$finalLG}};
    if (scalar(@snp) > 1) {
        for (my $j=1; $j<@snp; $j++) {
            $forward ++ if ($snp[$j] > $snp[$j-1]);
            $backward ++ if ($snp[$j] < $snp[$j-1]);
        }
    }
    if ($forward > $backward) {
        $orientation{$key} = "+";
    }elsif ($forward < $backward) {
        $orientation{$key} = "-";
    }else {
        $orientation{$key} = "NA";
    }
}

my %geno;
open (GENO, "$ARGV[1]") || die $!;
my $first = <GENO>;
chomp $first;
while (<GENO>) {
    chomp;
    next if ($_ !~ /^scaffold/);
    my @inf = split /\s+/, $_;
    my ($sca, $pos) = $inf[0] =~ /(\w+)_(\d+)/;
    $geno{$sca}{$pos} = $_;
}
close GENO;

my $dir = $ARGV[2];
my $out2;
foreach my $key1 (sort {$a cmp $b} keys %position) {
    my $output = "$dir/$key1.geno";
    open (OUT, ">$output") || die $!;
    print OUT "$first\n";
    foreach my $key2 (sort {$a <=> $b} keys %{$position{$key1}}) {
        my $sca = $position{$key1}{$key2};
        if ($orientation{$sca} =~ /\+|NA/) {
            $out2 .= "$key1\t$sca\t+\n" if ($orientation{$sca} =~ /\+/);
            $out2 .= "$key1\t$sca\tNA\n" if ($orientation{$sca} =~ /\NA/);
            foreach my $key3 (sort{$a <=> $b} keys %{$geno{$sca}}) {
                print OUT $geno{$sca}{$key3},"\n";
            }
        }else {
            $out2 .= "$key1\t$sca\t-\n";
            foreach my $key3 (sort{$b <=> $a} keys %{$geno{$sca}}) {
                print OUT $geno{$sca}{$key3},"\n";
            }
        }
    }
    close OUT;
}

my $output2 = "$dir/scaffold_anchor.info";
open (OUT2, ">$output2") || die $!;
print OUT2 "$out2\n";
close OUT2;

sub phase {
    my ($aa1, $aa2, $aa3, $aa4) = @_;
    my $min = $aa1;
    $min = $aa3 if (&distance(@$aa2) > &distance(@$aa4));
    return $min;
}

sub distance {
    my @aa = @_;
    my $dist = $aa[0];
    @aa = sort{$a <=> $b} @aa;
    if (scalar(@aa) > 1) {
        $dist = $aa[-1] - $aa[0];
    }
    return $dist;
}
