#!/usr/bin/perl -w 
#written by Lingbin

##################################################################
## This script uses multiple hash of perl to genrate the contig ##
## chromosome assignment (v3.0).                                ##
##################################################################

use List::Util qw(sum);

if (scalar(@ARGV) < 3) {
    die "Usage: $0 <INPUT> <OUTPUT> <DISTANCE>\n";
}

my $input = $ARGV[0];
my $output = $ARGV[1];
my $distance = $ARGV[2];

# put file to hash table
my %hash;
open (FH,"$input") || die "can not open the file $input\n";
while(<FH>){
    chomp;
    my @line = split("\t");
    push @{$hash{$line[0]}{$line[5]}},$line[10];
}
close(FH);

# calculate total alignment length
my %hash2;
foreach my $ctg (sort keys %hash){
    foreach my $chr (sort keys %{$hash{$ctg}} ){
        my $total = sum @{$hash{$ctg}{$chr}};
        $hash2{$ctg}{$chr} = $total;
    }
}

# find the max chr
my %hash3;
foreach my $ctg (sort keys %hash2) {
    my @key =sort {$hash2{$ctg}{$b} <=> $hash2{$ctg}{$a}} keys %{$hash2{$ctg}};
    my $max = $key[0];
    $hash3{$ctg}{$max} = 1;
}

# put file into hash
my %hash4;
open (FH,"$input") || die "can not open the file $input\n";
while(<FH>){
    chomp;
    my @line = split("\t");
    push @{$hash4{$line[0]}{$line[5]}},$_;
}
close(FH);

# merge genome coordinate
my %hash5;
foreach my $ctg (sort keys %hash3){
    foreach my $chr (sort keys %{$hash3{$ctg}}){
    my ($prev_start, $prev_end);
    my @sorted_array = sort { (split /\t/, $a)[7] <=> (split /\t/, $b)[7] } @{$hash4{$ctg}{$chr}};
        foreach my $fileline ( @sorted_array ){
        my @splitline = split("\t",$fileline);
        my $start = $splitline[7];
        my $end = $splitline[8];
        if (!defined $prev_start) {
            $prev_start = $start;
            $prev_end = $end;
            }
        elsif ($start <= $prev_end + $distance) {
            $prev_end = max ($prev_end, $end);
            }
        else {
            my $pos = join ("\t", $prev_start, $prev_end);
            my $length = $prev_end - $prev_start;
            push @{$hash5{$ctg}{$chr}{$pos}},$length;
            $prev_start = $start;
            $prev_end = $end;
            }
        }
        my $pos = join ("\t", $prev_start, $prev_end);
        my $length = $prev_end - $prev_start;
        push @{$hash5{$ctg}{$chr}{$pos}},$length;
    }
}
sub max { $_[0] > $_[1] ? $_[0] : $_[1] }

# keep longest contig in output
open (OUT,">$output") || die "can not open the file $output\n";
foreach my $ctg (sort keys %hash5) {
    foreach my $chr (sort keys %{$hash5{$ctg}}){
    my @key =sort {$hash5{$ctg}{$chr}{$b} <=> $hash5{$ctg}{$chr}{$a}} keys %{$hash5{$ctg}{$chr}};
    my $longest = $key[0];
    print OUT "$chr\t$longest\t$ctg\n"
    }
}
close(OUT);

# sort output
system("bedtools sort -i $output -faidx /net/eichler/vol28/projects/hic_cohorts/nobackups/00.DataBase/References/CHM13/reference/T2T-CHM13v2.fasta.fai > $output.sort");
system("mv $output.sort $output");
