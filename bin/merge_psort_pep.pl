#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

=head1 NAME

merge sort pep for cdhit

=head1 SYNOPSIS

perl merge_psort_pep.pl [options]

Sequence ID format required: "SampleID|SequenceID"

 Options:
   -fa              cdhit output fa, required
   -input           inputs, at least 1 pep
   -preprocessed    input pep preprocessed or not, default false
   -cpu             number of cpu, default 0
   -output          output file, required
   -help            brief help message

=cut

my $help;
my $preprocessed;
my $fa;
my @input;
my $output;
my $cpu = 0;

GetOptions('help|?' => \$help,
           'fa=s' => \$fa,
           'preprocessed' => \$preprocessed,
           'input=s{1,}' => \@input,
           'cpu=i' => \$cpu,
           'output=s' => \$output) or pod2usage();
pod2usage() if defined $help;
pod2usage() unless defined $output and defined $fa and $#input > -1;

# 排序
my $cmd = "parallel -j $cpu \"seqkit sort -r -l {} > $output-sort-{#}\" ::: ".join(" ", @input);
if ($preprocessed) {
    $cmd = "parallel -j $cpu \"ln -s \\\$(readlink -f {}) $output-sort-{#}\" ::: ".join(" ", @input);
}
!system($cmd) or die $!;

my %file;
for (my $i = 1; $i <= $#input + 1; ++$i) {
    my $fh;
    my $fid;
    if ($input[$i-1]=~/\.gz$/ and $preprocessed) {
        open $fh, "gunzip -c $input[$i-1] |" or die $!;
        chomp($fid = `gunzip -c $input[$i-1] | head -n 1 | cut -d '>' -f2 | cut -d '|' -f1`);
    } else {
        open $fh, "$output-sort-$i" or die $!;
        chomp($fid = `head -n 1 $output-sort-$i | cut -d '>' -f2 | cut -d '|' -f1`);
    }
    $file{$fid} = $fh;
}

$/ = "\n>";
if ($fa=~/\.gz$/) {
    open I, "gunzip -c $fa |" or die $!;
} else {
    open I, $fa or die $!;
}
open O, ">$output" or die $!;
while (<I>) {
    my $seq_i = $_;
    my $exists = 0;
    if ($seq_i =~ /^>?(\S+)\|(\S+)/) {
        my $fid = $1;
        my $id = "$1|$2";
        my $fh = $file{$fid};
        while (<$fh>) {
            my $seq = $_;
            if ($seq =~ /^>?(\S+)/) {
                if ($1 eq $id) {
                    $seq =~ s/^>?/>/;
                    $seq =~ s/>$//;
                    print O "$seq";
                    $exists = 1;
                    last;
                }
            }
        }
        die("not exists $id in $fid!") if $exists < 1;
    }
    die('Invalid Sequence ID format!') if $exists < 1;
}

`rm -f $output-sort*`;

