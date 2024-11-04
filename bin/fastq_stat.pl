#!/usr/bin/env perl

############################################################
############################################################
use strict;
use Getopt::Std;
use FindBin qw($RealBin);
use vars qw($opt_a $opt_b $opt_c $opt_o $opt_h $opt_s $input_s $cmd);
use POSIX;

sub usage{
    print STDERR "usage: $0 <option> <value>...\n";
    print STDERR "Options\n";
    print STDERR "\t-s : id, required\n";
    print STDERR "\t-a : read1, required\n";
    print STDERR "\t-b : read2, optional\n";
    print STDERR "\t-c : single, optional\n";
    print STDERR "\t-o : output, required\n";
    print STDERR "\t-h : show the help message\n";
    exit(1);
}

getopts('a:b:c:o:s:h');

&usage if defined $opt_h;
&usage unless defined $opt_s && defined $opt_a && defined $opt_o;

my($sum_r, $sum_b, $sum_gc);

sub format_bytes {
    my ($input, $level) = @_;
    my @units = ('', 'K', 'M', 'G', 'T', 'P');
    my $size = log($input)/log($level);
    my $out = sprintf "%.2f", $level ** ($size - floor($size));
    return $out.$units[floor($size)];
}

sub fast_stat {
    my ($input) = @_;
    die "not exists $input\n" unless -f $input;
    if ($input =~ /\.gz$/) {
        open FQ, "gzip -dc $input |" or die $!;
    } else { 
        open FQ, $input or die $!;
    }
    while (<FQ>) {
        die "reads file does not look like a FASTQ file\n" unless /^@/;
        $_ = <FQ>;
        chomp;
        ++$sum_r;
        my $seq = $_;
        $sum_b += length($seq);
        my $GC = $seq =~ tr/gcGC/gcGC/;
        $sum_gc += $GC;
        $_ = <FQ>;
        die "reads file does not look like a FASTQ file\n" unless /^\+/;
        <FQ>;
    }
    return $sum_r;
}

my $check_a = &fast_stat($opt_a);
if (defined $opt_b) {
    my $check_b = &fast_stat($opt_b);
    die "seq of reads files are not equal\n" if $check_a ne $check_b - $check_a;
}
&fast_stat($opt_c) if defined $opt_c;

open O, ">","$opt_o" or die $!;
print O "sample\treads\tbases\tmean_length\tGC(%)\n";
if ($sum_b) {
    print O "$opt_s\t$sum_r\t", &format_bytes($sum_b, 1000), "\t", $sum_b/$sum_r, "\t", sprintf "%.4g\n", $sum_gc/$sum_b * 100;
} else {
    print O "$opt_s\t\t\t\t\n";
}
close O;

