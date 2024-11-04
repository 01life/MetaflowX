#!/usr/bin/env perl

############################################################
############################################################
use strict;
use Getopt::Std;
use vars qw($opt_a $opt_b $opt_r $opt_o $opt_h $opt_s);
use JSON;
use POSIX;

sub usage{
    print STDERR "usage: $0 <option> <value>...\n";
    print STDERR "Options\n";
    print STDERR "\t-s : id, required\n";
    print STDERR "\t-a : raw fastq info, required\n";
    print STDERR "\t-b : clean fastq info, optional\n";
    print STDERR "\t-r : rmhost info, optional\n";
    print STDERR "\t-o : output, required\n";
    print STDERR "\t-h : show the help message\n";
    exit(1);
}

getopts('a:b:r:o:s:h');

&usage if defined $opt_h;
&usage unless defined $opt_b or defined $opt_r;
&usage unless defined $opt_s && defined $opt_a && defined $opt_o;

sub format_bytes {
    my ($input, $level) = @_;
    unless ($input =~ /^\d+$/) {
       return $input; 
    }
    my @units = ('', 'K', 'M', 'G', 'T', 'P');
    my $size = log($input)/log($level);
    my $out = sprintf "%.2f", $level ** ($size - floor($size));
    return $out.$units[floor($size)];
}

sub rjson {
    my ($f) = @_;
    open I, $f;
    my $json;
    while (<I>) {
        $json .= $_;
    }
    $json = JSON->new->utf8->decode($json);
    my $bf = $json->{"summary"}->{"before_filtering"};
    my $af = $json->{"summary"}->{"after_filtering"};
    return ($bf->{"total_reads"}, $bf->{"total_bases"}, $bf->{"gc_content"}*100, $af->{"total_reads"});
}

sub rtxt {
    my ($f) = @_;
    open I, $f;
    my $head = <I>;
    chomp $head;
    my @head = split/\t/, $head;
    $_ = <I>;
    chomp;
    my @s = split/\t/;
    my %h;
    for (my $i = 0; $i <= $#s; ++$i) {
        $h{$head[$i]} = $s[$i];
    }
    return ($h{"reads"}, $h{"bases"}, $h{"GC(%)"});
}

sub rrmhost {
    my ($f) = @_;
    open I, $f;
    my %h;
    while(<I>) {
        chomp;
        my @s = split/\t/;
        $h{$s[0]} = $s[1];
    }
    return ($h{"Input Reads"}, $h{"Surviving Reads"}, $h{"Surviving Bases"});
}

my ($read1, $read2, $readm, $base1, $base2, $gc1, $gc2, $effective);

if ($opt_a =~ /json$/) {
    ($read1, $base1, $gc1, $readm) = &rjson($opt_a);
} else {
    ($read1, $base1, $gc1) = &rtxt($opt_a);
}
if (defined $opt_b) {
    if ($opt_b =~ /json$/) {
        ($read2, $base2, $gc2) = &rjson($opt_b);
    } else {
        ($read2, $base2, $gc2) = &rtxt($opt_b);
    }
} else {
    ($readm, $read2, $base2) = &rrmhost($opt_r);
}
$base1 = &format_bytes($base1, 1000);
$base2 = &format_bytes($base2, 1000);
$effective = sprintf "%.4g", $read2/$read1*100;
$gc1 = sprintf "%.4g", $gc1;
open O, ">","$opt_o" or die $!;
if (defined $readm) {
    my $readt = $read1 - $readm;
    my $readh = $readm - $read2;
    print O "ID\traw_reads\traw_bases\tclean_reads\tclean_bases\teffective_rate(%)\tGC(%)\tlq_reads\thost_reads\n";
    print O "$opt_s\t$read1\t$base1\t$read2\t$base2\t$effective\t$gc1\t$readt\t$readh\n";
} else {
    print O "ID\traw_reads\traw_bases\tclean_reads\tclean_bases\teffective_rate(%)\tGC(%)\n";
    print O "$opt_s\t$read1\t$base1\t$read2\t$base2\t$effective\t$gc1\n";
}
close O;

