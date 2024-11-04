#!/usr/bin/perl 
use strict;
use warnings;
use Cwd 'abs_path';
die "perl $0 [input 1] [rn 1] [row 1] [input 2] [cn 2] [col 2] [rn] [output]\n" unless @ARGV == 8;

my ($input1, $rn1, $row1, $input2, $cn2, $col2, $rn, $output) = @ARGV;
my (@n, %g, @temp, $i, @index, $temp);

die "Overlap In-Output...\n" if abs_path($input1) eq abs_path($output);
die "not exists $input1\n" unless -f $input1;

if ($input1 =~ /\.gz$/) {
    open IN, "gzip -dc $input1 |" or die $!;
} else { 
    open IN, $input1 or die $!;
}
for ($i = 0; $i < $row1 - 1; ++$i) {
    <IN>;
}
chomp($_ = <IN>);
@index = split /\t/;
shift(@index) if $rn1 eq 1;
close IN;

open IN, "$input2" or die $!;
<IN> if $cn2 eq 1;
while (<IN>) {
    chomp;
    $temp = (split /\t/)[$col2 - 1];
    unless (exists $g{$temp}) {
        for ($i = 0; $i <= $#index; ++$i) {
            push(@n, $i + $rn1) if $index[$i] eq $temp;
        }
    } 
    $g{$temp} = 1;
}

if ($rn1 eq 1 and $rn eq 1) {
    unshift(@n, 0);
}

if ($input1 =~ /\.gz$/) {
    open IN, "gzip -dc $input1 |" or die $!;
} else { 
    open IN, $input1 or die $!;
}
open OU, ">$output" or die $!;
while (<IN>) {
    chomp;
    @temp = split /\t/;
    if ($#temp < $#n) {
        print OU "$_\n";
        next;
    }
    print OU "$temp[$n[0]]";
    for ($i = 1; $i <= $#n; ++$i) {
        print OU "\t$temp[$n[$i]]";
    }
    print OU "\n";
}
close OU;
close IN;
