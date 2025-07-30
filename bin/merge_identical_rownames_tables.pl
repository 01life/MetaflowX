#!/usr/bin/env perl

=head1 Program: merge_identical_rownames_tables.pl

=head1 Version: V1.0

=head1 Updated: Apr 7 2020

=head1 Description: This program use to merge identical rownames tables

=head1 
            
    Usage: perl merge_identical_rownames_tables.pl [options]

    Options: 
    -table <str>   tables list, "ID\ttable" 
    -k     <int>   which column to merge, default 2
    -o     <str>   output file
    Contact: xiehailiang@aimigene.com
   
=head1 
         
=cut

use strict;
use warnings;

use Getopt::Long;
use File::Temp;

#initialize some parameters
our ($table_list, $k, $outfile);  # input

GetOptions( 
    "table=s" => \$table_list,
    "k=i" => \$k,
    "o=s" => \$outfile,
);
#get the introduction information
die `pod2text $0` if (!$table_list || !$outfile);

$k ||= 2;

my $command = "";
open I, $table_list or die $!;
open O, ">$outfile" or die $!;
while (<I>) {
    chomp;
    my @s = split/\t/;
    print O "\t$s[0]";
    die "not exists $s[1]\n" unless -f $s[1];
    if ($s[1] =~ /\.gz$/) {
        if ($command eq "") {
            $command = "<(gzip -dc $s[1] | cut -f 1,$k) ";
        } else {
            $command .= "<(gzip -dc $s[1] | cut -f $k) ";
        }
    } else { 
        if ($command eq "") {
            $command = "<(cut -f 1,$k $s[1]) ";
        } else {
            $command .= "<(cut -f $k $s[1]) ";
        }
    }
}
print O "\n";
close I;
close O;

$command = "paste " . $command . " | sed 1d >> $outfile";
my $fh = File::Temp->new();
my $filename= $fh->filename;
print $fh "$command";
close $fh;
!system("bash $filename") or die $!;
#!system("bash", "-c", $command) or die $!;
