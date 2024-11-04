#!/usr/bin/env perl

=head1 Program: functional_profile.pl

=head1 Version: V2.0

=head1 Updated: Mar 22 2020

=head1 Description: This program use to get functional profile from gene profile and functional annotation(not 1 to 1)

=head1 
            
    Usage: perl functional_profile.pl [options]

    Options: 
    -k   <str>   corr[gene\tKO1;KO2;... or gene\tKO1,KO2,... or gene KO]
    -l   <str>   corr[KO1\tgene1;gene2;... or KO1\tgene1,gene2... or KO gene]
    -pr  <str>  filename of pr, optional
    -f   <str>   gene profile file
    -b   ID of gene profile is 1 by 1
    -n   normalization or not
    -o   <str>   output dir
    gene is the same as the fist col of gene profile, only one option -k -l is permitted
    Contact: xiehailiang@aimigene.com
   
=head1 
         
=cut

use strict;
use warnings;

use Getopt::Long;

#initialize some parameters
our ($corr_k, $corr_l, $pr, $profile_gene, $obo, $nor, $outdir);  # input
our (%anno, %profile, %ss, @sample);    # anno sample profile
our (@s, @t, $i);                           # else

GetOptions( 
    "k=s" => \$corr_k,
    "l=s" => \$corr_l,
    "pr=s" => \$pr,
    "f=s" => \$profile_gene,
    "b!" => \$obo,
    "n!" => \$nor,
    "o=s" => \$outdir,
);
#get the introduction information
die `pod2text $0` if ((!$corr_k && !$corr_l) || !$profile_gene || ($corr_k && $corr_l));

$outdir ||= ".";
$outdir =~ s/\/$//;

my $pwd = $ENV{'PWD'};

$outdir = "$pwd/$outdir" if ($outdir !~ /^[\/|~]/);

$corr_k = "$pwd/$corr_k" if ($corr_k && $corr_k !~ /^[\/|~]/);

$corr_l = "$pwd/$corr_l" if ($corr_l && $corr_l !~ /^[\/|~]/);

$pr = "accumulation.pr" if !$pr && !$nor;
$pr = "normalization.pr" if !$pr && $nor;

$profile_gene = "$pwd/$profile_gene" if ($profile_gene !~ /^[\/|~]/);

$pr = $outdir.'/'.$pr;


`mkdir -p $outdir` unless (-e $outdir);

#anno
if ($corr_k) {
    die "not exists $corr_k\n" unless -f $corr_k;
    if ($corr_k =~ /\.gz$/) {
        open FA, "gzip -dc $corr_k |" or die $!;
    } else { 
        open FA, $corr_k or die $!;
    }
    while (<FA>) {
        chomp;
        @s = split /\t/;
        @t = split/[,;]/, $s[1];
        for ($i = 0; $i <= $#t; ++$i) {
            $anno{$s[0]}{$t[$i]} = 1;
        }
    }
    close FA;
}
if ($corr_l) {
    die "not exists $corr_l\n" unless -f $corr_l;
    if ($corr_l =~ /\.gz$/) {
        open FA, "gzip -dc $corr_l |" or die $!;
    } else { 
        open FA, $corr_l or die $!;
    }
    while (<FA>) {
        chomp;
        @s = split /\t/;
        @t = split/[,;]/, $s[1];
        for ($i = 0; $i <= $#t; ++$i) {
            $anno{$t[$i]}{$s[0]} = 1;
        }
    }
    close FA;
}

#gene profile
die "not exists $profile_gene\n" unless -f $profile_gene;
if ($profile_gene =~ /\.gz$/) {
    open GP, "gzip -dc $profile_gene |" or die $!;
} else {
    open GP, $profile_gene or die $!;
}
$_ = <GP>;
chomp;
@sample = split(/\t/, $_);

unless (defined $obo) {
    while (<GP>) {
        chomp;
        @s = split /\t/;
        if (exists $anno{$s[0]}) {
            for my $anno (sort keys %{$anno{$s[0]}}) {
                for ($i = 1; $i <= $#s; ++$i) {  
                    $profile{$anno}{$i} += $s[$i];  
                }
            }
        }
    }
} else {
    while (<GP>) {
        if (exists $anno{$obo}) {
            chomp;
            @s = split /\t/;
            for my $anno (sort keys %{$anno{$s[0]}}) {
                for ($i = 1; $i <= $#s; ++$i) {  
                    $profile{$anno}{$i} += $s[$i];  
                }
            }
        }
        ++$obo;
    }
}
close GP;

die("not exists hash profile!") unless keys %profile;

unless (defined $nor) {
    open PP, ">$pr" or die $!;
    print PP join("\t", @sample), "\n";

    foreach my $anno (sort keys %profile) {
        my $anno_nospace = $anno;
        $anno_nospace =~ s/\s+/\_/g;
        print PP $anno_nospace;
        for ($i = 1; $i <= $#sample; ++$i) {
            print PP "\t", sprintf "%.6g", $profile{$anno}{$i};
        }
        print PP "\n";
    }
    close PP;
} else {
    foreach my $anno (sort keys %profile) {
        for ($i = 1; $i <= $#sample; ++$i) {
            $ss{k}{$i} += $profile{$anno}{$i};
        }
    }

    open PP, ">$pr" or die $!;
    print PP join("\t", @sample), "\n";

    foreach my $anno (sort keys %profile) {
        my $anno_nospace = $anno;
        $anno_nospace =~ s/\s+/\_/g;
        print PP $anno_nospace;
        for ($i = 1; $i <= $#sample; ++$i) {
            if ($ss{k}{$i} eq 0) {
                print PP "\t0";
            } else {
                print PP "\t", sprintf "%.6g", $profile{$anno}{$i}/$ss{k}{$i};
            }
        }
        print PP "\n";
    }
    close PP;
}
