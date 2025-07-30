#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

=head1 NAME

psort and div for cdhit

=head1 SYNOPSIS

perl psort2div.pl [options]

 Options:
   -input           inputs, at least 1 input
   -minlen          min length of fa, default 150
   -preprocessed    input fa preprocessed or not, default false
   -chunk           chunk size of split file, default 500000000
   -number          number of split file, if set, takes precedence over -chunk
   -cpu             number of cpu, default 0
   -output          output dir, required
   -prefix          output prefix, default div
   -help            brief help message

=cut

my $help;
my $preprocessed;
my $number;
my $chunk = 500000000;
my $minlen = 150;
my $cpu = 0;
my @input;
my $output;
my $prefix = "div";

GetOptions('help|?' => \$help,
           'preprocessed' => \$preprocessed,
           'input=s{1,}' => \@input,
           'minlen=i' => \$minlen,
           'chunk=i' => \$chunk,
           'number=i' => \$number,
           'cpu=i' => \$cpu,
           'prefix=s' => \$prefix,
           'output=s' => \$output) or pod2usage();
pod2usage() if defined $help;
pod2usage() unless defined $output and $#input > -1;

$prefix = "$output/$prefix";
# 检查是否已有结果
if (-d $output) {
    die("exists $output and not defined number of split file!") unless defined $number;
    chomp(my $number_exists = `ls $prefix-* | wc -l`);
    if ($number_exists > 0) {
        die("number of split file not match") if $number_exists ne $number;
        for (my $i = 0; $i < $number; ++$i) {
            die("not exists $prefix-$i") unless -e "$prefix-$i";
        }
        exit 0;
    }
}
# 长度过滤及排序
`mkdir -p $output`;
my $cmd = "parallel -j $cpu \"seqkit seq -m $minlen {} | seqkit sort -r -l > $prefix-sort-{#}\" ::: ".join(" ", @input);
if ($preprocessed) {
    $cmd = "parallel -j $cpu \"ln -s \\\$(readlink -f {}) $prefix-sort-{#}\" ::: ".join(" ", @input);
}
!system($cmd) or die $!;
# 排序后文件长度
$cmd = "parallel -j $cpu \"seqkit fx2tab -l -n -i $prefix-sort-{#} | cut -f 2 > $prefix-len-{#}\" ::: ".join(" ", @input);
!system($cmd) or die $!;

# 记录不同长度的fa都来自哪些文件
my %len;
my $size = 0;
for (my $i = 1; $i <= $#input + 1; ++$i) {
    open I, "$prefix-len-$i" or die $!;
    while (<I>) {
        chomp;
        push(@{$len{$_}}, $i);
        $size += $_;
    }
}

$chunk = $size/$number if defined $number;

# 按chunk大小进行切分文件
my %file;
for (my $i = 1; $i <= $#input + 1; ++$i) {
    my $fh;
    if ($input[$i-1]=~/\.gz$/ and $preprocessed) {
        open $fh, "gunzip -c $input[$i-1] |" or die $!;
    } else {
        open $fh, "$prefix-sort-$i" or die $!;
    }
    $file{$i} = $fh;
}
$/ = "\n>";
my $out_num = 0;
my $out_chunk = 0;
open O, ">$prefix-$out_num";
for my $k1(sort {$b <=> $a} keys %len) {
    for my $k2(@{$len{$k1}}) {
        $out_chunk += $k1;
        if ($out_chunk > $chunk and $out_chunk > $k1) {
            $out_chunk = $k1;
            ++$out_num;
            open O, ">$prefix-$out_num";
        }
        my $fh = $file{$k2};
        my $seq = <$fh>;
        $seq =~ s/^>?/>/;
        $seq =~ s/>$//;
        print O "$seq";
    }
}

`rm -f $prefix-len* $prefix-sort*`;

if (defined $number) {
    if ($number < $out_num + 1) {
        my $j = $number - 1;
        for (my $i = $j + 1; $i <= $out_num; ++$i) {
            `cat $prefix-$i >> $prefix-$j; rm -f $prefix-$i`;
        }
    }
}
