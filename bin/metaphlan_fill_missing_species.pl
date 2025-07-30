#!/usr/bin/env perl

open I, $ARGV[0] or die $!;
while(<I>) {
    print;
    unless(/^#/){
        chomp;
        @v=split/\t/;
        for ($i=1;$i<=$#v;++$i) {
            $v{$v[0]}{$i}=$v[$i];
        }
        if ($v[1]=~/\|/) {
            $taxid=1;
        }
        @t=split/\|/,$v[0];
        for ($i=0;$i<$#t;++$i) {
            $pid=join("|",@t[0..$i]);
            $sid="$pid|$t[$i+1]";
            $ps{$i}{$pid}{$sid}=1;
        }
    }
}

if ($taxid > 0) {
    for $k1(sort keys %v) {
        @t=split/\|/,$k1;
        @u=split/\|/,$v{$k1}{1};
        for ($i = 0; $i <= $#t; ++$i) {
            next if $u[$i] eq "";
            $tax{$u[$i]}{$t[$i]} = 1;
        }
    }
    for $k1(sort keys %tax) {
        @s = keys %{$tax{$k1}};
        if ($#s > 0) {
            print STDERR "Warning: $k1 has more than 1 tax_name!\n";
        }
    }
}

for $k1(sort {$b <=> $a} keys %ps) {
    for $k2(sort keys %{$ps{$k1}}) {
        unless(exists $v{$k2}){
            print "$k2";
            for $k3(sort keys %{$ps{$k1}{$k2}}) {
                for ($i=1+$taxid;$i<=$#v;++$i) {
                    $v{$k2}{$i}+=$v{$k3}{$i};
                }
                $taxid_k2 = $v{$k3}{1};
            }
            if ($taxid > 0) {
                $taxid_k2 =~s/\|[^\|]*$//;
                $v{$k2}{1} = $taxid_k2;
                print "\t$taxid_k2";
            }
            for ($i=1+$taxid;$i<=$#v;++$i) {
                print "\t$v{$k2}{$i}";
            }
            print "\n";
        }   
    }
}
