#!/usr/bin/env perl
# =============================================================================
# CD-HIT
# http://cd-hit.org/
# http://bioinformatics.burnham.org/cd-hi
#
# program written by
#                                      Weizhong Li
#                                      UCSD, San Diego Supercomputer Center
#                                      La Jolla, CA, 92093
#                                      Email liwz@sdsc.edu
# =============================================================================

use strict;
use File::Basename;
use FindBin qw($RealBin);
use Getopt::Long;
no strict "refs";

my $script_name = $0;
!system("which cd-hit > /dev/null") or die $!;

my @in;
my $out;
my $print_help;
my $para              = "";
my $host_no           = 0;
my @hosts             = ();
my $seg_no            = 64;
my $restart_in        = "";
my $queue             = 0;
my $local_cpu         = 0;
my $queue_type        = "PBS";
my $prog              = "cd-hit";
my $queue_options     = "";
my $prog_options      = "";
my %stable_files = ();

GetOptions("help|h" => \$print_help,
           "i=s{1,}" => \@in,
           "o=s" => \$out,
           "B=s" => \$para,
           "L=s" => \$local_cpu,
           "P=s" => \$prog,
           "S=s" => \$seg_no,
           "Q=s" => \$queue,
           "T=s" => \$queue_type,
           "R=s" => \$restart_in,
           "queue_options=s" => \$queue_options,
           "prog_options=s" => \$prog_options);

!$print_help || print_usage();
($#in > -1 and $out) || print_usage();

die("exits $out") if -s $out;

my $pwd            = `pwd`; chop($pwd);
my $work_dir       = "$out.tmp";
my $restart_file   = "$out.restart";
my $div            = "$work_dir/div";
my @commands       = ();
my @command_status = ();
my $command_no     = 0;
my $in = join(" ", @in);
my $cmd;
my ($i, $j, $k, $i1, $j1, $k1);

# readin a list of hosts
if ($para) {
  open(PARA, "$para") || die "can not open $para";
  while(my $ll= <PARA>){
    chop($ll); $ll =~ s/\s//g;
    next unless ($ll);
    push(@hosts, $ll); $host_no++;
  }
  close(PARA);
}

if ($queue) {
  for ($i=0; $i<$queue; $i++) {
    push(@hosts, "queue_host.$i");
  }
  $host_no = $queue;
}
if ($local_cpu) {
  for ($i=0; $i<$local_cpu; $i++) {
    push(@hosts, "localhost.$i");
  }
  $host_no = $local_cpu;
}
die "no host" unless $host_no;

`mkdir -p $work_dir`;

if (-e $restart_in) {
  read_restart();
}
else {
  assign_commands(0);
  write_restart();
}

#main runing loop
my $sleep_time = 1;
while(1) {
  #refresh job status by checking output files
  #check whether all jobs are done or not
  my $finish_flag = 1;
  my $status_change = 0;
  for ($i=0; $i<$command_no; $i++) {
    next if ($command_status[$i] eq "done");
    $finish_flag = 0;
    my $tcmd = $commands[$i];
    my $output = "";
    if ($i eq 0) {
      $output = $div;
      if (-e "$output.done") {
        $command_status[$i] = "done";
        $status_change = 1;
        chomp($seg_no = `ls $output/* | wc -l`);
        for (my $ii=0; $ii<$seg_no; $ii++) {
          `touch $output/div-$ii.done`;
        }
        assign_commands(1);
      }
      if (-e "$output.error") {
        die("error $output.sh !");
      }
    } elsif ($tcmd =~ / -o\s+(\S+)/) {
      $output = $1;
      if ((-e $output) or (-e "$output.clstr") or (-e "$output.done")) {
        $command_status[$i] = "done";
        $status_change = 1;
      }
      if (-e "$output.error") {
        die("error $output.sh !");
      }
    }
  }
  if ($status_change) {
    write_restart();
  }
  else {
    sleep($sleep_time); print ".";
  }
  last if $finish_flag;

  my $job_sent = 0;
  for ($i=0; $i<$command_no; $i++) {
    next if ($command_status[$i] eq "done");
    next if ($command_status[$i] eq "run");
    my $tcmd = $commands[$i];
    my $in1 = "";
    my $in2 = "";
    my $output = "";
    if ($tcmd =~ / -i\s+(\S+)/ ) {$in1 = $1;}
    if ($tcmd =~ / -i2\s+(\S+)/) {$in2 = $1;}
    if ($tcmd =~ / -o\s+(\S+)/ ) {$output = $1;}
    my $input_flag = 0;

    if ($i eq 0) {
      $input_flag = 1;
      $output = "$work_dir/div";
    } elsif (($in1 =~ /\S/) and ($in2 =~ /\S/)) {
      $input_flag = 1 if ((-e $in1) and (-e $in2));
    }
    elsif ($in1 =~ /\S/) {
      $input_flag = 1 if (-e $in1);
    }
    else {
      die "Error at $tcmd\n";
    }
    next unless $input_flag;

    #now input files are ready, wait
    wait_stable_file($in1) if ($in1 =~ /\S/);
    wait_stable_file($in2) if ($in2 =~ /\S/);

    my $thost_idx = wait_for_available_host();
    my $thost     = $hosts[$thost_idx];
    my $job   = basename($output);
    my $tsh   = "$work_dir/$job.sh";
    my $tlock = "$work_dir/$out.$$.$thost_idx.lock";
    open(TSH, "> $tsh") || die;
    $cmd = `date > $tlock`;
    print TSH "date > $tlock\n";
    print TSH "( $tcmd && touch $output.done ) || touch $output.error\n";
    print TSH "rm -f $tlock\n";
    close(TSH);
    if ($local_cpu) {
      $cmd = `sh $tsh  >/dev/null 2>&1 &`;
      $command_status[$i] = "run";
      print "run at $thost $tsh\n";
    }
    elsif ($queue) {
      my $t = "cd-hit-$job";
      my $dt = "$work_dir/$t";
      if ($queue_type eq "PBS") {
        open QUEUE, ">$dt.sh" or die $!;
        print QUEUE "cd $pwd\nsh $tsh\n";
        close QUEUE;
        `qsub $queue_options -N $t -o $dt.log -e $dt.err $dt.sh`;
      }
      elsif ($queue_type eq "SGE") {
        open QUEUE, ">$dt.sh" or die $!;
        print QUEUE "#!/bin/sh\n#$ -S /bin/bash\n#$ -v PATH\ncd $pwd\nsh $tsh\n";
        close QUEUE;
        `qsub $queue_options -N $t $dt.sh`;
      }
      elsif ($queue_type eq "slurm") {
        open QUEUE, ">$dt.sh" or die $!;
        print QUEUE "#!/bin/bash\ncd $pwd\nsh $tsh\n";
        close QUEUE;
        `sbatch $queue_options -J $t -e $dt.err -o $dt.log $dt.sh`;
      } 
      else {
        die "not correct queue!";
      }
      $command_status[$i] = "run";
    }
    else {
      $cmd = `ssh -xqf $thost 'cd $pwd; sh $tsh  >/dev/null 2>&1 &'`;
      $command_status[$i] = "run";
      print "run at $thost $tsh\n";
    }
    $sleep_time = 1;
    ++$job_sent;
    last if $job_sent eq $host_no;
  }

  if ((not $job_sent) and ($sleep_time < 60)) {
    $sleep_time +=5;
  }

} ############ main run loop 


######## merge all .clstr file
my $out_clstr = "$out.clstr";
if (not -s $out_clstr) {

  my @reps = ();
  for ($i=0; $i<$seg_no; $i++) {
    my $master_clstr = "$work_dir/$i-o.clstr";
#    die "No file $master_clstr\n" unless (-s $master_clstr);
    next unless (-s $master_clstr);

    my $this_rep = "$work_dir/$i-o";
    die "No rep $this_rep\n" unless (-e $this_rep);
    push(@reps, $this_rep);

    my @slave_clstr = ();
    for ($j=$i+1; $j<$seg_no; $j++) {
      my $tclstr = "$work_dir/$j.vs.$i.clstr";
      if (-s $tclstr) {push(@slave_clstr,$tclstr); }
      else {die "No file $tclstr\n";}
    }

    if (@slave_clstr) {
      my $tclstrs = join(" ", @slave_clstr);
      $cmd = "clstr_merge.pl $master_clstr $tclstrs >> $out_clstr";
      print  "$cmd\n";
      !system($cmd) or die $!;
    }
    else { #this is the last piece
      print  "cat $master_clstr >> $out_clstr";
      $cmd = `cat $master_clstr >> $out_clstr`;    
    }
  }

  my $out_clstr_ren = "$out.clstr.$$";
  open(TMP, $out_clstr) || die;
  open(OTMP, "> $out_clstr_ren") || die;
  my $no = 0;
  my $cno;
  my $ll;
  while($ll=<TMP>){
    if ($ll =~ /^>Cluster (\d+)/) {
      print OTMP ">Cluster $no\n"; $no++;
      $cno  = 0;
    }
    else {
      $ll =~ s/^\d+/$cno/;
      print OTMP $ll;
      $cno++;
    }
  }
  close(TMP);
  close(OTMP);
  sleep(10);
  $cmd = `mv $out_clstr_ren $out_clstr`;

  my $reps = join(" ", @reps);
  $cmd = `cat $reps > $out`;
}

if (1) {
  $cmd = `grep CPU $work_dir/*log`;
  my @lls = split(/\n/, $cmd);
  my $cpu = 0;
  my $ll;
  foreach $ll (@lls) {
    if ($ll =~ /CPU\s+time\s+(\d+)/) {
      $cpu += $1;
    }
  }
  print "Total CPU time: $cpu\n";
}
if (-s $out and -s $out_clstr) {
  `rm -rf $out.tmp $out.restart`;
}

sub wait_for_available_host {
  my ($i, $j, $k);
  my $sleep = 30;
  while(1) {

    for ($i=0; $i<$host_no; $i++) {
      my $thost = $hosts[$i];
      my $tlock = "$work_dir/$out.$$.$i.lock";
      next if (-e $tlock);
      return $i;
    }
    sleep($sleep);
    $sleep +=30;
    if ($sleep >= 300) { $sleep = 30; }
  }
}
########## END wait_for_available_host


sub wait_stable_file {
  my ($i, $j, $k);
  my $f = shift;
  return if ($stable_files{$f});
  return unless (-e $f);

  if (-e "$f.done") { $stable_files{$f} = 1; return; }
  my $size0 = -s $f;
  while(1) {
    sleep(10);
    my $size1 = -s $f;
    if ($size0 == $size1) { $stable_files{$f} = 1; last; }
    else {$size0 = $size1; }
  }
}
########## END wait_stable_file


sub write_restart {
  my ($i, $j, $k);
  open(RES, "> $restart_file") || die;

  for ($i=0; $i<$command_no; $i++) {
    print RES "$commands[$i]\n$command_status[$i]\n";
  }
  close(RES);
}
########## END write_restart

sub assign_commands {
  my ($i, $j, $k);
  my $cmd;
  my ($idb, $idbo, $jdb, $idbout, $idblog);
  my $type = shift;

  if ($type eq 0) {
    $command_no = 0;
    my $cpu = 1;
    if ($prog_options =~ /-T (\d+)/) {
        $cpu = $1;
    }
    $cmd = "perl $RealBin/psort2div.pl -input $in -output $div -number $seg_no -preprocessed -cpu $cpu";
    push(@commands,       $cmd);
    push(@command_status, "wait");
    $command_no++;
    return 0;
  }

  for ($i=0; $i<$seg_no; $i++) {
    $idb    = "$div/div-$i";
    $idblog = "$work_dir/$i.log";
    #compare to previous segs
    for ($j=0; $j<$i; $j++) {
      $jdb = "$work_dir/$j-o";
      $idbo = "$work_dir/$i.vs.$j";
      $cmd = "$prog-2d -i $jdb -i2 $idb -o $idbo $prog_options >> $idblog && rm -f $idb";
      push(@commands,       $cmd);
      push(@command_status, "wait");
      $command_no++;
      $idb = $idbo;
    }
    #self comparing
    $cmd = "$prog -i $idb -o $work_dir/$i-o $prog_options >> $idblog && rm -f $idb";
    push(@commands,       $cmd);
    push(@command_status, "wait");
    $command_no++;
  }
}


sub read_restart {
  $command_no = 0;
  open(RRRR, "$restart_in") || die;
  my $ll;
  while ($ll = <RRRR>) {
    chop($ll);
    push(@commands, $ll);
    $ll = <RRRR>;
    chop($ll);
    push(@command_status, $ll);
    $command_no++;
  }
  close(RRRR);
}
########## END read_restart


sub print_usage {
  print <<EOD;
Usage: $script_name options
        This script divide a big clustering job into pieces and submit
        jobs to remote computers over a network to make it parallel. 
        After all the jobs finished, the script merge the clustering
        results as if you just run a single cd-hit or cd-hit-est.

        You can also use it to divide big jobs on a single computer if
        your computer does not have enough RAM (with -L option).

Requirements:
      1 When run this script over a network, the directory where you 
        run the scripts and the input files must be available on 
        all the remote hosts with identical path.
      2 If you choose "ssh" to submit jobs, you have to have 
        passwordless ssh to any remote host, see ssh manual to
        know how to set up passwordless ssh.
      3 I suggest to use queuing system instead of ssh, 
        I currently support PBS and SGE
      4 cd-hit cd-hit-2d cd-hit-est cd-hit-est-2d
        cd-hit-div cd-hit-div.pl must be in same directory where
        this script is in.
Options 

     -i input filename in fasta format, required
     -o output filename, required
    --P program, "cd-hit" or "cd-hit-est", default "cd-hit"
    --prog_options options of program, default ""
    --B filename of list of hosts, 
        requred unless -Q or -L option is supplied
    --L number of cpus on local computer, default $local_cpu
        when you are not running it over a cluster, you can use 
        this option to divide a big clustering jobs into small
        pieces, I suggest you just use "--L 1" unless you have
        enough RAM for each cpu
    --S Number of segments to split input DB into, default $seg_no
    --Q number of jobs to submit to queue queuing system, default $queue
        by default, the program use ssh mode to submit remote jobs
    --T type of queuing system, "PBS", "SGE", "slurm" are supported, default $queue_type
    --queue_options options of queuing system, default ""
    --R restart file, used after a crash of run
     -h print this help

    Questions, bugs, contact Weizhong Li at liwz\@sdsc.edu

EOD
  exit;
}
#### END print_usage
