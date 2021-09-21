#!/usr/bin/env perl

package SHIMS3;

use FindBin qw($Bin);
use lib "$Bin/local/lib/perl5";

use strict;
use warnings;
use Getopt::Long;
use File::Which;
use Cwd;
use Cwd 'abs_path';
use Cwd 'getcwd';
use File::chmod::Recursive;
use File::Path 'rmtree';
use Bio::SeqIO;
use Bio::Seq::Quality;

sub main;
sub version;
sub usage;
sub extract_full_length;
sub identify_full_length($);
sub make_fl_consensus;
sub iterate_racon;
sub combine_all_reads($);
sub ingest_unpaired_reads($$);
sub ingest_paired_reads($$$);
sub align_and_polish($$$$);
sub align_single($);
sub align_paired($);
sub align_pacbio($);
sub align_nanopore($);
sub align_pacbio_for_consed($$);
sub align_nanopore_for_consed($$);
sub process_long_read_sam($$$);
sub revcomp($);
sub make_gap5_database($);
sub make_ace_file($$);

# enumerate package variables;
## Version information
our $VERSION = "1.0.0";

## Requirements
our $minimap2_major;
our $minimap2_minor;
our $minimap2_major_req = 2;
our $minimap2_minor_req = 17;

## Executables
our (
  $gzip_exec,
  $gunzip_exec,
  $minimap2_exec,
  $racon_exec,
  $samtools_exec,
  $tg_index_exec,
  $makeregionsfile_exec,
  $consed_exec);
## Arguments
our (
  $contigs,
  $iterations,
  $vector,
  @freads,
  @rreads,
  @sreads,
  @preads,
  @nreads,
  $output,
  $directory,
  $gap5,
  $consed,
  $help,
  $version);
## Other globals
our (
  $warning,
  $stem,
  $fl_reads,
  $longest_read,
  $longest_paf_one,
  $longest_intermediate,
  $longest_paf_two,
  $longest_consensus);

&main();


sub main{
  # Get command line options:
  GetOptions(
      'b=s' => \$vector,
      'n=s' => \@nreads,
      'h|help|u|usage' => \$help,
      'v|version' => \$version);

  # Deal with requests for help or version information

  if ($version){
    &version();
    exit -1;
  }

  if ($help){
    &usage();
    exit -1;
  }

  # Check on other arguments

  ## Is there a vector sequence file that we can read?
  if (defined($vector)){
    if (-e $vector && -f _ && -r _){
    }else{
      $help = 1;
      $warning .= "Cannot read vector sequence file: $vector\n";
    }
  }else{
    $help = 1;
    $warning .= "No file of vector sequence specified!\n";
  }

  ## Fix the arrays of reads to cope with the possibility that they can be
  ## specified independently and as a comma-separated list
  @nreads = split(/\,/,(join(",",@nreads)));

  ## Check that we have nanopore reads to work with
  unless (@nreads){
    $help = 1;
    $warning .= "No nanopore read files: use the -n flag to provide a file of reads\n";
  }

  ## Check whether read files exist and are readable
  foreach my $readfile (@freads, @rreads, @sreads, @preads, @nreads){
    unless (-e $readfile && -f _ && -r _ ){
      $help = 1;
      $warning .= "Can't access reads in: $readfile\n";
    }
  }

  ## gzip and gunzip
  $gzip_exec = which('gzip');
  unless (defined($gzip_exec) && -e $gzip_exec && -x _){
    $help = 1;
    $warning .= "Please install gzip: https://www.gnu.org/software/gzip/\n";
  }
  $gunzip_exec = "$gzip_exec -d ";

  ## minimap2
  $minimap2_exec = which('minimap2');
  if (defined($minimap2_exec) && -e $minimap2_exec && -x _){
    open (MMV, "$minimap2_exec --version |");
    while (<MMV>){
      if (m/^(\d+)\.(\d+)/){
        $minimap2_major = $1;
        $minimap2_minor = $2;
      }
    }
    close MMV;
    if (
      ($minimap2_major < $minimap2_major_req) ||
      ($minimap2_major == $minimap2_major_req &&
      $minimap2_minor < $minimap2_minor_req)){

      $help = 1;
      $warning .= "Please install minimap2 version >=$minimap2_major_req.$minimap2_minor_req: https://github.com/lh3/minimap2\n";
    }
  }else{
    $help = 1;
    $warning .= "Please install minimap2: https://github.com/lh3/minimap2\n";
  }

  # Warn if there's a problem

  if ($warning && $help){
    warn($warning);
    &usage();
    exit -1;
  }elsif ($warning){
    warn($warning);
  }

  # Do the work

  my $names;

  foreach my $nanopore_reads (@nreads){
    my %full_length_names = &identify_full_length($nanopore_reads);
    $names .= join("\n", keys(%full_length_names))."\n";
  }
  print $names;
  return 0;
}

sub version {
	print qq/
identify_full_length: list nanopore reads that start and end with vector sequence
(Version: $VERSION)
/;
	return 1;
}

sub usage {
	&version();
	print qq/

  USAGE: $0 -b vector.fa -n nanopore.fq

    Basic arguments

    -n <nanopore reads>      fasta or fastq file of nanopore reads

    -b <fasta file>          fasta file of vector backbone sequence

    Other:

    -v, --version            print version number and exit

    -h, --help, -u, --usage  print this helpful screen and exit

/;
	return 1;
}

sub dienice($) {
  my $warning = shift(@_);
  warn $warning;
  &usage();
  exit -1;
}

sub extract_full_length {
  my $fl_count;
  my $max_length = 0;
  my $max_seq;

  open (my $fl_handle, "| $gzip_exec > $fl_reads")
    || &dienice("can't write full-length read data to: $fl_reads\n");
  my $fl_io = Bio::SeqIO->new(	-fh => $fl_handle,
                                -format => "fastq");

  open (my $l_handle, ">$longest_read")
    || &dienice("can't write longest read to: $longest_read\n");
  my $l_io = Bio::SeqIO->new(	-fh => $l_handle,
                              -format => "fasta");

  foreach my $nanopore_reads (@nreads){
    my %full_length_names = &identify_full_length($nanopore_reads);
    if (scalar(keys(%full_length_names)) > 0){
      my $n_io;
      my $is_fasta = 0;
      ### detect format and setup Bio::SeqIO object
      ### detect format -- is it gzipped?
      if ($nanopore_reads =~ m/\.gz$/){
        open (my $n_handle, "$gunzip_exec -c $nanopore_reads |")
          || &dienice("can't read gzipped nanopore read data from: $nanopore_reads\n");
        ### detect format -- is it fasta?
        if ($nanopore_reads =~ m/\.(fa|fna|fasta|mfa)\.gz$/){
          $is_fasta = 1;
          $n_io = Bio::SeqIO->new(	-fh => $n_handle,
                                    -format => "fasta");
        }else{
          $n_io = Bio::SeqIO->new(	-fh => $n_handle,
                                    -format => "fastq");
        }
      }else{
        ### detect format -- is it fasta?
        if ($nanopore_reads =~ m/\.(fa|fna|fasta|mfa)$/){
          $is_fasta = 1;
          $n_io = Bio::SeqIO->new(	-file => "$nanopore_reads",
                                    -format => "fasta");
        }else{
          $n_io = Bio::SeqIO->new(	-file => "$nanopore_reads",
                                    -format => "fastq");
        }
      }
      while (my $seq = $n_io->next_seq) {
        if ($full_length_names{$seq->id()}){
          my $len = $seq->length();
          if ($is_fasta){
            my $qstring = "20 " x $len;
            $qstring =~ s/\s+$//;
            my $newseq = Bio::Seq::Quality->new(	-id => $seq->id(),
                                                  -seq => $seq->seq(),
                                                  -qual=> $qstring);
            $seq = $newseq;
          }
          $fl_io->write_seq($seq);
          $fl_count++;
          if ($max_length < $len){
            $max_seq = $seq;
            $max_length = $len;
          }
        }
      }
    }
  }
  close $fl_handle;

  if ($fl_count > 0){
    my $newseq = Bio::Seq->new(	-id => $output,
                                -seq => $max_seq->seq());
    $l_io->write_seq($newseq);
    close $l_handle;
    return 1;
  }else{
    close $l_handle;
    &dienice("There were no full-length nanopore reads.\n");
    return 0;
  }
}

sub identify_full_length($) {
  my $nanopore_reads = shift(@_);
  my %begvb;
  my %begve;
  my %endvb;
  my %endve;
  my %list;

  # read PAF file produced by minimap2
  ###
  # we're looking for read names (col 0)
  # that match vector (col 5)
  # with quality = 60 (col 11)
  # less than vector length (col 6)
  # distance from the start (col 3)
  # and the end (col 1 - col 2)
  # on the same strand (col 4)
  #
  # These hits shouldn't overlap in vector (col 7 and col 8)

  open(PAF, "$minimap2_exec -x map-ont $vector $nanopore_reads |")
    || &dienice("Alignment of $nanopore_reads to $vector with $minimap2_exec failed.\n");
  while (<PAF>){
    my @paf = split(/\t/, $_);
    if ($paf[11] == 60){
      if ($paf[3] < $paf[6]){
        $begvb{$paf[0]}{$paf[5]}{$paf[4]} = $paf[7];
        $begve{$paf[0]}{$paf[5]}{$paf[4]} = $paf[8];
      }elsif($paf[1] - $paf[2] < $paf[6]){
        $endvb{$paf[0]}{$paf[5]}{$paf[4]} = $paf[7];
        $endve{$paf[0]}{$paf[5]}{$paf[4]} = $paf[8];
      }
    }
  }
  close PAF;

  foreach my $r (keys(%begvb)){
    foreach my $v (keys(%{$begvb{$r}})){
      foreach my $s (keys(%{$begvb{$r}{$v}})){
        if (defined($endvb{$r}{$v}{$s})
          && ($begve{$r}{$v}{$s} < $endvb{$r}{$v}{$s}
          || $begvb{$r}{$v}{$s} > $endve{$r}{$v}{$s})){
          $list{$r} = 1;
        }
      }
    }
  }

  return %list;
}

sub make_fl_consensus {
  system("$minimap2_exec -x map-ont $longest_read $fl_reads > $longest_paf_one");
  system("$racon_exec $fl_reads $longest_paf_one $longest_read > $longest_intermediate");
  system("$minimap2_exec -x map-ont $longest_intermediate $fl_reads > $longest_paf_two");
  system("$racon_exec $fl_reads $longest_paf_two $longest_intermediate > $longest_consensus");
  return 1;
}

sub iterate_racon {
  my $all_reads = "$stem.all_reads.fq.gz";
  my $old_contigs = $longest_consensus;
  &combine_all_reads($all_reads);
  foreach my $i (1 .. $iterations){
    my $new_contigs = join(".", ($stem, "polish", $i, "fa"));
    my $paf = join(".", ($stem, "polish", $i, "paf"));
    &align_and_polish($old_contigs, $paf, $all_reads, $new_contigs);
    $old_contigs = $new_contigs;
  }
  return $old_contigs;
}

sub combine_all_reads($) {
  my $all_reads = shift(@_);

  #setup filehandle and Bio::SeqIO object
  open (my $all_handle, "| $gzip_exec > $all_reads") || &dienice("can't write all read data to: $all_reads\n");
  my $all_io = Bio::SeqIO->new(	-fh => $all_handle,
                                -format => "fastq");

  if (@sreads){
    foreach my $i (0 .. $#sreads){
      &ingest_unpaired_reads($all_io, $sreads[$i]);
    }
  }

  if (@freads){
    foreach my $i (0 .. $#freads){
      &ingest_paired_reads($all_io, $freads[$i], $rreads[$i]);
    }
  }

  if (@preads){
    foreach my $i (0 .. $#preads){
      &ingest_unpaired_reads($all_io, $preads[$i]);
    }
  }

  if (@nreads){
    foreach my $i (0 .. $#nreads){
      &ingest_unpaired_reads($all_io, $nreads[$i]);
    }
  }

  close $all_handle;
  return 1;
}

sub ingest_unpaired_reads($$) {
  my ($all_io, $in_name) = @_;
  my $in_io;
  my $is_fasta = 0;
  my $count = 0;
  ### detect format and setup Bio::SeqIO object
  ### detect format -- is it gzipped?
  if ($in_name =~ m/\.gz$/){
    open (my $in_handle, "$gunzip_exec -c $in_name |") || die "can't read gzipped single read data from: $in_name\n";
    ### detect format -- is it fasta?
    if ($in_name =~ m/\.(fa|fna|fasta|mfa)\.gz$/){
      $is_fasta = 1;
      $in_io = Bio::SeqIO->new(	-fh => $in_handle,
                                -format => "fasta");
    }else{
      $in_io = Bio::SeqIO->new(	-fh => $in_handle,
                                -format => "fastq");
    }
  }else{
    ### detect format -- is it fasta?
    if ($in_name =~ m/\.(fa|fna|fasta|mfa)$/){
      $is_fasta = 1;
      $in_io = Bio::SeqIO->new(	-file => "$in_name",
                                -format => "fasta");
    }else{
      $in_io = Bio::SeqIO->new(	-file => "$in_name",
                                -format => "fastq");
    }
  }
  while (my $seq = $in_io->next_seq) {
    if ($is_fasta){
      my $len = $seq->length();
      my $qstring = "20 " x $len;
      $qstring =~ s/\s+$//;
      my $newseq = Bio::Seq::Quality->new(	-id => $seq->id(),
                                            -seq => $seq->seq(),
                                            -qual=> $qstring);
      $seq = $newseq;
    }
    $all_io->write_seq($seq);
    $count++
  }
  return $count;
}

sub ingest_paired_reads($$$) {
  my ($all_io, $fin_name, $rin_name) = @_;
  my $fin_io;
  my $rin_io;
  my $is_fasta = 0;
  my $count = 0;
  ### detect format and setup Bio::SeqIO objects
  ### detect format -- is it gzipped?
  if ($fin_name =~ m/\.gz$/){
    open (my $fin_handle, "$gunzip_exec -c $fin_name |") || die "can't read gzipped paired read data from: $fin_name\n";
    open (my $rin_handle, "$gunzip_exec -c $rin_name |") || die "can't read gzipped paired read data from: $rin_name\n";
    ### detect format -- is it fasta?
    if ($fin_name =~ m/\.(fa|fna|fasta|mfa)\.gz$/){
      $is_fasta = 1;
      $fin_io = Bio::SeqIO->new(	-fh => $fin_handle,
                                  -format => "fasta");

      $rin_io = Bio::SeqIO->new(	-fh => $rin_handle,
                                  -format => "fasta");

    }else{
      $fin_io = Bio::SeqIO->new(	-fh => $fin_handle,
                                  -format => "fastq");

      $rin_io = Bio::SeqIO->new(	-fh => $rin_handle,
                                  -format => "fastq");
    }
  }else{
    ### detect format -- is it fasta?
    if ($fin_name =~ m/\.(fa|fna|fasta|mfa)$/){
      $is_fasta = 1;
      $fin_io = Bio::SeqIO->new(	-file => "$fin_name",
                                  -format => "fasta");
      $rin_io = Bio::SeqIO->new(	-file => "$rin_name",
                                  -format => "fasta");
    }else{
      $fin_io = Bio::SeqIO->new(	-file => "$fin_name",
                                  -format => "fastq");
      $rin_io = Bio::SeqIO->new(	-file => "$rin_name",
                                  -format => "fastq");
    }
  }
  while (my ($fseq, $rseq) = ($fin_io->next_seq, $rin_io->next_seq)) {
    if ($is_fasta){

      my $flen = $fseq->length();
      my $fqstring = "20 " x $flen;
      $fqstring =~ s/\s+$//;

      my $newfseq = Bio::Seq::Quality->new(	-id => $fseq->id(),
                                            -seq => $fseq->seq(),
                                            -qual=> $fqstring);
      $fseq = $newfseq;

      my $rlen = $rseq->length();
      my $rqstring = "20 " x $rlen;
      $rqstring =~ s/\s+$//;
      my $newrseq = Bio::Seq::Quality->new(	-id => $rseq->id(),
                                            -seq => $rseq->seq(),
                                            -qual=> $rqstring);
      $rseq = $newrseq;

    }

    ### check to see whether fwd and rev reads have the same id
    ### add /1 and /2 to differentiate

    if ($fseq->id() eq $rseq->id()){
      my $readid = $fseq->id();
      $fseq->id($readid."/1");
      $rseq->id($readid."/2");
    }
    $all_io->write_seq($fseq);
    $all_io->write_seq($rseq);
    $count++;
  }
  return $count;
}

sub align_and_polish($$$$) {
  my ($old_contigs, $paf, $all_reads, $new_contigs) = @_;

  if (-e $paf){
    print STDERR "removing previous alignment in: $paf\n";
    unlink $paf;
  }
  if (-e $new_contigs){
    print STDERR "removing previously polished contigs in: $new_contigs\n";
    unlink $new_contigs;
  }

  if (@sreads){
    foreach my $i (0 .. $#sreads){
      system("$minimap2_exec -x asm20 $old_contigs $sreads[$i] >> $paf");
    }
  }

  if (@freads){
    foreach my $i (0 .. $#freads){
      system("$minimap2_exec -x sr $old_contigs $freads[$i] $rreads[$i] >> $paf");
    }
  }

  if (@nreads){
    foreach my $i (0 .. $#nreads){
      system("$minimap2_exec -x map-ont $old_contigs $nreads[$i] >> $paf");
    }
  }

  if (@preads){
    foreach my $i (0 .. $#preads){
      system("$minimap2_exec -x map-pb $old_contigs $preads[$i] >> $paf");
    }
  }

  my $racon_cmd = "$racon_exec $all_reads $paf $old_contigs >$new_contigs";
  print STDERR "$racon_cmd\n";
  system($racon_cmd);
  return 1;
}

sub align_single($) {
  my $target = shift(@_);
  my @sams;
  foreach my $i (0 .. $#sreads){
    my $rginfo = '@RG\tID:S'.$i.'\tSM:S'.$i.'\tPL:ILLUMINA';
    my $samfile = $stem.".single.$i.sam";
    push (@sams, $samfile);
    system("$minimap2_exec -x asm20 -a -L --sam-hit-only -R '$rginfo' $target $sreads[$i] | $samtools_exec sort -O SAM - >$samfile");
  }
  return (@sams);
}

sub align_paired($) {
  my $target = shift(@_);
  my @sams;
  foreach my $i (0 .. $#freads){
    my $rginfo = '@RG\tID:FR'.$i.'\tSM:FR'.$i.'\tPL:ILLUMINA';
    my $samfile = $stem.".paired.$i.sam";
    push (@sams, $samfile);
    system("$minimap2_exec -x sr -a -L --sam-hit-only -R '$rginfo' $target $freads[$i] $rreads[$i] | $samtools_exec sort -O SAM - >$samfile");
  }
  return (@sams);
}

sub align_pacbio($) {
  my $target = shift(@_);
  my @sams;
  foreach my $i (0 .. $#preads){
    my $rginfo = '@RG\tID:P'.$i.'\tSM:P'.$i.'\tPL:PACBIO';
    my $samfile = $stem.".pacbio.$i.sam";
    push (@sams, $samfile);
    system("$minimap2_exec -x map-pb -a -L --sam-hit-only -R '$rginfo' $target $preads[$i] | $samtools_exec sort -O SAM - >$samfile");
  }
  return (@sams);
}

sub align_nanopore($) {
  my $target = shift(@_);
  my @sams;
  foreach my $i (0 .. $#nreads){
    my $rginfo = '@RG\tID:N'.$i.'\tSM:N'.$i.'\tPL:PACBIO';
    my $samfile = $stem.".nanopore.$i.sam";
    push (@sams, $samfile);
    system("$minimap2_exec -x map-ont -a -L --sam-hit-only -R '$rginfo' $target $nreads[$i] | $samtools_exec sort -O SAM - >$samfile");
  }
  return (@sams);
}

sub align_pacbio_for_consed($$) {
  my ($target, $subread_length) = @_;
  my @sams;
  foreach my $i (0 .. $#preads){
    my $rginfo = '@RG\tID:P'.$i.'\tSM:P'.$i.'\tPL:PACBIO';
    my $samfile = $stem.".pacbio.$i.consed.sam";
    push (@sams, $samfile);
    my ($minimap2, $samtools);
    open ($minimap2, "$minimap2_exec -x map-pb -a -L --sam-hit-only -R '$rginfo' $target $preads[$i] |") || &dienice("can't align $preads[$i] to $target with minimap2\n");
    open ($samtools, "| $samtools_exec sort -O SAM -o $samfile -") || &dienice("can't sort alignment with samtools\n");
    my $rcfq = &process_long_read_sam($minimap2,$samtools,$subread_length);
    close $minimap2;
    close $samtools;
    if ($rcfq){
      my $rcsamfile = $stem.".pacbio.rc.$i.consed.sam";
      push (@sams, $rcsamfile);
      my ($rcminimap2, $rcsamtools);
      my $rcfqfile = $stem."pacbio.$i.rc.fq";
      open (my $rcfqh, ">$rcfqfile") || die "can't write to $rcfqfile\n";
      my $rcfq_io = Bio::SeqIO->new(	-fh => $rcfqh,
                                      -format => "fastq");
      foreach my $rcfq_seq (@{$rcfq}){
        $rcfq_io->write_seq($rcfq_seq);
      }
      close $rcfqh;
      open ($rcminimap2, "$minimap2_exec -x map-pb -a -L --sam-hit-only -R '$rginfo' $target $rcfqfile |") || &dienice("can't align $rcfqfile to $target with minimap2\n");
      open ($rcsamtools, "| $samtools_exec sort -O SAM -o $rcsamfile -") || &dienice("can't sort alignment with samtools\n");
      &process_long_read_sam($rcminimap2,$rcsamtools,$subread_length);
      close $rcminimap2;
      close $rcsamtools;
    }
  }
  return (@sams);
}

sub align_nanopore_for_consed($$) {
  my ($target, $subread_length) = @_;
  my @sams;
  foreach my $i (0 .. $#nreads){
    my $rginfo = '@RG\tID:N'.$i.'\tSM:N'.$i.'\tPL:PACBIO';
    my $samfile = $stem.".nanopore.$i.consed.sam";
    push (@sams, $samfile);
    my ($minimap2, $samtools);
    open ($minimap2, "$minimap2_exec -x map-ont -a -L --sam-hit-only -R '$rginfo' $target $nreads[$i] |") || &dienice("can't align $nreads[$i] to $target with minimap2\n");
    open ($samtools, "| $samtools_exec sort -O SAM -o $samfile -") || &dienice("can't sort alignment with samtools\n");
    my $rcfq = &process_long_read_sam($minimap2,$samtools,$subread_length);
    close $minimap2;
    close $samtools;
    if ($rcfq){
      my $rcsamfile = $stem.".nanopore.rc.$i.consed.sam";
      push (@sams, $rcsamfile);
      my ($rcminimap2, $rcsamtools);
      my $rcfqfile = $stem."nanopore.$i.rc.fq";
      open (my $rcfqh, ">$rcfqfile") || die "can't write to $rcfqfile\n";
      my $rcfq_io = Bio::SeqIO->new(	-fh => $rcfqh,
                                      -format => "fastq");
      foreach my $rcfq_seq (@{$rcfq}){
        $rcfq_io->write_seq($rcfq_seq);
      }
      close $rcfqh;
      open ($rcminimap2, "$minimap2_exec -x map-ont -a -L --sam-hit-only -R '$rginfo' $target $rcfqfile |") || &dienice("can't align $rcfqfile to $target with minimap2\n");
      open ($rcsamtools, "| $samtools_exec sort -O SAM -o $rcsamfile -") || &dienice("can't sort alignment with samtools\n");
      &process_long_read_sam($rcminimap2,$rcsamtools,$subread_length);
      close $rcminimap2;
      close $rcsamtools;
    }
  }
  return (@sams);
}

sub process_long_read_sam($$$) {
  my ($inhandle,$outhandle,$subread_length) = @_;
  # define flags
  my $reverse_complement = 0x0010;
  my $read_unmapped = 0x0004;
  my $not_primary_alignment = 0x0100;
  my $supplementary_alignment = 0x0800;
  # collection of Bio::Seq::Quality objects
  my $return_fastq;
  # some counters
  my %alignment_number;
  # process sam file from inhandle
  while(!eof($inhandle)){
    defined($_ = readline $inhandle) || die "readline failed: $!\n";
    if (m/^\@/){
      # print out header lines without changes
      print $outhandle $_;
    }else{
      # process alignment lines
      chomp;
			# split the sam fields
			my ($read_base_name,
					$flag,
					$rname,
					$pos,
					$mapq,
					$cigar,
					$rnext,
					$pnext,
					$tlen,
					$seq,
					$qual,
					@info) = split(/\t/,$_);
      # make sure we have a valid alignment
      if (  defined($read_base_name) &&
            defined($flag) &&
            defined($rname) &&
            defined($pos) &&
            defined($mapq) &&
            defined($cigar) &&
            defined($rnext) &&
            defined($pnext) &&
            defined($tlen) &&
            defined($seq) &&
            defined($qual)){

        # take the length; we need it a couple times
        my $seqlen = length($seq);

        if ($flag & $not_primary_alignment || $flag & $supplementary_alignment){
          # ignore these
        }elsif ($flag & $reverse_complement){
          # flip reads around
          if ($qual eq "*"){
            $qual = "20 " x $seqlen;
            $qual =~ s/\s+$//;
          }else{
            $qual = reverse($qual);
          }
          my $rcseq = Bio::Seq::Quality->new(
            -id   => $read_base_name."_rc",
            -seq  => &revcomp($seq),
            -qual => $qual
          );
          push(@{$return_fastq}, $rcseq);
        }else{
          # keep track of how many alignments we have
          $alignment_number{$read_base_name}++;

          # split the CIGAR string into its various operations
          $cigar =~ s/(\D)/$1\t/g;
          my @cigarettes = split(/\t/,$cigar);

          # initialize a bunch of counters
          my $sub_cigar = "";
          my $alignment_length = 0;
          my $read_count = 1;
          my $sequence_length = 0;
          my $reference_length = 0;
          my $cumulative_seq = 0;
          my $cumulative_ref = 0;

          # loop over the CIGAR operations
					foreach my $i (@cigarettes){

            # get the number of bases and the operation
            $i =~ m/(\d+)(\D)/;
            my ($n,$o) = ($1,$2);

            # set up some variables
            my $permit_split = 0;
            my $left_clip = 0;
            my $right_clip = 0;

            # check whether we're:
            # 1. in a match
            # 2. in hard/soft clipped sequence

            if ($o =~ m/[M=X]/ && $n > 1){
              $permit_split = $n;
            }elsif ($o =~ m/[HS]/ && $n > 1){
              if ($cumulative_seq == 0){
                $left_clip = 1;
              }else{
                $right_clip = 1;
              }
            }


            if ($left_clip && $o eq "S"){
              # skip over clipped sequence
              $cumulative_seq += $n;
            }elsif ($right_clip){
              # ignore this clipped sequence
            }else{
              # initialize a counter
              my $m = 0;
              # loop over bases in the operation
              while ($n){
                # check whether we can split a subread
                if ($permit_split){
                  # if we hit the max subread length, its time to print
                  unless (  $alignment_length < $subread_length ||
                            $seqlen - $cumulative_seq < $subread_length){
                    # append current cigar op and reset counter
                    if ($m){
                      $sub_cigar .= $m.$o;
                    }
                    $m = 0;
                    # print subreads if they have matching bases
                    if ($sub_cigar =~ m/[M=X]/){
                      # get quality information
                      my $qual_out;
                      if ($qual eq "*"){
                        $qual_out = $qual;
                      }else{
                        $qual_out = substr($qual,$cumulative_seq,$sequence_length);
                      }
                      # print the subread
                      print $outhandle join("\t",(
                        $read_base_name."_".$alignment_number{$read_base_name}."_".$read_count,
                        $flag,
                        $rname,
                        $pos+$cumulative_ref,
                        $mapq,
                        $sub_cigar,
                        $rnext,
                        $pnext,
                        $sequence_length,
                        substr($seq,$cumulative_seq,$sequence_length),
                        $qual_out,
                        @info))."\n";
                    }
                    # clear counters for next subread
                    $sub_cigar = "";
                    $cumulative_seq += $sequence_length;
                    $cumulative_ref += $reference_length;
                    $alignment_length = 0;
                    $read_count += 1;
                    $sequence_length = 0;
                    $reference_length = 0;
                  }
                }
                # decrement read operation count
                $n--;
                # increment subread operation count
                $m++;
                # increment alignment length
                $alignment_length++;
                # check for valid CIGAR operation and
                # increment the right counters
                if ($o eq "M" || $o eq "=" || $o eq "X"){
                  # these consume both sequence and reference
                  $sequence_length++;
                  $reference_length++;
                }elsif ($o eq "I" || $o eq "S"){
                  # these only consume sequence
                  $sequence_length++;
                }elsif ($o eq "D" || $o eq "N"){
                  # these only consume reference
                  $reference_length++;
                }elsif ($o eq "H" || $o eq "P"){
                  # these consume neither
                }else{
                  die "invalid CIGAR operation: $i\n";
                }
              }
              # add this op to the subread CIGAR string
  						$sub_cigar .= $m.$o;
            }


          }

          # print out the last subread
					if ($sub_cigar =~ m/[M=X]/){
						my $qual_out;
						if ($qual eq "*"){
							$qual_out = $qual;
						}else{
							$qual_out = substr($qual,$cumulative_seq,$sequence_length);
						}
						print $outhandle join("\t",(
							$read_base_name."_".$alignment_number{$read_base_name}."_".$read_count,
							$flag,
							$rname,
							$pos+$cumulative_ref,
							$mapq,
							$sub_cigar,
							$rnext,
							$pnext,
							$sequence_length,
							substr($seq,$cumulative_seq,$sequence_length),
							$qual_out,
							@info))."\n";
          }
        }
      }
    }
  }
  return $return_fastq;
}

sub revcomp($) {
	my $ff = shift(@_);
	my $fc = $ff;
	$fc =~ tr/ACGTacgt/TGCAtgca/;
	my $rc = reverse($fc);
	return $rc;
}

sub make_gap5_database($){
  my $sams = shift(@_);
  my $db_created = 0;
  my $merged_sam = $stem.".merged.sam";
  system("$samtools_exec merge -O sam $merged_sam ".join(" ", @{$sams}));
  system("$tg_index_exec -o $stem -p -9 -s $merged_sam");
  return 1;
}

sub make_ace_file($$){
  my ($assembly, $sams) = @_;

  my $regions = $assembly;
  $regions =~ s/\.\S+/Regions.txt/;
  system("$makeregionsfile_exec $assembly");

  my $merged_bam = $stem.".merged.bam";
  my $nosupp_bam = $stem.".merged.nosupplement.bam";

  system("$samtools_exec merge -O bam $merged_bam ".join(" ", @{$sams}));
  system("$samtools_exec view -F 2048 -b -o $nosupp_bam $merged_bam");
  system("$samtools_exec index $nosupp_bam");

  my $consed_dir = "$directory/consed";

  if (-e "$consed_dir") {
    rmtree([ "$consed_dir" ]) || print "$! : for $consed_dir\n";
  }
  system("$consed_exec -bam2ace -bamFile $nosupp_bam -regionsFile $regions -dir $consed_dir");

  my $cwd = getcwd();
  chdir("$consed_dir/edit_dir/");
  system("$consed_exec -fixConsensus -ace bam2Ace.ace");
  system("$consed_exec -removeColumnsOfPads -ace bam2Ace.ace.1");
  chdir($cwd);


  return 1;
}
