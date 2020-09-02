#!/usr/bin/env perl

use FindBin qw($Bin);
use lib "$Bin/local/lib/perl5";

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Seq::Quality;

our $VERSION = "1.0.0";

sub main;
sub version;
sub usage;
sub dienice($);


&main();

sub main(){

  #set up variables
  my ($input_fasta, $output_fasta, $definition, $number_of_ns, $start_position, $reverse_complement, $help, $version);
  GetOptions(
      'i=s' => \$input_fasta,
      'o=s' => \$output_fasta,
      'd=s' => \$definition,
      'n=i' => \$number_of_ns,
      's=i' => \$start_position,
      'r|reverse|rc' => \$reverse_complement,
      'h|help|u|usage' => \$help,
      'v|version' => \$version);

  my $warning;

  my $in_handle;
  my $out_handle;
  my $in_io;
  my $out_io;


  if (defined($input_fasta)){
    if (-e $input_fasta && -f _ && -r _){
      open($in_handle, $input_fasta) || &dienice("Cannot read fasta input file: $input_fasta\n");
    }else{
      $help = 1;
      $warning .= "Cannot read fasta input file: $input_fasta\n";
    }
  }else{
    $warning .= "No fasta input file specified; using STDIN\n";
    $in_handle = \*STDIN;
  }

  if (defined($output_fasta)){
    if (-e $output_fasta && -f _ && -w _){
      open($out_handle, ">$output_fasta") || &dienice("Cannot write to fasta output file: $output_fasta\n");
    }elsif (-e $output_fasta){
      $help = 1;
      $warning .= "Cannot write to fasta output file: $output_fasta\n";
    }else{
      open($out_handle, ">$output_fasta") || &dienice("Cannot write to fasta output file: $output_fasta\n");
    }
  }else{
    $warning .= "No fasta output file specified; using STDOUT\n";
    $out_handle = \*STDOUT;
  }

  #warn user if options aren't valid

  if ($warning){
    warn($warning);
  }

  #deal with requests for help or version information

  if ($version){
    &version();
    exit -1;
  }

  if ($help){
    &usage();
    exit -1;
  }

  #initialize Bio::SeqIO objects

  $in_io = Bio::SeqIO->new(	-fh => $in_handle,
                            -format => "fasta");

  $out_io =  Bio::SeqIO->new(	-fh => $out_handle,
                      				-format => "fasta");

  while (my $seq = $in_io->next_seq) {
    my $len = $seq->length();
    my $sequence_name = $seq->id();
    my $padding = "";
    if (defined($start_position)){
      if ($start_position == 0){
        #generate a random position -- bioperl starts are 1-based
        $start_position = int(rand($len)) + 1;
      }elsif ($start_position < 0 || $start_position > $len){
        $warning .= "Invalid start position $start_position is not in range 1 - $len for $sequence_name\n";
        &dienice($warning);
      }
    }else{
      #generate a random position -- bioperl starts are 1-based
      $start_position = int(rand($len)) + 1;
    }

    my $def;

    if (defined($definition)){
      $def = "$definition";
    }else{
      $def = $seq->id();
    }


    if (defined($number_of_ns) && $number_of_ns > 0){
      $padding = "N" x $number_of_ns;
    }

    my $new_sequence;
    my $end_position;

    if (defined($reverse_complement) && $start_position < $len){
      $end_position = $start_position + 1;
      my $beg = $seq->trunc(1,$start_position);
      my $begrc = $beg->revcom();
      my $end = $seq->trunc($end_position,$len);
      my $endrc = $end->revcom();
      $new_sequence = Bio::Seq->new(	-id => $def,
                                      -seq => $begrc->seq() . $padding . $endrc->seq());

    }elsif (defined($reverse_complement)){
      my $full = $seq->trunc(1,$len);
      my $fullrc = $full->revcom();
      $new_sequence = Bio::Seq->new(	-id => $def,
                                      -seq => $fullrc->seq() . $padding);
    }elsif ($start_position > 1){

      $end_position = $start_position - 1;
      $new_sequence = Bio::Seq->new(	-id => $def,
                                      -seq => $seq->subseq($start_position,$len) . $padding . $seq->subseq(1,$end_position));
    }else{
      $new_sequence = Bio::Seq->new(	-id => $def,
                                      -seq => $seq->seq() . $padding);
    }
    $out_io->write_seq($new_sequence);
  }

};


sub dienice($) {
  my $warning = shift(@_);
  warn $warning;
  &usage();
  exit -1;
}

sub usage {
	&version();
	print qq/

USAGE: $0 -i input.fa -o output.fa

  Basic arguments

  -i <fasta file>          input fasta file

  -o <fasta file>          output fasta file

  Optional arguments

  -d <string>              definition line for rotated sequence

  -s <integer>             start position of rotated sequence

                           omitting this option or specifying 0 results in
                           sequences rotated to a random start position

  -n <integer>             number of Ns at breakpoint

  -r, --rc, --reverse      reverse and complement

  Other:

  -v, --version            print version number and exit

  -h, --help, -u, --usage  print this helpful screen and exit

/;
	return 1;
}

sub version {
	print qq/
rotate_fasta: rotate a fasta sequence
(Version: $VERSION)
/;
	return 1;
}
