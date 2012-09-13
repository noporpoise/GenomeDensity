#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;

# Duples pipes info from: http://docstore.mik.ua/orelly/perl/prog3/ch16_03.htm
use IPC::Open2;
local (*Reader, *Writer);

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_density_around_recomb.pl <bin_size> <max_bins> " .
               "<recomb.txt> [in.vcf]\n";
  exit;
}

if(@ARGV < 3 || @ARGV > 4)
{
  print_usage();
}

my $bin_size = shift;
my $max_bins = shift;
my $recomb_file = shift;
my $vcf_file = shift;

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-")
{
  open($vcf_handle, $vcf_file)
    or print_usage("Cannot open VCF file '$vcf_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($vcf_handle, "<&=STDIN") or print_usage("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a VCF file");
}

my $vcf = new VCFFile($vcf_handle);

#
# Open recomb.txt
#
my $recomb_handle;

open($recomb_handle, $recomb_file)
  or print_usage("Cannot open recomb file '$recomb_file'");

my $recomb_line;

my $vcf_entry = $vcf->read_entry();
my $recomb_entry = read_recomb_entry();

while(defined($vcf_entry) && defined($recomb_entry))
{
  while(defined($vcf_entry) && defined($recomb_entry) &&
        $vcf_entry->{'CHROM'} ne $recomb_entry->{'CHROM'})
  {
    while($vcf_entry->{'CHROM'} lt $recomb_entry->{'CHROM'} &&
          defined($vcf_entry = $vcf->read_entry())) {}

    if(!defined($vcf_entry))
    {
      last;
    }

    while($vcf_entry->{'CHROM'} gt $recomb_entry->{'CHROM'} &&
          defined($recomb_entry = read_recomb_entry())) {}
  }

  if(defined($vcf_entry) && defined($recomb_entry))
  {
    my $curr_chrom = $vcf_entry->{'CHROM'};

    # Open two-way pipe to density_around
    my $pid = open2(\*Reader, \*Writer,
                    "/Users/isaac/c/density_around/density_around --density ");

    # read and print variants
    do
    {
      print Writer $vcf_entry->{'POS'} . "," .
                   $vcf_entry->{'POS'}+length($vcf_entry->{'true_REF'}) . "\n";
    }
    while(defined($vcf_entry = $vcf->read_entry()) &&
          $vcf_entry->{'CHROM'} eq $curr_chrom);

    # Separator
    print "\n\n";

    # read and print recombination hotspots
    do
    {
      print Writer $recomb_entry->{'START'} . "," . $vcf_entry->{'END'} . "\n";
    }
    while(defined($recomb_entry = read_recomb_entry()) &&
          $recomb_entry->{'CHROM'} eq $curr_chrom);
  
    # Now read
    

    # Now close
    close(Writer);
    close(Reader);
    waitpid($pid, 0);
  }
}

close($vcf_handle);
close($recom_handle);

sub read_recomb_entry
{
  # DEV
  return {'CHROM' => '', 'POS' => -1};
}
