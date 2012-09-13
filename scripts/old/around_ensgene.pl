#!/usr/bin/perl

use warnings;
use strict;

use File::Path qw(mkpath);

my $density_cmd = "~/c/density_around/scripts/vcf_density_around.pl";
my $filter_info_cmd = "~/bioinf-perl/vcf_scripts/vcf_filter_by_info.pl";
my $filter_slippage_cmd = "~/bioinf-perl/vcf_scripts/vcf_filter_slippage.pl";

my @bin_sizes = (100, 1000);
my @bin_numbers = (200, 200);

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./around_ensgene.pl <ensGene.txt> <SNP|max_indel> <chr_sizes.csv> <in.vcf> <out_dir>\n";
  print STDERR "  Prints commands to produce density graphs for dupe/non-dupe, INS/DEL, indel sizes\n";
  exit;
}

if(@ARGV != 5)
{
  print_usage();
}

my $ens_gene_file = shift;
my $indel_or_snps = shift;
my $chr_sizes = shift;
my $vcf_file = shift;
my $out_dir = shift;

# Check input
my $max_indel;
my $is_snps = 0;

if($indel_or_snps =~ /^SNPs?$/i)
{
  $is_snps = 1;
}
else
{
  if($indel_or_snps !~ /^\d+$/)
  {
    print_usage("max_indel is not a positive integer and isn't 'SNP'");
  }

  $max_indel = $indel_or_snps;
}

# Loop through different bin sizes / number of bins
for(my $b = 0; $b < @bin_sizes; $b++)
{
  my $bin_size = $bin_sizes[$b];
  my $num_bins = $bin_numbers[$b];

  run_cmd($bin_size, $num_bins);

  # 1=>slippage, 2=>no slippage, 3=>both
  my @slippage_states = $is_snps ? (3) : (1,2,3);
  for my $slippage (@slippage_states)
  {
    # With no polarity
    run_cmd($bin_size, $num_bins, $slippage);

    if(!$is_snps)
    {
      # With no polarity and indel size
      for(my $indel_size = 1; $indel_size <= $max_indel; $indel_size++)
      {
        run_cmd($bin_size, $num_bins, $slippage, undef, $indel_size);
      }
    }

    # With polarity
    for my $polarity (qw(INS DEL))
    {
      # all indel sizes
      run_cmd($bin_size, $num_bins, $slippage, $polarity);

      if(!$is_snps)
      {
        for(my $indel_size = 1; $indel_size <= $max_indel; $indel_size++)
        {
          run_cmd($bin_size, $num_bins, $slippage, $polarity, $indel_size);
        }
      }
    }
  }
}

sub run_cmd
{
  my ($bin_size, $num_bins, $slippage, $polarity, $indel_size) = @_;

  for my $ens_gene_type (qw(ENSGENE_TX ENSGENE_CDS))
  {
    my $out_file = $out_dir . '/density.' . $ens_gene_type;
    my $run = "";
    my $in = $vcf_file;
    
    if(defined($slippage))
    {
      if($slippage == 1)
      {
        # slippage
        $out_file .= '.slippage';
        $run .= "$filter_slippage_cmd $in | ";
        $in = "-";
      }
      elsif($slippage == 2)
      {
        # non-slippage
        $out_file .= '.nonslip';
        $run .= "$filter_slippage_cmd --invert $in | ";
        $in = "-";
      }
      else
      {
        # both - do nothing
      }
    }

    if(defined($polarity) && $polarity ne "")
    {
      $out_file .= '.'.$polarity;
      $run .= "$filter_info_cmd $in INDEL '^$polarity\$' | ";
      $in = "-";
    }

    if(defined($indel_size) && $indel_size ne "")
    {
      $out_file .= '.'.$indel_size.'bp';
      $run .= "$filter_info_cmd $in SVLEN '^-?$indel_size\$' | ";
      $in = "-";
    }

    $out_file .= '.'.$bin_size.'width.'.$num_bins.'bins.csv';

    $run .= "$density_cmd $bin_size $num_bins $ens_gene_type $ens_gene_file " .
            "$chr_sizes $out_file $in";

    print "$run\n";
    #print `$run`;
  }
}

