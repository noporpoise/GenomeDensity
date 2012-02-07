#!/usr/bin/perl

use warnings;
use strict;

use File::Path qw(mkpath);

my $density_cmd = "~/c/density_around/scripts/vcf_density_around.pl";
my $filter_info_cmd = "~/perl/vcf_scripts/vcf_filter_by_info.pl";
my $filter_slippage_cmd = "~/perl/vcf_scripts/vcf_filter_slippage.pl";

my @bin_sizes = (100, 1000);
my @bin_numbers = (200, 200);

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./around_rmsk.pl <rmsk_dir> <max_indel> <chr_sizes.csv> <in.vcf> <out_dir>\n";
  print STDERR "  Prints commands to produce density graphs for dupe/non-dupe, INS/DEL, indel sizes\n";
  exit;
}

if(@ARGV != 5)
{
  print_usage();
}

my $rmsk_dir = shift;
my $max_indel = shift;
my $chr_sizes = shift;
my $vcf_file = shift;
my $out_dir = shift;

# get repeat classes
opendir(RMSK_DIR, $rmsk_dir) or die("Couldn't open dir '$rmsk_dir': $!");
my @rmsk_files = grep(/^rmsk\.(.*)\.txt$/i, readdir(RMSK_DIR));
closedir(RMSK_DIR);

#print "RMSK files: " . join(",", @rmsk_files) . "\n";

my @repeat_classes = @rmsk_files;

for (@repeat_classes)
{
  s/^rmsk\.(.*)\.txt$/$1/gi;
}

#print "Repeat classes: " . join(",", @repeat_classes) . "\n\n";

# Make output directory
#mkpath($out_dir) or die("Output dir couldn't be created ($out_dir)\n");

# Loop through different bin sizes / number of bins
for(my $b = 0; $b < @bin_sizes; $b++)
{
  my $bin_size = $bin_sizes[$b];
  my $num_bins = $bin_numbers[$b];

  for(my $r = 0; $r < @repeat_classes; $r++)
  {
    my $repeat_class = $repeat_classes[$r];
    my $rmsk_file = $rmsk_dir.'/'.$rmsk_files[$r];

    run_cmd($bin_size, $num_bins, $repeat_class, $rmsk_file);

    # 1=>slippage, 2=>no slippage, 3=>both
    for my $slippage (1,2,3)
    {
      # With no polarity
      run_cmd($bin_size, $num_bins, $repeat_class, $rmsk_file, $slippage);

      # With no polarity and indel size
      for(my $indel_size = 1; $indel_size <= $max_indel; $indel_size++)
      {
        run_cmd($bin_size, $num_bins, $repeat_class, $rmsk_file,
                $slippage, undef, $indel_size);
      }

      # With polarity
      for my $polarity (qw(INS DEL))
      {
        # all indel sizes
        run_cmd($bin_size, $num_bins, $repeat_class, $rmsk_file,
                $slippage, $polarity);

        for(my $indel_size = 1; $indel_size < $max_indel; $indel_size++)
        {
          run_cmd($bin_size, $num_bins, $repeat_class, $rmsk_file,
                  $slippage, $polarity, $indel_size);
        }
      }
    }
  }
}

sub run_cmd
{
  my ($bin_size, $num_bins, $repeat_class, $rmsk_file,
      $slippage, $polarity, $indel_size) = @_;

  my $out_file = $out_dir . '/density.' . $repeat_class;
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

  $run .= "$density_cmd $bin_size $num_bins RMSK $rmsk_file $chr_sizes $out_file $in";

  print "$run\n";
  #print `$run`;
}
