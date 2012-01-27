#!/usr/bin/perl

use strict;
use warnings;

use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use DensityAround;

## Prefs
my $csvsep = ",";
my $tmp_dir = "tmp_density";
my $tmp_vcf_dir = "tmp_vcf";
my $tmp_rmsk_dir = "tmp_rmsk";
my $tmp_out_dir = "tmp_rmsk_vcf";
##

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_density_around_rmsk.pl <bin_size> <num_bins> " .
               "<rmsk.txt> <chr_sizes.csv> <out.csv> [in.vcf]\n";
  print STDERR "  VCF and rmsk don't need to be sorted.  \n";
  
  exit;
}

if(@ARGV < 5 || @ARGV > 6)
{
  print_usage();
}

my $bin_size = shift;
my $num_bins = shift;
my $rmsk_file = shift;
my $chr_sizes_file = shift;
my $out_csv = shift;
my $vcf_file = shift;

if($bin_size !~ /\d+/)
{
  print_usage("Invalid bin size: +ve integers plz");
}
elsif($num_bins !~ /\d+/)
{
  print_usage("Invalid number of bins: +ve integers plz");
}

if(-e $tmp_dir)
{
  print_usage("Output directory already exists '$tmp_dir'");
}

mkdir($tmp_dir);

#
# Create output dirs
#
$tmp_vcf_dir = $tmp_dir."/".$tmp_vcf_dir;
$tmp_rmsk_dir = $tmp_dir."/".$tmp_rmsk_dir;
$tmp_out_dir = $tmp_dir."/".$tmp_out_dir;
mkdir($tmp_vcf_dir) or die("Cannot create dir '$tmp_vcf_dir'");
mkdir($tmp_rmsk_dir) or die("Cannot create dir '$tmp_rmsk_dir'");
mkdir($tmp_out_dir) or die("Cannot create dir '$tmp_out_dir'");

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

# Load chromosome sizes
my $chr_sizes = load_chr_sizes($chr_sizes_file);
my @chroms = sort {$a cmp $b} keys %$chr_sizes;

print "Chromosomes: " . join(", ", @chroms) . "\n";

#
# Split rmsk.txt file into separate chr positions
#
my $rmsk_handle;
open($rmsk_handle, $rmsk_file) or die("Cannot open rmsk file '$rmsk_file'\n");
print "Dumping rmsk coordinates..\n";
rmsk_dump_positions($rmsk_handle, $tmp_rmsk_dir, $chr_sizes, $csvsep);
close($rmsk_handle);

# Split up VCF file
print "Dumping VCF positions..\n";
vcf_dump_positions($vcf, $tmp_vcf_dir, $chr_sizes);
close($vcf_handle);

# Now run density around
my @resulting_chroms = ();

for my $chrom (@chroms)
{
  my $events_file_fw = $tmp_vcf_dir."/vcf_".$chrom."_fw.csv";
  my $objects_file_fw = $tmp_rmsk_dir."/rmsk_".$chrom."_fw.csv";
  my $counts_out_fw = $tmp_out_dir."/".$chrom."_fw.out";

  my $events_file_rv = $tmp_vcf_dir."/vcf_".$chrom."_rv.csv";
  my $objects_file_rv = $tmp_rmsk_dir."/rmsk_".$chrom."_rv.csv";
  my $counts_out_rv = $tmp_out_dir."/".$chrom."_rv.out";

  if(-e $events_file_fw && -e $objects_file_fw &&
     -e $events_file_rv && -e $objects_file_rv)
  {
    density_around($events_file_fw, $objects_file_fw, $bin_size, $num_bins,
                   $chr_sizes->{$chrom}, $counts_out_fw);

    density_around($events_file_rv, $objects_file_rv, $bin_size, $num_bins,
                   $chr_sizes->{$chrom}, $counts_out_rv);
  
    push(@resulting_chroms, $chrom);
  }

  
}

# Now merge
my @merge_files = ((map {$tmp_out_dir."/".$_."_fw.out"} @resulting_chroms),
                   (map {$tmp_out_dir."/".$_."_rv.out"} @resulting_chroms));

my $handle;

print "Merging results into file '$out_csv'...\n";

open($handle, ">$out_csv") or die("Cannot open file $out_csv\n");
merge_similar_csvs($handle, $csvsep, @merge_files);
close($handle);

print "Done!\n";
