#!/usr/bin/perl

use strict;
use warnings;

use FindBin;
use lib $FindBin::Bin;

use File::Path qw(remove_tree); # rmdir for none-empty directories

use VCFFile;
use DensityAround;

## Prefs
my $csvsep = ",";
my $tmp_dir_base = "tmp_density";
my $tmp_vcf_dir = "tmp_vcf";
my $tmp_objects_dir = "tmp_obj";
my $tmp_out_dir = "tmp_obj_vcf";
##

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_density_around.pl <bin_size> <num_bins> " .
               "<FILE_TYPE> <objects.txt> <chr_sizes.csv> <out.csv> [in.vcf]\n";
  print STDERR "  VCF and objects file don't need to be sorted.  \n";

  exit;
}

if(@ARGV < 6 || @ARGV > 7)
{
  print_usage();
}

my $bin_size = shift;
my $num_bins = shift;
my $file_type = shift;
my $objects_file = shift;
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

my $file_codes_hash = get_file_codes_hash();
my $file_codes_hash_ss = get_file_codes_hash_ss();
$file_type = uc($file_type);
my $file_code = $file_codes_hash->{$file_type};
my $single_stranded = 0;

if(!defined($file_code))
{
  $file_code = $file_codes_hash_ss->{$file_type};
  $single_stranded = 1;
}

if(!defined($file_code))
{
  my @codes = (keys %$file_codes_hash, keys %$file_codes_hash_ss);
  print_usage("FILE_TYPE not one of (" . join(", ", @codes) . ")");
}

#
# Create tmp directory
#
my $tmp_dir = $tmp_dir_base;
my $created_tmp_dir = 0;

if(-e $tmp_dir)
{
  for(my $i = 1; ; $i++)
  {
    if(-e $tmp_dir)
    {
      $tmp_dir = $tmp_dir_base."_".$i;
    }
    else
    {
      mkdir($tmp_dir);
      $created_tmp_dir = 1;
      last;
    }
  }
}
else
{
  mkdir($tmp_dir);
  $created_tmp_dir = 1;
}

if(!$created_tmp_dir)
{
  print_usage("Many tmp directories already exist: '$tmp_dir..'");
}

#
# Create output dirs
#
$tmp_vcf_dir = $tmp_dir."/".$tmp_vcf_dir;
$tmp_objects_dir = $tmp_dir."/".$tmp_objects_dir;
$tmp_out_dir = $tmp_dir."/".$tmp_out_dir;
mkdir($tmp_vcf_dir) or die("Cannot create dir '$tmp_vcf_dir'");
mkdir($tmp_objects_dir) or die("Cannot create dir '$tmp_objects_dir'");
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

#
# Split objects file into separate chromosomes and strands
#
my $obj_handle;
open($obj_handle, $objects_file)
  or die("Cannot open objects file '$objects_file'\n");

print "Dumping object coordinates..\n";

dump_object_positions($file_code, $obj_handle, $tmp_objects_dir,
                      $chr_sizes, $csvsep);

close($obj_handle);

# Split up VCF file
print "Dumping VCF positions..\n";
vcf_dump_positions($vcf, $tmp_vcf_dir, $chr_sizes, $single_stranded);
close($vcf_handle);

# Now run density around
my @resulting_chroms = ();

if($single_stranded)
{
  for my $chrom (@chroms)
  {
    my $events_file = $tmp_vcf_dir."/vcf_".$chrom.".csv";
    my $objects_file = $tmp_objects_dir."/obj_".$chrom.".csv";
    my $counts_out = $tmp_out_dir."/".$chrom.".out";
  
    if(-e $events_file && -e $objects_file)
    {
      density_around($events_file, $objects_file, $bin_size, $num_bins,
                     $chr_sizes->{$chrom}, $counts_out);

      push(@resulting_chroms, $chrom);
    }
  }
}
else
{
  for my $chrom (@chroms)
  {
    my $events_file_fw = $tmp_vcf_dir."/vcf_".$chrom."_fw.csv";
    my $objects_file_fw = $tmp_objects_dir."/obj_".$chrom."_fw.csv";
    my $counts_out_fw = $tmp_out_dir."/".$chrom."_fw.out";

    my $events_file_rv = $tmp_vcf_dir."/vcf_".$chrom."_rv.csv";
    my $objects_file_rv = $tmp_objects_dir."/obj_".$chrom."_rv.csv";
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
}

print "Chromosomes: " . join(", ", @resulting_chroms) . "\n";

# Now merge
my @merge_files;

if($single_stranded)
{
  @merge_files = map {$tmp_out_dir."/".$_.".out"} @resulting_chroms;
}
else
{
  @merge_files = ((map {$tmp_out_dir."/".$_."_fw.out"} @resulting_chroms),
                  (map {$tmp_out_dir."/".$_."_rv.out"} @resulting_chroms));
}

my $handle;

print "Merging results into file '$out_csv'...\n";

open($handle, ">$out_csv") or die("Cannot open file $out_csv\n");
merge_similar_csvs($handle, $csvsep, @merge_files);
close($handle);

print "Removing temp directory '$tmp_dir'...\n";

remove_tree($tmp_dir) or die("Couldn't remove tmp directory");

print "Done!\n";
