#!/usr/bin/perl

use strict;
use warnings;

use FindBin;
use lib $FindBin::Bin;

use Cwd 'abs_path';
use File::Basename;

use File::Path qw(remove_tree); # rmdir for none-empty directories
#use File::Path qw(rmtree); # rmdir for none-empty directories

use VCFFile;

## Prefs
my $csvsep = ",";
my $density_cmd = dirname(abs_path($0))."/../density_around";
my $tmp_dir_base = 'tmp';
##

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_density_around2.pl [--lengths <in.csv>] <bin_size> <num_bins> " .
                 "<out.csv> <objects.bed> [in.vcf]\n";
  print STDERR "  Input files do not need to be sorted.  \n";

  exit;
}

if(@ARGV == 0)
{
  print_usage();
}

my $use_chrom_lengths = 0;
my %chr_lengths = ();

if($ARGV[0] =~ /^-?-l(engths)?$/i)
{
  shift;
  my $file = shift;

  open(FILE, $file) or die("Cannot read chr lengths from '$file'");

  while(defined(my $line = <FILE>))
  {
    chomp($line);
    my ($chr,$len) = split($csvsep,$line);
    $chr_lengths{$chr} = $len;
  }

  close(FILE);

  $use_chrom_lengths = 1;
}

# Get args
my ($bin_size, $num_bins, $out_csv, $objects_bed, $vcf_file, $unused) = @ARGV;

if(!defined($objects_bed))
{
  print_usage("Not enough arguments");
}
elsif(defined($unused))
{
  print_usage("Too many arguments");
}

if($bin_size !~ /\d+/)
{
  print_usage("Invalid bin size: +ve integers plz");
}
elsif($num_bins !~ /\d+/)
{
  print_usage("Invalid number of bins: +ve integers plz");
}

#
# Open VCF handle (do this first in case it fails)
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
# Create tmp directory
#
my $tmp_dir = $tmp_dir_base;
my $created_tmp_dir = 0;

# 20 attempts
for(my $i = 1; $i < 20; $i++)
{
  $tmp_dir = $tmp_dir_base."_".int(rand()*1000000);

  if(!(-e $tmp_dir))
  {
    mkdir($tmp_dir) or die("Cannot create directory '$tmp_dir'");
    $created_tmp_dir = 1;
    last;
  }
}

if(!$created_tmp_dir)
{
  print_usage("Many tmp directories already exist beginning with '$tmp_dir'");
}

#
# Create output sub-directories
#
my $tmp_vcf_dir = $tmp_dir."/vcf/";
my $tmp_objects_dir = $tmp_dir."/obj/";
my $tmp_out_dir = $tmp_dir."/out/";

mkdir($tmp_vcf_dir) or die("Cannot create dir '$tmp_vcf_dir'");
mkdir($tmp_objects_dir) or die("Cannot create dir '$tmp_objects_dir'");
mkdir($tmp_out_dir) or die("Cannot create dir '$tmp_out_dir'");

#
# 1) Create .csv files for each chrom from object bed files
#

print "Dumping object positions..\n";

my $has_reverse_strand = 0;

my $bed_line;
open(BED, $objects_bed) or die("Cannot open objects file '$objects_bed'");

my %object_chr_to_handle = ();

while(defined($bed_line = <BED>))
{
  chomp($bed_line);

  my ($bed_chrom, $bed_start, $bed_end, $bed_strand) = split("\t", $bed_line);

  if(!defined($bed_end))
  {
    die("Bed file line too short $bed_line");
  }

  my $key = $bed_chrom;

  if(!defined($bed_strand))
  {
    $bed_strand = '+';
  }

  if($bed_strand eq '+')
  {
    $key .= '_forward';
  }
  elsif($bed_strand eq '-')
  {
    $key .= '_reverse';
    ($bed_start, $bed_end) = (-$bed_end, -$bed_start);
    $has_reverse_strand = 1;
  }
  else
  {
    die("Bed file: don't understand line '$bed_line'");
  }

  if(!defined($object_chr_to_handle{$key}))
  {
    my $file = $tmp_objects_dir.$key;

    open($object_chr_to_handle{$key}, ">$file") or die("Cannot open '$file'");
  }

  my $fh = $object_chr_to_handle{$key};
  print $fh $bed_start.$csvsep.$bed_end."\n";
}

# Close handles
for my $handle (values %object_chr_to_handle)
{
  close($handle);
}

close(BED);

#
# 2) Create .csv files for each chrom from variants
#

print "Dumping VCF positions..\n";

my %vcf_chr_to_handle = ();

while(defined(my $vcf_entry = $vcf->read_entry()))
{
  my $chr = $vcf_entry->{'CHROM'};

  if(!defined($vcf_chr_to_handle{$chr.'_forward'}))
  {
    my $file_name = $chr.'_forward';
    my $file = $tmp_vcf_dir.'/'.$file_name;

    open($vcf_chr_to_handle{$file_name}, ">$file")
      or die("Cannot open file '$file'");

    if($has_reverse_strand)
    {
      $file_name = $chr.'_reverse';
      $file = $tmp_vcf_dir.'/'.$file_name;

      open($vcf_chr_to_handle{$file_name}, ">$file")
        or die("Cannot open file '$file'");
    }
  }

  my $fh_fw = $vcf_chr_to_handle{$chr.'_forward'};
  print $fh_fw $vcf_entry->{'POS'}."\n";

  if($has_reverse_strand)
  {
    my $fh_rv = $vcf_chr_to_handle{$chr.'_reverse'};
    print $fh_rv (-$vcf_entry->{'POS'})."\n";
  }
}

# Close handles
for my $handle (values %vcf_chr_to_handle)
{
  close($handle);
}

$vcf->close_vcf();

# Now run density around

print "Getting density around...\n";

for my $chrom_and_dir (keys %vcf_chr_to_handle)
{
  my $events_file = $tmp_vcf_dir.$chrom_and_dir;
  my $objects_file = $tmp_objects_dir.$chrom_and_dir;

  my $counts_out = $tmp_out_dir.$chrom_and_dir.".out";

  # Resolve chrom start/end
  my ($chrom_start, $chrom_end);

  if($use_chrom_lengths)
  {
    my ($chr,$dir) = ($chrom_and_dir =~ /^(.*)(_forward|_reverse)$/);

    if(!defined($chr_lengths{$chr}))
    {
      print STDERR "Missing chr length: $chr\n";
    }
    else
    {
      $chrom_start = 1;
      $chrom_end = $chr_lengths{$chr};

      if($dir eq "_reverse")
      {
        ($chrom_start, $chrom_end) = (-$chrom_end, -$chrom_start);
      }
    }
  }

  if(-e $events_file && -e $objects_file)
  {
    density_around($events_file, $objects_file, $counts_out,
                   $bin_size, $num_bins,
                   $chrom_start, $chrom_end);
  }
}

# Now merge
my @merge_files;

opendir(DIR, $tmp_out_dir);
my @output_files = grep {$_ =~ /\.out$/} readdir(DIR);
closedir(DIR);

if(@output_files > 0)
{
  print "Merging results into file '$out_csv'...\n";

  @output_files = map {$tmp_out_dir.$_} @output_files;

  my $handle;
  open($handle, ">$out_csv") or die("Cannot open file '$out_csv'\n");
  merge_csvs($handle, @output_files);
  close($handle);
}
else
{
  print STDERR "No results :-(\n";
  exit(1);
}

print "Removing temp directory '$tmp_dir'...\n";
remove_tree($tmp_dir) or die("Couldn't remove tmp directory");

print "Done!\n";


sub density_around
{
  my ($events_file, $objects_file, $counts_out, $bin_size, $num_bins,
      $chrom_start, $chrom_end) = @_;

  my $cmd = "$density_cmd";

  if(defined($chrom_start))
  {
    $cmd .= " --density $chrom_start $chrom_end";
  }

  $cmd .= " $bin_size $num_bins \"$events_file\" \"$objects_file\" \"$counts_out\"";

  #print "$cmd\n";

  `$cmd` or die("Cannot exec: `$cmd`");
}

sub merge_csvs
{
  my ($out_handle, $first_file, @csv_files) = @_;

  my @row_names = ();
  my @counts = ();
  my @densities = ();

  print "merge_csvs: '$first_file'\n";

  open(FILE, $first_file) or die("Cannot open file '$first_file'");

  while(defined(my $line = <FILE>))
  {
    chomp($line);
    my ($row, $value, $density) = split($csvsep, $line);
    push(@row_names, $row);
    push(@counts, $value);

    if(defined($density))
    {
      push(@densities, $density);
    }
  }

  close(FILE);

  for my $file (@csv_files)
  {
    open(FILE, $file) or die("Cannot open file '$file'");

    for(my $i = 0; defined(my $line = <FILE>); $i++)
    {
      my ($row, $value, $density) = split($csvsep, $line);
      $counts[$i] += $value;

      if(defined($density))
      {
        $densities[$i] += $density;
      }
    }

    close(FILE);
  }

  # Print merge csv header
  if($use_chrom_lengths)
  {
    # Density is the number of times we have seen a bin
    print $out_handle join($csvsep, 'bin', 'count', 'density')."\n";
  }
  else
  {
    print $out_handle join($csvsep, 'bin', 'count')."\n";
  }

  for(my $i = 0; defined($counts[$i]); $i++)
  {
    print $out_handle $row_names[$i] . $csvsep . $counts[$i] .
                      (defined($densities[$i]) ? $csvsep.$densities[$i] : '') .
                      "\n";
  }
}
