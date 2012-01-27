package DensityAround;

use strict;
use warnings;

use Carp;

use base 'Exporter';
our @EXPORT = qw(density_around load_chr_sizes
                 vcf_dump_positions rmsk_dump_positions
                 merge_similar_csvs);

my $density_cmd = "density_around";

sub density_around
{
  my ($events_file, $objects_file, $bin_size, $num_bins,
      $chrom_size, $counts_out, $density_out) = @_;

  my $cmd = "$density_cmd " .
            "--density 1 $chrom_size " .
            "$bin_size $num_bins " .
            "\"$events_file\" \"$objects_file\" \"$counts_out\"";

  print "$cmd\n";

  `$cmd`;
}

sub _get_chr_name
{
  my ($chr) = @_;

  if($chr =~ /^chr([xy])$/i)
  {
    return 'chr'.uc($1);
  }
  else
  {
    return lc($chr);
  }
}

sub load_chr_sizes
{
  my ($chr_sizes_file) = @_;

  my $chr_sizes = {};

  open(CHROMS, $chr_sizes_file)
    or print_usage("Cannot open chr sizes file '$chr_sizes_file'");

  my $chr_line = <CHROMS>;

  if($chr_line !~ /^.*\s*[,\t]\s*\d+$/)
  {
    # Read a header line
    $chr_line = <CHROMS>;
  }

  do
  {
    if($chr_line =~ /^(.*)\s*[,\t]\s*(\d+)$/)
    {
      my ($chr, $length) = ($1,$2);
      $chr = _get_chr_name($chr);
      
      $chr_sizes->{$chr} = $length;
    }
    else
    {
      chomp($chr_line);
      die("Chromosome sizes error: '$chr_line'");
    }
  }
  while(defined($chr_line = <CHROMS>));

  close(CHROMS);
  
  return $chr_sizes;
}

# Returns chromosomes
sub vcf_dump_positions
{
  my ($vcf, $out_dir, $chr_sizes) = @_;

  my $vcf_entry;
  my $curr_chrom = "";
  my ($fw_handle, $rv_handle);

  my %chrs = ();

  while(defined($vcf_entry = $vcf->read_entry()))
  {
    my $chr = _get_chr_name($vcf_entry->{'CHROM'});

    if(!defined($chr_sizes->{$chr}))
    {
      # Alternatively, freak out and die
      next;
    }
  
    if($chr ne $curr_chrom)
    {
      # Change output file
      $curr_chrom = $chr;
      $chrs{$curr_chrom} = 1;

      if(defined($fw_handle))
      {
        close($fw_handle);
        close($rv_handle);
      }

      my $forward_file = $out_dir."/vcf_".$curr_chrom."_fw.csv";
      my $reverse_file = $out_dir."/vcf_".$curr_chrom."_rv.csv";

      open($fw_handle, ">$forward_file") or die("Cannot open '$forward_file'");
      open($rv_handle, ">$reverse_file") or die("Cannot open '$reverse_file'");
    }

    my $end = $vcf_entry->{'true_POS'} + length($vcf_entry->{'true_REF'});

    print $fw_handle $vcf_entry->{'true_POS'}."\n";
    print $rv_handle ($chr_sizes->{$curr_chrom} - $end + 1)."\n";
  }

  if(defined($fw_handle))
  {
    close($fw_handle);
    close($rv_handle);
  }

  return sort {$a cmp $b} keys %chrs;
}

# Returns chromosomes
sub rmsk_dump_positions
{
  my ($rmsk_handle, $out_dir, $chr_sizes, $csvsep) = @_;

  if(!defined($csvsep))
  {
    $csvsep = ",";
  }

  my $rmsk_line;
  my $curr_chrom = "";
  my ($fw_handle, $rv_handle);

  my %chrs = ();

  while(defined($rmsk_line = <$rmsk_handle>))
  {
    my (undef,undef,undef,undef,undef,$chr,$start,$end,undef,$strand)
      = split(/\t/, $rmsk_line);

    $chr = _get_chr_name($chr);

    if(!defined($chr_sizes->{$chr}))
    {
      # Alternatively, freak out and die
      print STDERR "Chrom size not given for '$chr'\n";
      next;
    }

    if($chr ne $curr_chrom)
    {
      $curr_chrom = $chr;
      $chrs{$curr_chrom} = 1;

      if(defined($fw_handle))
      {
        close($fw_handle);
        close($rv_handle);
      }

      my $fw_file = $out_dir."/rmsk_".$curr_chrom."_fw.csv";
      my $rv_file = $out_dir."/rmsk_".$curr_chrom."_rv.csv";

      open($fw_handle, ">$fw_file") or die("Cannot open '$fw_file'\n");
      open($rv_handle, ">$rv_file") or die("Cannot open '$rv_file'\n");
    }

    if($strand eq "+")
    {
      print $fw_handle $start.$csvsep.$end."\n";
    }
    elsif($strand eq "-")
    {
      my $rv_start = $chr_sizes->{$curr_chrom} - $end + 1;
      my $rv_end = $chr_sizes->{$curr_chrom} - $start + 1;

      if($rv_start < 0 || $rv_end < 0)
      {
        chomp($rmsk_line);
        print STDERR "rmsk entry outside of chromosome size\n";
        print STDERR "$rmsk_line\n";
        die();
      }

      print $rv_handle $rv_start.$csvsep.$rv_end."\n";
    }
    else
    {
      chomp($rmsk_line);
      die("Unknown strand on line: '$rmsk_line'\n");
    }
  }

  if(defined($fw_handle))
  {
    close($fw_handle);
    close($rv_handle);
  }

  return sort {$a cmp $b} keys %chrs;
}

sub merge_similar_csvs
{
  my ($out_handle, $csvsep, @csv_files) = @_;

  my $bins_offsets = [];
  my $counts = [];
  my $densities = [];

  for my $file (@csv_files)
  {
    load_similar_csvs($bins_offsets, $counts, $densities, $file);
  }

  # dump
  for(my $i = 0; $i < @$counts; $i++)
  {
    print $out_handle join($csvsep, $bins_offsets->[$i],
                           $counts->[$i], $densities->[$i]) . "\n";
  }
}

sub load_similar_csvs
{
  my ($bins_offsets, $counts, $densities, $file) = @_;

  open(FILE, $file) or croak("Cannot open file '$file'");

  my $line = <FILE>;

  # Check if we just read a header
  if(defined($line) && $line !~ /^\s*[0-9\.\-e\^]*\s*[,\t]/i)
  {
    # read next line
    $line = <FILE>;
  }

  my $line_num = 0;

  do
  {
    chomp($line);

    if($line =~ /^\s*(.*?)\s*[,\t]\s*(.*?)\s*[,\t]\s*(.*?)\s*$/)
    {
      if(!defined($bins_offsets->[$line_num]))
      {
        $bins_offsets->[$line_num] = $1;
      }
      elsif($bins_offsets->[$line_num] != $1)
      {
        carp("Files aren't the same! (file: $file; entry: $line_num): '$line'\n");
      }

      $counts->[$line_num] += $2;
      $densities->[$line_num] += $3;
      $line_num++;
    }
    else
    {
      carp("Unexpected line (file: $file; entry: $line_num): '$line'");
    }
  }
  while(defined($line = <FILE>));

  close(FILE);
}
