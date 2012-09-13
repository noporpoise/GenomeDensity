package DensityAround;

use strict;
use warnings;

use Carp;

use base 'Exporter';
our @EXPORT = qw(density_around load_chr_sizes
                 get_file_codes_hash get_file_codes_hash_ss
                 vcf_dump_positions dump_object_positions
                 merge_output_csvs);

my $density_cmd = "density_around";

use constant {
  RMSK_FILE => 1,
  ENSGENE_FILE_TX => 2,
  ENSGENE_FILE_CDS => 3,
  MOTIF_FILE => 4,
  HOTSPOT_FILE => 5,
  COLDSPOT_FILE => 6
  };

my %file_types = ('RMSK' => RMSK_FILE,
                  'ENSGENE_TX' => ENSGENE_FILE_TX,
                  'ENSGENE_CDS' => ENSGENE_FILE_CDS,
                  'MOTIF' => MOTIF_FILE);

# Single stranded file types
my %file_types_ss = ('HOTSPOT' => HOTSPOT_FILE,
                     'COLDSPOT' => COLDSPOT_FILE);

sub get_file_codes_hash
{
  return \%file_types;
}

sub get_file_codes_hash_ss
{
  return \%file_types_ss;
}

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
    or die("Cannot open chr sizes file '$chr_sizes_file'");

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
  my ($vcf, $out_dir, $chr_sizes, $single_stranded) = @_;

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
        
        if(!$single_stranded) {
          close($rv_handle);
        }
      }

      if($single_stranded)
      {
        my $file = $out_dir."/vcf_".$curr_chrom.".csv";
        open($fw_handle, ">$file") or die("Cannot open '$file'");
      }
      else
      {
        my $forward_file = $out_dir."/vcf_".$curr_chrom."_fw.csv";
        my $reverse_file = $out_dir."/vcf_".$curr_chrom."_rv.csv";

        open($fw_handle, ">$forward_file") or die("Cannot open '$forward_file'");
        open($rv_handle, ">$reverse_file") or die("Cannot open '$reverse_file'");
      }
    }

    my $end = $vcf_entry->{'true_POS'} + 1 + length($vcf_entry->{'true_REF'});

    print $fw_handle $vcf_entry->{'true_POS'}."\n";
    
    if(!$single_stranded)
    {
      print $rv_handle ($chr_sizes->{$curr_chrom} - $end + 1)."\n";
    }
  }

  if(defined($fw_handle))
  {
    close($fw_handle);

    if(!$single_stranded) {
      close($rv_handle);
    }
  }

  return sort {$a cmp $b} keys %chrs;
}

sub dump_hotspot_positions2
{
  my ($filetype, $handle, $out_dir, $chr_sizes, $csvsep) = @_;

  if($filetype != HOTSPOT_FILE && $filetype != COLDSPOT_FILE)
  {
    die("Not hotspot or coldspot file type! ($filetype)\n");
  }

  my $line;
  my $curr_chrom = "";
  my $obj_handle;

  my %chrs = ();

  while(defined($line = <$handle>))
  {
    my @cols = split(/\t/, $line);

    my $chr = _get_chr_name($cols[0]);
    my $start = $cols[1];
    my $end = $cols[2];

    if(!defined($chr_sizes->{$chr}))
    {
      # Ignore, warn or freak out and die
      #print STDERR "Chrom size not given for '$chr'\n";
      #die("Chrom size not given for '$chr'\n");
      next;
    }

    if($chr ne $curr_chrom)
    {
      $curr_chrom = $chr;
      $chrs{$curr_chrom} = 1;

      if(defined($obj_handle))
      {
        close($obj_handle);
      }

      my $file = $out_dir."/obj_".$curr_chrom.".csv";
      #print "opening '$file'\n";
      #exit;
      open($obj_handle, ">$file") or die("Cannot open '$file'\n");
    }

    print $obj_handle $start.$csvsep.$end."\n";
  }

  if(defined($obj_handle))
  {
    close($obj_handle);
  }

  return sort {$a cmp $b} keys %chrs;
}

sub dump_hotspot_positions
{
  my ($filetype, $handle, $out_dir, $chr_sizes, $csvsep) = @_;

  if($filetype != HOTSPOT_FILE && $filetype != COLDSPOT_FILE)
  {
    die("Not hotspot or coldspot file type! ($filetype)\n");
  }

  my $header = <$handle>;

  if(!defined($header))
  {
    die("Empty hotspot / coldspot file\n");
  }

  my @cols = split(/\t/, $header);

  my ($header_hot_start) = grep {$cols[$_] =~ /HotStart/i} 0..$#cols;
  my ($header_hot_end) = grep {$cols[$_] =~ /HotEnd/i} 0..$#cols;

  my ($header_cold_start) = grep {$cols[$_] =~ /ColdStart/i} 0..$#cols;
  my ($header_cold_end) = grep {$cols[$_] =~ /ColdEnd/i} 0..$#cols;

  my $line;
  my $curr_chrom = "";
  my $obj_handle;

  my %chrs = ();

  while(defined($line = <$handle>))
  {
    @cols = split(/\t/, $line);

    my $chr = _get_chr_name($cols[0]);

    my ($start, $end);

    if(!defined($chr_sizes->{$chr}))
    {
      # Ignore, warn or freak out and die
      #print STDERR "Chrom size not given for '$chr'\n";
      #die("Chrom size not given for '$chr'\n");
      next;
    }
    
    if($filetype == HOTSPOT_FILE)
    {
      $start = $cols[$header_hot_start];
      $end = $cols[$header_hot_end];
    }
    else
    {
      $start = $cols[$header_cold_start];
      $end = $cols[$header_cold_end];
    }

    # Hotspot / coldspots are in a silly 123.456 notation meaning 123,456
    $start =~ s/\.//g;
    $end =~ s/\.//g;

    if($chr ne $curr_chrom)
    {
      $curr_chrom = $chr;
      $chrs{$curr_chrom} = 1;

      if(defined($obj_handle))
      {
        close($obj_handle);
      }

      my $file = $out_dir."/obj_".$curr_chrom.".csv";
      #print "opening '$file'\n";
      #exit;
      open($obj_handle, ">$file") or die("Cannot open '$file'\n");
    }

    print $obj_handle $start.$csvsep.$end."\n";
  }

  if(defined($obj_handle))
  {
    close($obj_handle);
  }

  return sort {$a cmp $b} keys %chrs;
}

# Returns chromosomes
sub dump_object_positions
{
  my ($filetype, $handle, $out_dir, $chr_sizes, $csvsep) = @_;

  if(!defined($csvsep))
  {
    $csvsep = ",";
  }

  if($filetype == COLDSPOT_FILE || $filetype == HOTSPOT_FILE)
  {
    #return dump_hotspot_positions($filetype, $handle, $out_dir,
    #                              $chr_sizes, $csvsep);
    return dump_hotspot_positions2($filetype, $handle, $out_dir,
                                   $chr_sizes, $csvsep);
  }
  elsif(!(grep {$_ == $filetype} values %file_types))
  {
    croak("Invalid filetype '$filetype' - not one of " .
          "(" . join(", ", keys %file_types) . ")\n");
  }

  my $line;
  my $curr_chrom = "";
  my ($fw_handle, $rv_handle);

  my %chrs = ();

  while(defined($line = <$handle>))
  {
    my ($chr, $start, $end, $strand);

    if($filetype == RMSK_FILE)
    {
      (undef,undef,undef,undef,undef,$chr,$start,$end,undef,$strand)
        = split(/\t/, $line);
    }
    elsif($filetype == ENSGENE_FILE_TX)
    {
      (undef,undef,$chr,$strand,$start,$end) = split(/\t/, $line);
    }
    elsif($filetype == ENSGENE_FILE_CDS)
    {
      (undef,undef,$chr,$strand,undef,undef,$start,$end) = split(/\t/, $line);
    } 
    elsif($filetype == MOTIF_FILE)
    {
      ($chr,$start,$end,$strand) = split(/\t/, $line);
    }

    $chr = _get_chr_name($chr);

    if(!defined($chr_sizes->{$chr}))
    {
      # Ignore, warn or freak out and die
      #print STDERR "Chrom size not given for '$chr'\n";
      #die("Chrom size not given for '$chr'\n");
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

      my $fw_file = $out_dir."/obj_".$curr_chrom."_fw.csv";
      my $rv_file = $out_dir."/obj_".$curr_chrom."_rv.csv";

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
        chomp($line);
        print STDERR "object entry outside of chromosome size (" .
                     $chr_sizes->{$curr_chrom} . ")\n";
        print STDERR "$line\n";
        die();
      }

      print $rv_handle $rv_start.$csvsep.$rv_end."\n";
    }
    else
    {
      chomp($line);
      die("Unknown strand on line: '$line'\n");
    }
  }

  if(defined($fw_handle))
  {
    close($fw_handle);
    close($rv_handle);
  }

  return sort {$a cmp $b} keys %chrs;
}

sub merge_output_csvs
{
  my ($out_handle, $csvsep, $headings, @csv_files) = @_;

  my $bins_offsets = [];
  my $counts = [];
  my $densities = [];

  my $rows = [];

  for my $file (@csv_files)
  {
    load_similar_csvs($bins_offsets, $counts, $densities, $rows, $file);
  }

  if(defined($headings))
  {
    print $out_handle join($csvsep, @$headings, "rate.kbp")."\n";
  }

  # dump
  for(my $i = 0; $i < @$counts; $i++)
  {
    print $out_handle join($csvsep, $bins_offsets->[$i], @{$rows->[$i]},
                           $counts->[$i], $densities->[$i],
                           sprintf("%.3f", 1000*$counts->[$i] / $densities->[$i])) .
                      "\n";
  }
}

sub load_similar_csvs
{
  my ($bins_offsets, $counts, $densities, $rows, $file) = @_;

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
        # Create row
        $rows->[$line_num] = [];
      }
      elsif($bins_offsets->[$line_num] != $1)
      {
        carp("Files aren't the same! (file: $file; entry: $line_num): '$line'\n");
      }

      $counts->[$line_num] += $2;
      $densities->[$line_num] += $3;
      push(@{$rows->[$line_num]}, $2, $3);

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
