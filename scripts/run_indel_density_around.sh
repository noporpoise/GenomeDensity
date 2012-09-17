vcf_file=$1
chrom_sizes_file=$2
bed_file=$3
num_of_bins=$4
bin_size=$5
out_csv_base=$6

echo "run_indel_density_around.sh"
echo "  VCF file: $vcf_file"
echo "  Chromosome lengths file: $chrom_sizes_file"
echo "  BED file: $bed_file"
echo "  Number of bins: $num_of_bins"
echo "  Bin size: $bin_size"
echo "  Out CSV base: $out_csv_base"

BIONINF_PATH=~/bioinf-perl
DENSITY_PATH=~/c/density_around

# SNPs
csv_file=$out_csv_base"."$num_of_bins"bins."$bin_size"width.snps.csv"

$BIONINF_PATH/vcf_scripts/vcf_filter_variants.pl SNP $vcf_file | \
$DENSITY_PATH/scripts/vcf_density_around.pl --lengths $chrom_sizes_file $num_of_bins $bin_size $csv_file $bed_file || { echo 'Error' ; exit 1; }

# Insertions
for i in {1..10}
do
  csv_file=$out_csv_base"."$num_of_bins"bins."$bin_size"width."$i"ins.csv"

  $BIONINF_PATH/vcf_scripts/vcf_filter_ins_del.pl +$i $vcf_file | \
  $DENSITY_PATH/scripts/vcf_density_around.pl --lengths $chrom_sizes_file $num_of_bins $bin_size $csv_file $bed_file || { echo 'Error' ; exit 1; }
done

# Deletions
for i in {1..10}
do
  csv_file=$out_csv_base"."$num_of_bins"bins."$bin_size"width."$i"del.csv"

  $BIONINF_PATH/vcf_scripts/vcf_filter_ins_del.pl -$i $vcf_file | \
  $DENSITY_PATH/scripts/vcf_density_around.pl --lengths $chrom_sizes_file $num_of_bins $bin_size $csv_file $bed_file || { echo 'Error' ; exit 1; }
done
