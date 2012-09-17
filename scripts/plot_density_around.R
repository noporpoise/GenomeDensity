#
# Plot distribution of variants around genes
#
set_pars = function()
{
  par("cex.main"=1);
  par("cex.axis"=1);
  par("cex.lab"=1);
  par("font.main"=1);
  par("font.axis"=1);
  par("font.lab"=1);
  par("font.sub"=1);
  par("cex"=1.1);
}


plot_around = function(file,title,xaxis,yaxis,out_file,do_rates=F)
{
  data = read.csv(file=file, sep=',', as.is=T, header=T);
  bin_width = data[2,'bin'] - data[1,'bin'];

  pdf(file=out_file, width=12, height=6);
  set_pars();

  if(do_rates)
  {
    field='rate'
    ybounds=c(0,max(data[,field]));

    data[,'rate'] = data[,'count'] / data[,'density'];
    #data[,'rate'] = (data[,'count'] * rate_per) / (bin_width * data[,'density']);
  }
  else
  {
    field='count'
    ybounds=c(0, max(data[bins_before_zero | bins_after_zero,'count']))
  }

  xbounds=c(data[1,'bin'],data[nrow(data),'bin'])

  plot(0, 0, type='n', main=title, xlab=xaxis, ylab=yaxis,
       xlim=xbounds, ylim=ybounds);

  points(data[bins_before_zero,'bin'], data[bins_before_zero,field], type='l');
  points(data[bins_after_zero,'bin'], data[bins_after_zero,field], type='l');

  # Draw dot on zero
  zero_index = which(data[,'bin']==0);
  points(0, data[zero_index,field], type='b');
  # and line
  abline(v=0, col="gray", lty=3);

  dev.off();
}

# Internal function
do_plot_variants = function(num_of_bins, bin_size, csv_base_name, title_name,
                            output_dir, output_prefix, short_type, long_type)
{
  file=paste(csv_base_name,'.',num_of_bins,'bins.',bin_size,'width.',short_type,'.csv',sep='');
  title=paste(title_name,' ',long_type,', ',num_of_bins,' x ',bin_size,'bp bins',sep='');
  xaxis=paste(long_type,' per ',bin_size,'bp',sep='');
  out_file_counts=paste(output_dir,'/counts/',output_prefix,'.',num_of_bins,'bins.',bin_size,'width.',short_type,'.counts.pdf',sep='');
  out_file_rates=paste(output_dir,'/rates/',output_prefix,'.',num_of_bins,'bins.',bin_size,'width.',short_type,'.rates.pdf',sep=''); 

  plot_counts_around(file, title, 'Dist from gene TSS/TSE', xaxis, out_file_counts);
  plot_rates_around(file, bin_size, title, 'Dist from gene TSS/TSE', xaxis, out_file_rates);
}


plot_variants = function(num_of_bins, bin_size, csv_base_name, title_name,
                         output_dir, output_prefix, max_indel_size)
{
  short_types=c('ins','del');
  long_types=c('insertions', 'deletions');

  for(indel_size in 1:10)
  {
    for(type in 1:2)
    {
      short_type = paste(indel_size, short_types[type], sep='');
      long_type = paste(indel_size, 'bp ', long_types[type], sep='');

      do_plot_variants(num_of_bins, bin_size, csv_base_name, title_name,
                       output_dir, output_prefix, short_type, long_type);
    }
  }

  do_plot_variants(num_of_bins, bin_size, csv_base_name, title_name,
                   output_dir, output_prefix, 'snps', 'SNPs');
}
