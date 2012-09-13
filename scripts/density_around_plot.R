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

plot_counts_around = function(file,title,xaxis,yaxis,out_file)
{
  data = read.csv(file=file, sep=',', as.is=T, header=T);
  #bin_width = data[2,'bin'] - data[1,'bin'];

  pdf(file=out_file, width=12, height=6);
  set_pars();

  bins_before_zero=data[,'bin']<0
  bins_after_zero=data[,'bin']>0

  plot(0, 0, type='n', main=title, xlab=xaxis, ylab=yaxis,
       xlim=c(data[1,'bin'],data[nrow(data),'bin']),
       ylim=c(0, max(data[bins_before_zero | bins_after_zero,'count'])));

  points(data[bins_before_zero,'bin'], data[bins_before_zero,'count'], type='l');
  points(data[bins_after_zero,'bin'], data[bins_after_zero,'count'], type='l');

  # Draw dot on zero
  zero_index = which(data[,'bin']==0);
  points(0, data[zero_index,'count'], type='b');
  # and line
  abline(v=0, col="lightgray", lty=3);

  dev.off();
}


plot_rates_around = function(file,rate_per,title,xaxis,yaxis,out_file)
{
  data = read.csv(file=file, sep=',', as.is=T, header=T);
  bin_width = data[2,'bin'] - data[1,'bin'];

  pdf(file=out_file, width=12, height=6);
  set_pars();

  data[,'rate'] = (data[,'count'] * rate_per) / (bin_width * data[,'density']);

  plot(data[,'bin'], data[,'rate'], type='l',
       main=title, xlab=xaxis, ylab=yaxis,
       ylim=c(0,max(data[,'rate'])));

  # Draw dot on zero
  zero_index = which(data[,'bin']==0);
  points(0, data[zero_index,'rate'], type='b');
  # and line
  abline(v=0, col="lightgray", lty=3);

  dev.off();
}

plot_indels = function(num_of_bins, bin_size, csv_base_name, title_name,
                       output_dir, output_prefix, max_indel_size)
{
  short_types=c('ins','del');
  long_types=c('insertions', 'deletions');

  for(indel_size in 1:10)
  {
    for(type in 1:2)
    {
      short_type=short_types[type];
      long_type=long_types[type];

      file=paste(csv_base_name,'.',num_of_bins,'bins.',bin_size,'width.',indel_size,short_type,'.csv',sep='');
      title=paste(title_name,' ',indel_size,'bp ',long_type,', ',num_of_bins,' x ',bin_size,'bp bins',sep='');
      xaxis=paste(indel_size,'bp ',long_type,' per ',bin_size,'bp',sep='');
      out_file_counts=paste(output_dir,'/counts/',output_prefix,'.',indel_size,short_type,'.counts.pdf',sep='');
      out_file_rates=paste(output_dir,'/rates/',output_prefix,'.',indel_size,short_type,'.rates.pdf',sep='');

      plot_counts_around(file, title, 'Dist from gene TSS/TSE', xaxis, out_file_counts);
      plot_rates_around(file, bin_size, title, 'Dist from gene TSS/TSE', xaxis, out_file_rates);
    }
  }
}
