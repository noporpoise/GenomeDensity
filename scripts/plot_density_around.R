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
  # Check path to out_file exists
  out_dir=dirname(out_file);
  if(!file.exists(out_dir)) {
    dir.create(out_dir, showWarnings=FALSE, recursive=TRUE);
  }

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
