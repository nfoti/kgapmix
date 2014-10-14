function plot_sim_mtrx(locns, cl_inf)

  clf
  cpr = (permute(cl_inf(:,:,ones(1,20)),[1,3,2]) == cl_inf(:,:,ones(1,20)));
  colormap(gray);
  v_val = subplot(3,3,1);
  barh(fliplr(locns));
  set(gca,'ylim',[0 length(locns)]+.5)
  im1 = subplot(3,3,2);
  image(256.*squeeze(mean(cpr,1)));
  h_val = subplot(3,3,5);
  bar(locns);
  set(gca,'xlim',[0 length(locns)]+.5)
  im2 = subplot(3,3,9);
  image(256.*squeeze(var(cpr,1)));

  % Make scatter plot bigger, histograms smaller
  set(v_val,'Position',[0.05 .6 .1 .35],'tag','yhist');
  set(h_val,'Position',[.2 .41 .35 .15],'tag','xhist');
  set(im1,'Position',[0.2 0.6 0.35 0.35],'tag','scatter');
  set(im2,'Position',[0.2 0.02 0.35 0.35],'tag','scatter2');

