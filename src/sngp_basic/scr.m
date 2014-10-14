for olp = 1:50
  keep olp;
  gen_data;
  [rslt,it_times,p, sm_stat, mv_stat, lik,m_s_h,v_s_h] = gibbs_no_params(X);
  func_plot_noparams(m_c,n_c,clstr_subrg,v_s,clstr_mean,doc_subreg,rslt);
  str = sprintf('density%d.eps',olp)
  print('-depsc', str)
  clf;
  plot(lik)
  str = sprintf('lik%d.eps',olp)
  print('-depsc', str)
end;
