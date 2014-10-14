function wt = triang_wt(x)

  mu = 2;
  sigma = 1;

  mu_eff = mu - x;
 
  wt = mu_eff + sigma*normpdf(-mu_eff/sigma)/(1 - normcdf(-mu_eff/sigma));
