function [pll pl pl2 YY XX] = pred_llik_sngp(Ytest, rslt, var_x)
% PRED_LLIK_SNGP Compute predictive log-likelihood and likelihood for 
% held-out points for sngp.
%
% [pll pl pl2] = pred_llik_sngp(Ytest, rslt, var_x)

YY = cat(2, Ytest{:});
XX = [];
for i = 1:numel(Ytest)
  XX = [XX i*ones(1,numel(Ytest{i}))];
end
N = numel(YY);

pll = zeros(1,N);
pl = zeros(1,N);
pl2 = zeros(1,N);

nzinds = find(~cellfun(@isempty, rslt(:,1)'));

for ii = nzinds % unique numbers of clusters exhibited
  for jj = 1:size(rslt{ii,1},1) % samples with ii clusters
    means = rslt{ii,1}(jj,:);
    %m_c = rslt{ii,2}(jj,:);
    n_c = squeeze(rslt{ii,3}(jj,:,:));
    %clsubreg = rslt{ii,4}(jj,:);

    nnani = ~isnan(means);
    means = means(nnani);
    n_c = n_c(:,nnani);

    for n = 1:N
      xi = XX(n);
      pp = n_c(xi,:) ./ sum(n_c(xi,:));
      if any(isinf(pp))
        fprintf('Not a measure!\n');
        keyboard;
      end
      
      tmp = sum( pp.*normpdf(YY(n), means, sqrt(var_x)) );

      pll(n) = pll(n) + log(tmp);
      pl(n) = pl(n) + tmp;
      pl2(n) = pl2(n) + tmp.^2;
    end
  end
end
