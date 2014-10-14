function func_plot_full_noparams(freq_tot,freq_doc,clstr_subrg,v_s,clstr_mean,doc_subrg,rslt)

% Usage:  func_plot_full(m_c,n_c,clstr_subrg,v_s,doc_subrg,rslt)
%
% The first plot has m_c (the overall)
% The rest have n_c for each document
%
% [mean, m_c, n_c, clstr_subrg, v_s]

  rg = [-25:.1:25];
  num_docs = size(freq_doc,1);
  var_x = 1;
  var_p = 10;
  clf
  for dc_lp = 1:num_docs+1
    dc_lp
    subplot(num_docs+1,1,dc_lp)
    hold on;

    num_clstr_range = size(rslt,1);
    for c = 1:num_clstr_range
      for g = 1:size(rslt{c,1},1)
        dw = rg .* 0;
        if dc_lp == 1
          efreq = rslt{c,2}(g,:)./sum(rslt{c,2}(g,:));
          efreq_dbg = rslt{c,2}(g,:)./sum(rslt{c,2}(g,:));
        else
          % Weight of each subregion in the current document
          wt   = rslt{c,5}(g,doc_subrg{dc_lp-1});
          wt   = wt ./ sum(wt);

          efreq = squeeze(rslt{c,2}(g,:) .* 0);

          efreq_dbg = squeeze(rslt{c,2}(g,:) .* 0);

          for s = 1:length(doc_subrg{dc_lp-1})
            
            % 2nd line:Find all clusters belonging to subregion s
            tmp = squeeze(rslt{c,2}(g,:).* ...
                    (rslt{c,4}(g,:) == doc_subrg{dc_lp-1}(s)));

            if length(tmp)> 0 & sum(tmp)>0
              efreq = efreq + wt(s) .* tmp ./ sum(tmp);
            end;

            if doc_subrg{dc_lp-1}(s) == 2
              efreq_dbg = efreq_dbg + wt(s) .* tmp ./ sum(tmp);
            end;
          end;
        end;
        for lp = 1:c 
          dw = dw + efreq(lp) * normpdf(rg , rslt{c,1}(g,lp), var_x+var_p/rslt{c,2}(g,lp)); 
        end;
        plot(rg,dw)
%        dc_sum = dc_sum + wt;
      end;
    end;
%    dc_sum

    if dc_lp == 1
      freq = freq_tot;
      freq = freq ./ sum(freq);
    else
      wt   = v_s(doc_subrg{dc_lp-1});
      wt   = wt ./ sum(wt);

      freq = freq_tot .* 0;
      for s = 1:length(doc_subrg{dc_lp-1})
        tmp = freq_tot .* (clstr_subrg == (doc_subrg{dc_lp-1}(s)));
        if (length(tmp)> 0 & sum(tmp)> 0) 
          freq = freq + wt(s) .* tmp ./ sum(tmp);
        end;
      end;
    end;

    dw = rg .* 0;
    for lp = 1:length(clstr_mean) 
      dw = dw + freq(lp) * normpdf(rg, clstr_mean(lp), var_x); 
      plot(rg,freq(lp) * normpdf(rg, clstr_mean(lp), var_x),'g')
    end;
    plot(rg,dw,'r')
  end;
