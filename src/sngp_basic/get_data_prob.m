function p = get_data_prob(num_docs, doc_subreg, measure, var_x, var_theta, ...
                        clstr_subrg, v_s, m_s, X, c, m_c, sum_pts)

  p = 0;

  for d=1:num_docs
%    lg_v_sum = log(sum(v_s(doc_subreg{d})))  % Changed ;
    lg_v_sum = logsumexp(v_s(doc_subreg{d}));

    for i=1:length(X{d})
    
      x_curr = X{d}(i);
      c_curr = c{d}(i);
      s_curr = clstr_subrg(c_curr);  

%      p = p + log(v_s(s_curr)) - lg_v_sum;  % Changed 
      p = p + v_s(s_curr) - lg_v_sum;
      m_c(c_curr) = m_c(c_curr) - 1;
      if(m_c(c_curr) == 0)
        p = p + log(measure(s_curr));
      else
        p = p + log(m_c(c_curr));
      end;

      m_s(s_curr) = m_s(s_curr) - 1;
      if(m_s(s_curr) >= 0)
        p = p - log(measure(s_curr) + m_s(s_curr));
      elseif(m_s(s_curr) < 0)
        error "what?!"
      end;
        
      sum_pts(c_curr) = sum_pts(c_curr) - x_curr;

      p = p + get_lik(x_curr, m_c(c_curr), sum_pts(c_curr), var_theta, var_x);

    end;
  end;
