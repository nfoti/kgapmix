function p = get_liklhd(num_docs, doc_subreg, measure, var_x, var_theta, ...
                        clstr_subrg, v_s, m_s, X, c, m_c, sum_pts)

  p = 0;

  for d=1:num_docs
    lg_v_sum = log(sum(v_s(doc_subreg{d})));

    for i=1:length(X{d})
    
      x_curr = X{d}(i);
      c_curr = c{d}(i);
      s_curr = clstr_subrg(c_curr);  

      p = p + log(v_s(s_curr)) - lg_v_sum;
      m_c(c_curr) = m_c(c_curr) - 1;
      if(m_c(c_curr) == 0)
        p = p + log(measure(s_curr));
      else
        p = p + log(m_c(c_curr));
      end;

      m_s(s_curr) = m_s(s_curr) - 1;
      if(m_s(s_curr) > 0)
        p = p - log(measure(s_curr) + m_s(s_curr));
      elseif(m_s(s_curr) < 0)
        error "what?!"
      end;
        
%     if(p>=0)
%       error "what?!"
%     end;

      sum_pts(c_curr) = sum_pts(c_curr) - x_curr;

      mean_theta_post = sum_pts(c_curr) ./ (m_c(c_curr) + var_x/var_theta);
      var_theta_post  = var_x*var_theta ./ (var_x + m_c(c_curr).*var_theta);

      mean_post =(x_curr.*var_theta_post + mean_theta_post.*var_x)./ ...
                  (var_theta_post + var_x);
      var_post  = var_theta_post .* var_x ./ (var_theta_post + var_x);

      p = p + 0.5 * (log(var_post) -log(var_theta_post) + ...
                    mean_post^2/var_post - mean_theta_post^2/var_theta_post );

    end;
  end;
