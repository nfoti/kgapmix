function p = get_trans_prob(var_x, var_theta, measure, cl_sr, v_s, ...
                            m_s_orig, x_list, c_orig, num_pts_orig, ...
                            sum_pts_orig, c_final)

  c_updt = c_orig;
  num_pts = num_pts_orig;
  sum_pts = sum_pts_orig;
  m_s     = m_s_orig;
  p = 0;

  for i = 1:length(c_orig)
    c_o = c_orig(i);
    c_f = c_final(i);

    s_o = cl_sr(c_o);
    s_f = cl_sr(c_f);

    num_pts(c_o) = num_pts(c_o) - 1;
    sum_pts(c_o) = sum_pts(c_o) - x_list(i);
    m_s(s_o)     = m_s(s_o) - 1;

    mean_theta_post = sum_pts ./ (num_pts + var_x/var_theta);
    var_theta_post  = var_x*var_theta ./ (var_x + num_pts(c_f).*var_theta);

    mean_post =(x_list(i).*var_theta_post + mean_theta_post.*var_x)./ ...
                (var_theta_post + var_x);
    var_post  = var_theta_post .* var_x ./ (var_theta_post + var_x);

    wt = 0.5 * (log(var_post) -log(var_theta_post) + ...
               mean_post.^2/var_post - mean_theta_post.^2/var_theta_post );

    if num_pts_orig(c_f) == 0
       num = measure(s_f);
    else
       num = num_pts_orig(c_f);
    end;

    if(m_s(s_f) == 0)
      cl_sel_log_prob = 0;
    else
      cl_sel_log_prob = log(num) - log(measure(s_f) + m_s(s_f));
    end;
    wt = wt + cl_sel_log_prob + log(v_s(s_f));

    wt = exp(wt);
    wt = (wt + 1e-30)./sum(wt);  % Avoid underflow errors if the reverse 
                                 % probability is really low

    p = p + log(wt(c_f));
    c_updt(i) = c_final(i);

    num_pts(c_f) = num_pts(c_f) + 1;
    sum_pts(c_f) = sum_pts(c_f) + x_list(i);
    m_s(s_f)     = m_s(s_f) + 1;
    if(sum(isinf(p)))
      error "Huh??"
    end;
    if(length(p)> 1)
      error "Huh??"
    end;
  end;
