function p = get_join_prob(var_x, var_theta, measure, v_s, cl_sr, m_s, X, c,...
                           num_pts, sum_pts)

% We are given the set of data values and one of two cluster assignments. We 
% want the probability that the data points joined the restaurant the way they 
% did. We get this prob. upto a constant which depends on the other data points
% NOTE: all indexed points must belong to the same document !!!

  p = 0;
  
  for i = 1:2
    num_pts(i)    = num_pts(i) - sum(c == i);
    sum_pts(i)    = sum_pts(i) - sum(X(c == i));
    m_s(cl_sr(i)) = m_s(cl_sr(i)) - sum(c == i);
  end;

  for i=1:length(X)
    c_curr = c(i);
    s_curr = cl_sr(c_curr);

    mean_theta_post = sum_pts(c_curr) ./ (num_pts(c_curr) + var_x/var_theta);
    var_theta_post  = var_x*var_theta ./ (var_x + num_pts(c_curr).*var_theta);

    mean_post =(X(i).*var_theta_post + mean_theta_post.*var_x)./ ...
                (var_theta_post + var_x);
    var_post  = var_theta_post .* var_x ./ (var_theta_post + var_x);

    p = p + 0.5 * (log(var_post) -log(var_theta_post) + ...
                  mean_post^2/var_post - mean_theta_post^2/var_theta_post );

    if(num_pts(c_curr) == 0)
      num = measure(s_curr);
    else
      num = num_pts(c_curr);
    end;

    num = log(num) - log((measure(s_curr) + m_s(s_curr)).*(m_s(s_curr) > 0)...
          + (m_s(s_curr) == 0)) ;
    p = p + num + log(v_s(s_curr));
    num_pts(c_curr) = num_pts(c_curr) + 1;
    sum_pts(c_curr) = sum_pts(c_curr) + X(i);
    m_s(s_curr) = m_s(s_curr) + 1;

    if(isinf(p))
      error "Huh!!"
    end;

    if(length(p)> 1)
      error "Huh!!"
    end;
  end;
