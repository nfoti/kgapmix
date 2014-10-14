function [cl_updt, m_s_updt, num_pts_updt, sum_pts_updt] = ...
     restr_gibbs_samp(var_x, var_theta, msr, v_s_c, x_lst, num_pts, cls_sr,...
                      sum_pts, i, cl, m_s)

     num_pts(cl) = num_pts(cl) - 1;
     sum_pts(cl) = sum_pts(cl) - x_lst(i);
     m_s(cls_sr(cl)) = m_s(cls_sr(cl)) - 1;
     if(sum(m_s < 0) ~= 0)
       error "Dude??!!";
     end;

     mean_theta_post = sum_pts ./ (num_pts + var_x/var_theta);
     var_theta_post  = var_x*var_theta ./ (var_x + num_pts.*var_theta);

     mean_post =(x_lst(i).*var_theta_post + mean_theta_post.*var_x)./ ...
                 (var_theta_post + var_x);
     var_post  = var_theta_post .* var_x ./ (var_theta_post + var_x);

     wt = ( (mean_post.^2)./var_post - ...
            (mean_theta_post.^2)./var_theta_post );

     wt = 0.5 .* (log(var_post) - log(var_theta) + wt);

     num_pts_updt = num_pts;
     num_pts  = num_pts + msr.*(num_pts == 0);

     num_pts = num_pts ./( (msr + m_s).*(m_s > 0) + (m_s == 0)) ;

     p_vec    = log(num_pts) + log(v_s_c) + wt;  
     p_vec    = p_vec - logsumexp(p_vec);
     sel_clst = sampleDiscrete(exp(p_vec),1);

     sum_pts(sel_clst)      = sum_pts(sel_clst) + x_lst(i);
     num_pts_updt(sel_clst) = num_pts_updt(sel_clst) + 1;
     m_s(cls_sr(sel_clst)) = m_s(cls_sr(sel_clst)) + 1;

     if(sum(m_s < 0) ~= 0)
       error "Dude????";
     end;
     cl_updt = sel_clst;
     sum_pts_updt = sum_pts;
     m_s_updt = m_s;

