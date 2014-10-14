% Metropolis-Hastings step to exchange points between to clusters in the same
% document. Includes restricted gibbs sampling to improve acceptance 
% probability

function [c,m_c,n_c,m_s,n_s,sum_pts_all] = ...
         mh_exchg(num_docs,num_clstrs, doc_subreg, measure, var_theta,var_x,...
                  clstr_subrg,v_s,c,X,m_c,n_c,m_s,n_s,sum_pts_all)

 num_iters = 10;
 %%%%%%%%%%%%%%%
 % Pick a document at random
 %
 doc_id = 1 + floor(num_docs*rand*(1-1e-15));
 %%%%%%%%%%%%%%%

 %%%%%%%%%%%%%%%
 % Pick a pair of clstrs proportional to their combined sizes
 %
 %%%
 % Find which of all clusters lie in subregions encompassed by doc the current
 %
 vld_clstrs = ...
     find(ismember(clstr_subrg(1:num_clstrs),doc_subreg(doc_id,:)) > 0);

 clstr_perms = combnk(vld_clstrs,2);

 if(rand > 0.5)
%   disp('Weighted selection');
   clstr_wts   = sum(m_c(clstr_perms),2);
   clstr_wts   = clstr_wts ./ sum(clstr_wts);
   cl_id = clstr_perms(sampleDiscrete(clstr_wts,1),:);
 else
%   disp('Selecting uniformly');
   cl_id = clstr_perms((1+floor(size(clstr_perms,1)*rand*(1-1e-15))),:);
 end;

 msr     = measure(clstr_subrg(cl_id));
 v_s_c   = v_s(clstr_subrg(cl_id));
 num_pts = m_c(cl_id);
 sum_pts = sum_pts_all(cl_id);

 % Pointers to the variables that are going to be flipped around
 indx1 = find(c{doc_id} == cl_id(1));
 indx2 = find(c{doc_id} == cl_id(2));

 % Data points that will NOT be flipped
 anc1 = 1 + floor(length(indx1)*rand*(1-1e-15));
 anc2 = 1 + floor(length(indx2)*rand*(1-1e-15));

 indx = [indx1(1:(anc1-1)), indx1((anc1+1):length(indx1)), ...
         indx2(1:(anc2-1)), indx2((anc2+1):length(indx2))];

 % Local copies of the x-values and the cluster ids. These shall be written
 % back at the end

 x_lst        = X{doc_id}(indx);

 cl_orig = [ones(1,length(indx1)-1), 2.*ones(1,length(indx2)-1)]; 
 num_pts_orig = num_pts;
 sum_pts_orig = sum_pts;

 cl_updt = sampleDiscrete([0.5,0.5],length(indx)); 
 num_pts(1) = num_pts(1) + (sum(cl_updt == 1) - sum(cl_orig == 1));
 num_pts(2) = num_pts(2) + (sum(cl_updt == 2) - sum(cl_orig == 2));

 sum_pts(1) = sum_pts(1) + (sum(x_lst(cl_updt == 1)) - ...
                            sum(x_lst(cl_orig == 1)));
 sum_pts(2) = sum_pts(2) + (sum(x_lst(cl_updt == 2)) - ...
                            sum(x_lst(cl_orig == 2)));

% disp('New round');
% [num_pts_orig; sum_pts_orig]
% [num_pts; sum_pts]
 for iters = 1:num_iters

   if iters == num_iters
     cl_launch      = cl_updt; 
     num_pts_launch = num_pts;
     sum_pts_launch = sum_pts;
   end;

   % One round of Gibbs updating
   for i = 1:length(x_lst)
     [cl_updt(i), num_pts, sum_pts] = ...
       restr_gibbs_samp(var_x, var_theta, msr, v_s_c, x_lst, num_pts, ...
                        sum_pts, i ,cl_updt(i));

   end;
   [num_pts; sum_pts]
 end;
 [num_pts; sum_pts]

 p_fwd =  get_trans_prob(var_x, var_theta, msr, v_s_c, x_lst, cl_launch, ...
                      num_pts_launch, sum_pts_launch, cl_updt);

 p_bck =  get_trans_prob(var_x, var_theta, msr, v_s_c, x_lst, cl_launch, ...
                      num_pts_launch, sum_pts_launch, cl_orig);
 
 p_orig =  get_join_prob(var_x, var_theta, msr, v_s_c, x_lst, cl_orig, ...
                      num_pts_orig, sum_pts_orig);
 
 p_fin =  get_join_prob(var_x, var_theta, msr, v_s_c, x_lst, cl_updt, ...
                      num_pts, sum_pts);


% [p_bck , p_fin , p_fwd , p_orig]
% [p_bck - p_fwd , p_fin - p_orig]
 acc = exp(min(0,p_bck + p_fin - p_fwd - p_orig));
 
 if (rand > acc)
%   disp('Reject!');
   cl_updt = cl_orig;
   num_pts = num_pts_orig;
   sum_pts = sum_pts_orig;
 else
%   disp('Accept!');
 end;

% Now use cl_orig for another purpose

 cl_orig(cl_updt == 1) = cl_id(1);
 cl_orig(cl_updt == 2) = cl_id(2);
 
 c{doc_id}(indx)       = cl_orig;
 m_c(cl_id(1))         = num_pts(1);
 m_c(cl_id(2))         = num_pts(2);
 sum_pts_all(cl_id(1)) = sum_pts(1);
 sum_pts_all(cl_id(2)) = sum_pts(2);

% sum(cnt_diff)

 cnt_diff = [num_pts(1) - num_pts_orig(1), num_pts(2) - num_pts_orig(2)];
 n_c(doc_id,cl_id) = n_c(doc_id,cl_id) + cnt_diff;

 sr = clstr_subrg(cl_id);

 for i=1:2
   m_s(sr(i))        = m_s(sr(i)) + cnt_diff(i);
   n_s(doc_id,sr(i)) = n_s(doc_id,sr(i)) + cnt_diff(i);
 end;
