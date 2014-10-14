function [rslt,diagn,KK,ittime] = mcmc_silent(X,tt,L,burnin,num_samp,var_x,var_theta,measure)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  doc_subreg{doc} Cell array indicating which subregions belong to doc. 
%  subreg_doc{sub} Cell array indicating which docs subregion s belongs to. 
%  c{doc}(word) contains the cluster assignment of word
%
% State of the Markov Chain

%  n_s(doc,sub) is the number of words in subregion sub of document doc
%  m_s(sub) is the number of words in subregion sub 
%     (equals sum(n_s(doc,sub),1))
%  n_c(doc,clstr) is the number of words in cluster clstr of document doc
%  m_c{clstr} is the number of words in cluster clstr 
%     (equals sum(n_c(doc,sub),1))
%  clstr_subrg indicates to which subgroup each cluster belongs
%  sum_pts_cls sum of all points assigned to that cluster
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The returned structure consists of 5 columns corresponding to
% [mean, m_c, n_c, clstr_subrg, v_s, iter_info]
%
% The i'th row corresponds to the case of i clusters
% The corresponding elements of rslt are arrays of the values 
% of theta/m_c/n_c in all Gibbs sampling iterations with i clusters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set variables/parameters that don't change

%{
% foti
% This is the original hard coded example.  More flexible code below.
  num_docs = 4;
% num_subreg = 6;
% doc_subreg = {
%                [1,2,4];
%                [4,5,2,3,1];
%                [5,6,2,3,1];
%                [6,3,1]
%               };

% measure(1:num_subreg) = [.9,.3,.3,.1,.1,.1];

   num_subreg = 5;
   doc_subreg = { 
                 [1,2];
                 [2,3];
                 [3,4];
                 [4,5]
                };
   measure(1:num_subreg) = .3;

  for i = 1:num_subreg
    tmp = [];
    for j = 1:num_docs
      if sum(doc_subreg{j} == i)
        tmp = [tmp,j];
      end;
    end;
    subreg_doc{i} = tmp;
  end;
%}
  
  num_iter = burnin + num_samp;

  max_num_clstrs = 1000;
  tsampleits = 100;

  num_docs = numel(X);
  [doc_subreg, subreg_doc] = build_regions_sngp(tt, L);
  num_subreg = max(cat(2,doc_subreg{:}));

  %measure(1:num_subreg) = .3; % shape parameter of GP on region R (loose
                              % correspondance with \alpha in kgap I guess)
  measure = measure.*ones(1,num_subreg); % only allow scalar measure for now


  % passed in now
  %var_x      = 1;   % Variance of each mixture component (1/phi in kgap)
  %var_theta  = 10;    % Variance of parameter. This defines the measure on the 
                     % space (1/thtau1 in kgap)


  % Diagnostics
  lik = [0];
  m_s_hist = [];
  v_s_hist = [];

  smp_pts = 0;
  %{
  ct = 0;
  while(1)
    smp_pts  = randperm(400);
    smp_pts  = smp_pts(1:20);
    if(var(X{2}(smp_pts)) > 10) break; end;
    ct = ct + 1;
  end;
  %}

  % passed in now
  %burnin   = 0;%2000;
  %num_iter = 120;%6000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize empty/0 cells/arrays

  num_samples = ones(1,max_num_clstrs);

  % Each of the following variables are cell arrays, with element i
  % corresponding to the history of configurations with i clusters

  ret_mean  = cell(max_num_clstrs,1);
  ret_mc    = cell(max_num_clstrs,1);
  ret_nc    = cell(max_num_clstrs,1);
  ret_cl_sb = cell(max_num_clstrs,1);
  ret_v     = cell(max_num_clstrs,1);

  KK = zeros(1,num_iter-burnin);
  ittime = [];
  num_clstrs = 0;         % Number of clusters OVERALL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize remaining stuff randomly

  v_s = gamrnd(measure,1);

  %%%%%%%%%%%%%%%%%%%%%
  % Clusters in each subregion

  clstr_subrg = [];
  for i=1:num_subreg
    cnt =  unidrnd(20,1,1)+3;          % Number of clusters in each subgroup 
                                       % (is there a better init??)

    clstr_subrg = [clstr_subrg,i(ones(1,cnt))]; % Assign clusters to subgrp

    num_clstrs = num_clstrs + cnt;
  end;

  for doc = 1:num_docs
   % Get clstr IDs of all clstrs in subregions encompassed by the current doc 
    vld_clstrs = find(ismember(clstr_subrg(1:num_clstrs),doc_subreg{doc}) > 0);

    % Randomly assign words in doc from among these clusters
    c(doc) = {vld_clstrs(unidrnd(length(vld_clstrs),1,length(X{doc})))};
  end;

  % Given the cluster assignments, now calculate counts
  for doc = 1:num_docs
    % n_c(doc,i) gives the number of elements in the ith overall cluster in doc,
    %            if cluster i isn't present in doc, n_c(doc,i) = 0

    n_c(doc,1:num_clstrs)    = histc(c{doc},1:num_clstrs);

    subreg = clstr_subrg(c{doc});
    n_s(doc,1:num_subreg) = histc(subreg,1:num_subreg);
  end;

  for cl = 1:num_clstrs

    % Find all data points that belong to this cluster
    pts = [];
    for doc = 1:num_docs
      pts = [pts,X{doc}(find(c{doc} == cl))];
    end;

    sum_pts_cls(cl)    = sum(pts);

  end;
 
  % m_c/m_s give cluster/subregion counts across documents
  m_c = sum(n_c,1);
  m_s = sum(n_s,1);

  % Posterior mean of cluster parameter given data points is the 
  % empirical mean with some shrinkage
  %
  mean_pts_cls = sum_pts_cls ./ (m_c + var_x/var_theta);

  % Posterior variance of cluster parameter given data points
  var_pts_cls = var_x*var_theta ./ (var_x + m_c.*var_theta);

  sm_stat = zeros(2,3);
  mv_stat = zeros(1,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  fprintf('Sampling...\n');

  for iter = 1:num_iter
  
    %foti - restart timer if just collected a time (every tsampleits iterations)
    if mod(iter,tsampleits) == 1
      tstart = tic;
    end

    lik = [lik,get_data_prob(num_docs, doc_subreg, measure, var_x, var_theta,...
                             clstr_subrg, v_s, m_s, X, c, m_c, sum_pts_cls)];
    if lik(end-1)-lik(end) > 5000 && iter > 2
      error "Something weird happened!"
    end;
    %if(rem(iter,100) == 0)
    %  iter
    %  num_clstrs
    %  clstr_subrg
    %  v_s
    %  m_s + measure
    %  1./shp
    %end;
    %iter
    
    if mod(iter,100) == 1
      fprintf('iter: %d / %d\n', iter, num_iter);
    end

    if(rem(iter,5) ~= 0)
      %fprintf('Metropolis-Hastings split ahoy!\n');

      [num_clstrs,clstr_subrg,c,m_c,n_c,m_s,n_s,sum_pts_cls,mean_pts_cls,...
        var_pts_cls,sm_stat] = ...
         mh_spl_mrg(num_docs,num_clstrs,doc_subreg,measure,var_theta,var_x,...
                  clstr_subrg,v_s,c,X,m_c,n_c,m_s,n_s,sum_pts_cls,...
                  mean_pts_cls,var_pts_cls,sm_stat);
      
      %sm_stat

%    elseif(rem(iter,2) == 1)
    else

      m_s_hist = [m_s_hist,m_s'];
      v_s_hist = [v_s_hist,v_s'];

      %fprintf('Metropolis-Hastings move ahoy!\n');
      
      [clstr_subrg, m_s, n_s, mv_stat] = mh_move(num_docs,num_clstrs, ...
                doc_subreg, measure, clstr_subrg, v_s, c, m_c, m_s, n_s, ...
                mv_stat);
      %mv_stat
    end;

    if(rem(iter,25) == 0)
    %if true
      for g_iter = 1:20
        for doc = 1:num_docs  %VINAYAK
          for i = 1:length(X{doc})

          %%%%%%%%%%%%%%%%%%%% 
          % Get current point, its cluster assignment and the subregion this
          % cluster lies in
          %%%%%%%%%%%%%%%%%%%% 
            x_curr = X{doc}(i);
            c_curr = c{doc}(i);
            s_curr = clstr_subrg(c_curr);

          %%%%%%%%%%%%%%%%%%%% 
          % Counts excluding the current point
          %%%%%%%%%%%%%%%%%%%% 
            m_c(c_curr)     = m_c(c_curr) - 1;
            n_c(doc,c_curr) = n_c(doc,c_curr) - 1;
            m_s(s_curr)     = m_s(s_curr) - 1;
            n_s(doc,s_curr) = n_s(doc,s_curr) - 1;
        
            sum_pts_cls(c_curr) = sum_pts_cls(c_curr) - x_curr;

            % Posterior mean of cluster parameter given data points is the 
            % empirical mean with some shrinkage
            %
            mean_pts_cls(c_curr) = sum_pts_cls(c_curr) ./ ...
                                   (m_c(c_curr) + var_x/var_theta);

            % Posterior variance of cluster parameter given data points
            var_pts_cls(c_curr) = var_x*var_theta ./ ...
                                  (var_x + m_c(c_curr).*var_theta);

            mean_post = (x_curr.*var_pts_cls + mean_pts_cls.*var_x)...
                            ./(var_pts_cls + var_x);
        
            var_post = var_pts_cls .* var_x ./ (var_pts_cls + var_x);

            if(n_s(doc,s_curr) < 0)
              error "No!";
            end;
          %%%%%%%%%%%%%%%%%%%% 
          % Cluster assignment probabilities for the data point
          %%%%%%%%%%%%%%%%%%%% 
          
          % Only components in subregions associated with current doc

            s_r    = clstr_subrg(1:num_clstrs);
            p_indx = find(                                                 ...
                           ((ismember(s_r,doc_subreg{doc})) .* (m_c > 0))  ...
                          > 0);
            ln_p_vec  = log(m_c) + log(v_s(s_r)) - log(measure(s_r) + m_s(s_r));

            t_indx = find(ismember(1:num_subreg,doc_subreg{doc}));

            tail   = log(v_s(1:num_subreg)) + log(measure(1:num_subreg)) - ...
                      log(measure(1:num_subreg) + m_s(1:num_subreg)); 

            theta = ( (mean_post.^2)./var_post - ...
                      (mean_pts_cls.^2)./var_pts_cls );

            theta = 0.5 .* (log(var_post) - log(var_pts_cls) + theta);

            var_post_tail  = var_x*var_theta/(var_x + var_theta);  % and 
            tail_theta = x_curr*x_curr.*var_theta./(var_x*(var_x+var_theta))...
                          .* ones(1,num_subreg);

            tail_theta = 0.5 .* (log(var_post_tail) - log(var_theta) + ...
                                 tail_theta);


            % Select a cluster to assign x_curr to
            p_vec = [ln_p_vec(p_indx), tail(t_indx)] + ...
                    [theta(p_indx), tail_theta(t_indx)];
            p_vec = p_vec - logsumexp(p_vec);
            sel_clst = sampleDiscrete(exp(p_vec),1);

            % If the selected cluster is a new cluster

            if sel_clst > length(p_indx)
              sel_clst = t_indx(sel_clst - length(p_indx));
              num_clstrs = num_clstrs + 1;
              sel_clstr_corr = num_clstrs;
              sel_subreg     = sel_clst;

              n_c(doc,sel_clstr_corr) = 1;
              m_c(sel_clstr_corr)     = 1;
              sum_pts_cls(sel_clstr_corr) = x_curr;
            else
              sel_clst = p_indx(sel_clst);
              sel_clstr_corr = sel_clst;
              sel_subreg     = clstr_subrg(sel_clst);

              sum_pts_cls(sel_clstr_corr) = sum_pts_cls(sel_clstr_corr)+x_curr;
              n_c(doc,sel_clstr_corr) = n_c(doc,sel_clstr_corr) + 1;
              m_c(sel_clstr_corr)     = m_c(sel_clstr_corr) + 1;
            end;

            mean_pts_cls(sel_clstr_corr) = sum_pts_cls(sel_clstr_corr) ./ ...
                                   (m_c(sel_clstr_corr) + var_x/var_theta);

            var_pts_cls(sel_clstr_corr) = var_x*var_theta ./ ...
                                  (var_x + m_c(sel_clstr_corr).*var_theta);


            c{doc}(i) = sel_clstr_corr;
            clstr_subrg(sel_clstr_corr) = sel_subreg;

            n_s(doc,sel_subreg) = n_s(doc,sel_subreg) + 1;
            m_s(sel_subreg)     = m_s(sel_subreg) + 1;

            % Clean up old cluster if it has 0 members by exchanging it with the
            % last cluster
            if(m_c(c_curr) == 0)
              m_c(:,c_curr) = m_c(:,num_clstrs);
              n_c(:,c_curr) = n_c(:,num_clstrs);
              clstr_subrg(c_curr)  = clstr_subrg(num_clstrs);
              sum_pts_cls(c_curr)  = sum_pts_cls(num_clstrs);
              mean_pts_cls(c_curr) = mean_pts_cls(num_clstrs);
              var_pts_cls(c_curr)  = var_pts_cls(num_clstrs);

            % Repoint all words that belonged to the last cluster to its new
            % cluster value  
              for d = 1:num_docs
                c{d}(c{d} ==  num_clstrs) = c_curr;
              end;

              num_clstrs = num_clstrs - 1;
              m_c          = m_c(:,1:num_clstrs);
              n_c          = n_c(:,1:num_clstrs);
              clstr_subrg  = clstr_subrg(1:num_clstrs);
              sum_pts_cls  = sum_pts_cls(1:num_clstrs);
              mean_pts_cls = mean_pts_cls(1:num_clstrs);
              var_pts_cls  = var_pts_cls(1:num_clstrs);
            end;

            %if(sum((m_c == 0)))
            %  error "no!";
            %end;

          end;

        end;  %for doc
      end;
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gibbs sampling for V_s
    
    for doc = 1:num_docs
      n_l      = sum(n_s(doc,:));
      z_l(doc) = gamrnd(n_l,1./sum(v_s(doc_subreg{doc})));
    end;

    for s = 1:num_subreg
      l_s  = subreg_doc{s}; 
      shp      = 1 + sum(z_l(l_s));
      v_s(s) = gamrnd(m_s(s) + measure(s),1./shp);
    end;

    %foti
    if mod(iter,tsampleits) == 0
      telapsed = toc(tstart);
      ittime = [ittime telapsed];
    end

    if(iter == burnin) 
      sm_stat = sm_stat .* 0;
      mv_stat = mv_stat .* 0;
    end;

    if(iter > burnin) 

      %foti
      KK(iter-burnin) = num_clstrs;
      ret_mean {num_clstrs}(num_samples(num_clstrs),:) = (sum_pts_cls./m_c); 
      ret_mc   {num_clstrs}(num_samples(num_clstrs),:) = (m_c); 
      ret_nc   {num_clstrs}(num_samples(num_clstrs),:,:) = (n_c); 
      ret_cl_sb{num_clstrs}(num_samples(num_clstrs),:) = (clstr_subrg); 
      ret_v    {num_clstrs}(num_samples(num_clstrs),:) = (v_s); 
    
      % [The number clusters currently, and the count of how many times this
      %  has occurred so far] (the latter is redundant)
      ret_times(iter-burnin,:) = [num_clstrs,num_samples(num_clstrs)];

      % The cluster assignment of smp_pts
      %foti - had to set this to 0
      ret_clstr_hist(iter-burnin,:) = 0;%c{2}(smp_pts);
      num_samples(num_clstrs) = num_samples(num_clstrs) + 1;

%      rslt = [ret_mean, ret_mc, ret_nc, ret_cl_sb, ret_v];
    end;


  end;
  rslt  = [ret_mean, ret_mc, ret_nc, ret_cl_sb, ret_v];
  diagn = {sm_stat, mv_stat, lik, m_s_hist, v_s_hist, ret_times, smp_pts, ...
           ret_clstr_hist};
