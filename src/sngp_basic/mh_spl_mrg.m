% Metropolis-Hastings step to exchange points between two clusters in the same
% document. We first pick 2 points at random. If these belong to the same 
% cluster, then we create a new cluster with some points of the associated
% document (using restricted gibbs sampling to improve acceptance probability)
% If they belong to diffent clusters, one of them is a "pure" cluster, then
% we merge it into the other


function [num_clstrs, clstr_subrg, c,m_c,n_c,m_s,n_s,sum_pts_all, ...
            mean_pts_cls, var_pts_cls,stats] = ...
          mh_spl_mrg ...
            (num_docs,num_clstrs, doc_subreg,measure, var_theta, ...
               var_x, clstr_subrg,v_s,c,X,m_c,n_c,m_s,n_s,sum_pts_all,...
               mean_pts_cls,var_pts_cls,stats)

  if(length(mean_pts_cls) ~= num_clstrs)
    error "What??" ;
  end;
  num_iters = 10;
  %%%%%%%%%%%%%%%
  % Pick a document at random
  %
  doc_id = 1 + floor(num_docs*rand*(1-1e-15));
  %%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%
  % Pick a pair of words in this doc at random
  %
  i1 = (1+floor((length(X{doc_id})  )*rand*(1-1e-15)));
  i2 = (1+floor((length(X{doc_id})-1)*rand*(1-1e-15)));

  if(i2 >= i1) 
    i2 = i2 + 1;
  end;

  c1 = c{doc_id}(i1);
  c2 = c{doc_id}(i2);

%  disp('********************');
  if(c1 == c2)
%    disp('Split');
    sp_mrg = 1;      % split
    stats(2,2) = stats(2,2) + 1;
  else 
%    disp('Merge');
    sp_mrg = 0;
    stats(1,2) = stats(1,2) + 1;
  end;
%  disp('********************');
%
  % Store initial assignments and generate launch state

  if(sp_mrg == 1)

    % Pointers to the variables that are going to be flipped around
    indx1 = find(c{doc_id} == c1);
    indx2 = [];

    cl_orig = ones(1,length(indx1));
    num_pts_orig = [m_c(c1),0];
    sum_pts_orig = [sum_pts_all(c1),0];

%   % Is this picked cluster clean (ie no points from other documents?)
%   if(num_pts_orig(1) == length(indx1))
%     clean = 1;
%   else
%     clean = 0;
%   end;
    old_sub = 0;
    m_s_orig = m_s(clstr_subrg(c1));
    msr_orig = measure(clstr_subrg(c1));
    v_s_orig = v_s(clstr_subrg(c1));
    cl_sr_orig = [1,1];

  else
    % Pointers to the variables that are going to be flipped around
    indx1 = find(c{doc_id} == c1);
    indx2 = find(c{doc_id} == c2);

    if(length(indx1) ~= m_c(c1) && length(indx2) ~= m_c(c2))
%      disp('Reject: no pure clusters');
      % debug
%      length(indx1)
%      m_c(c1)
%      length(indx2)
%      m_c(c2)
      stats(1,3) = stats(1,3) + 1;
      % debug
 
      return;
    else
      % Atleast 1 cluster is clean

      if(length(indx1) == m_c(c1) && length(indx2) == m_c(c2))
        % Both are clean
%        clean = 1;
      else
%        clean = 0;
        % Make c2 the clean cluster
        if(length(indx1) == m_c(c1))
          tmp = indx1;
          indx1 = indx2;
          indx2 = tmp;
         
          tmp = c1;
          c1  = c2;
          c2  = tmp;
        end;
 
      end;
    end;

    old_sub = clstr_subrg(c2);
    cl_orig = [ones(1,length(indx1)), 2.*ones(1,length(indx2))];
    num_pts_orig = m_c([c1,c2]);
    sum_pts_orig = sum_pts_all([c1,c2]);

    m_s_orig = m_s(clstr_subrg([c1,c2]));
    msr_orig = measure(clstr_subrg([c1,c2]));
    v_s_orig = v_s(clstr_subrg([c1,c2]));
    cl_sr_orig = [1,2];

  end;

  % Local copies of the x-values and the cluster ids. These shall be written
  % back at the end
  indx   = [indx1,indx2];
  x_lst  = X{doc_id}(indx);

  % New subregion where c2 lies
  new_sub = doc_subreg{doc_id}(1 + floor(length(doc_subreg{doc_id})* ...
                                        rand*(1-1e-15)));

  if(clstr_subrg(c1) == new_sub)
    sub_reg = new_sub;
    cl_sr = [1,1];
  else
    sub_reg = [clstr_subrg(c1), new_sub];
    cl_sr = [1,2];
  end;

  msr     = measure(sub_reg);
  v_s_c   = v_s(sub_reg);
  m_s_loc = m_s(sub_reg);

  if new_sub ~= old_sub
    m_s_loc(cl_sr(2)) = m_s_loc(cl_sr(2)) + num_pts_orig(2);
  end;

  cl_updt = sampleDiscrete([0.5,0.5],length(indx));

  m_s_updt = m_s_loc;
  for i = 1:2
    num_pts(i) = num_pts_orig(i) + (sum(cl_updt == i) - sum(cl_orig == i));

    sum_pts(i) = sum_pts_orig(i) + (sum(x_lst(cl_updt == i)) - ...
                                    sum(x_lst(cl_orig == i)));
    
    m_s_updt(cl_sr(i)) = m_s_updt(cl_sr(i)) + ...
                         (sum(cl_updt == i) - sum(cl_orig == i));
    %(sum(cl_updt == i) - sum(cl_orig == i));
  end;


  for iters = 1:num_iters-1

    % One round of Gibbs updating
    for i = 1:length(x_lst)
      [cl_updt(i), m_s_updt, num_pts, sum_pts] = ...
        restr_gibbs_samp(var_x, var_theta, msr, v_s_c, x_lst, num_pts, ...
                         cl_sr, sum_pts, i, cl_updt(i), m_s_updt);

       if(sum(num_pts) ~= sum(num_pts_orig));
         error "Feck!";
       end;

       if(sum(m_s_updt < 0) ~= 0);
         error "Fork!";
       end;
    end;
  end;

  cl_launch      = cl_updt; 
  num_pts_launch = num_pts;
  sum_pts_launch = sum_pts;
  m_s_launch     = m_s_updt;

  % Transition from launch state to final state
  if(sp_mrg == 1)
    % One round of Gibbs updating
    for i = 1:length(x_lst)
      [cl_updt(i), m_s_updt, num_pts, sum_pts] = ...
        restr_gibbs_samp(var_x, var_theta, msr, v_s_c, x_lst, num_pts, ...
                         cl_sr, sum_pts, i, cl_updt(i), m_s_updt);

    end;
  else 
    % everything goes to the dirty state
    num_pts = [sum(num_pts),0];
    sum_pts = [sum(sum_pts),0];

    if(cl_sr(2) == 2)
      m_s_updt = m_s_updt + [(num_pts(1) - sum(cl_updt == 1)), ...
                              -sum(cl_updt == 2)];
    end;

    cl_updt(:) = 1;
  end;

%  [sum(num_pts),sum(num_pts_orig)]
  if(sum(num_pts) ~= sum(num_pts_orig));
    error "Heck!";
  end;

  % Calculate probabilities
  if(sp_mrg == 1)
    % reject splits which just move a cluster
    if(sum(cl_updt == 1) == 0 || sum(cl_updt == 2) == 0)
%      disp('Reject: Just moving a cluster');
      return;
    else
      % from launch to final state
      p_fwd =  get_trans_prob(var_x, var_theta, msr, cl_sr, v_s_c, ...
                         m_s_launch, x_lst, cl_launch, num_pts_launch, ...
                         sum_pts_launch, cl_updt);

      % from launch to original state
      p_bck = 0;

%   1
%   if(clean == 0)
%     % If either of the clusters were 'dirty' (ie has words from other 
%     % documents), then you pick the clean one with probability 1
%     p_bck = 0;
%   else
%     % Else you place the new cluster in a random subset
%     p_bck = - log(length(doc_subreg(doc_id,:)));
%     %  p_bck = log(v_s_c(1)) - log(sum(v_s(doc_subreg(doc_id,:))));
%   end;
    end;

  else
    % from launch to final state
    p_fwd = 0;
%   1
%   if(clean == 0)
%     p_fwd = 0;
%   else
%     p_fwd = - log(length(doc_subreg(doc_id,:)));
%     %  p_fwd = log(v_s_c(2)) - log(sum(v_s(doc_subreg(doc_id,:))));
%   end;

    % from launch to original state
    p_bck =  get_trans_prob(var_x, var_theta, msr, cl_sr, v_s_c, ...
                       m_s_launch, x_lst, cl_launch, num_pts_launch, ...
                       sum_pts_launch, cl_orig);
  end;

  p_orig =  get_join_prob(var_x, var_theta, msr_orig, v_s_orig, cl_sr_orig, ...
                       m_s_orig, x_lst, cl_orig, num_pts_orig, sum_pts_orig);

  p_fin  =  get_join_prob(var_x, var_theta, msr, v_s_c, cl_sr, m_s_updt, ...
                       x_lst, cl_updt, num_pts, sum_pts);

%  [p_bck , p_fin , p_fwd , p_orig]'
%  [p_bck - p_fwd , p_fin - p_orig]'
  if isinf([p_bck , p_fin , p_fwd , p_orig])
    error "!!!"
  end;
  acc = exp(min(0,p_bck + p_fin - p_fwd - p_orig));

  if (rand > acc)
%    disp('Reject!');
    cl_updt = cl_orig;       %#ok
    num_pts = num_pts_orig;  %#ok
    sum_pts = sum_pts_orig;  %#ok
  else
%    disp('Accept!');
    stats(sp_mrg+1,1) = stats(sp_mrg+1,1) + 1;

%   n_c 
%   n_s

    chng = num_pts(1) - num_pts_orig(1);
    m_c(c1) = m_c(c1) + chng;
    n_c(doc_id,c1) = n_c(doc_id,c1) + chng;
    sum_pts_all(c1)  = sum_pts(1);
    m_s(clstr_subrg(c1)) = m_s(clstr_subrg(c1)) + chng;
    n_s(doc_id,clstr_subrg(c1)) = n_s(doc_id,clstr_subrg(c1)) + chng;

    if(n_c(doc_id,c1) < 0 || n_s(doc_id,clstr_subrg(c1)) < 0)
      error "No!";
    end;


    mean_pts_cls(c1) = sum_pts_all(c1) ./ (m_c(c1) + var_x/var_theta);
    var_pts_cls(c1) = var_x*var_theta ./ (var_x + m_c(c1).*var_theta);

    % DOESNT CHANGE IF 1 not implemented
    % clstr_subrg(c1) 

    if(sp_mrg == 1)
      num_clstrs = num_clstrs + 1;
      c2 = num_clstrs;

      % Now use cl_orig for another purpose
      cl_orig(cl_updt == 1) = c1;
      cl_orig(cl_updt == 2) = c2;

      c{doc_id}(indx)       = cl_orig;

%      % Repoint all words that belonged to the last cluster to its new
%      % cluster value

%      for d = 1:num_docs
%        c{d}(c{d} ==  c2) = num_clstrs;
%      end;

      clstr_subrg(c2) = new_sub;
      m_c(c2) = num_pts(2);
      % VINAYAK: make sure I don't have to init n_c for other docs
      n_c(doc_id,c2) = num_pts(2);
      sum_pts_all(c2)  = sum_pts(2);
      m_s(new_sub) = m_s(new_sub) + num_pts(2);
      n_s(doc_id,new_sub) = n_s(doc_id,new_sub) + num_pts(2);

      if(num_pts(2) < 0)
        error('What??');
      end;

      mean_pts_cls(c2) = sum_pts_all(c2) ./ (m_c(c2) + var_x/var_theta);
      var_pts_cls(c2) = var_x*var_theta ./ (var_x + m_c(c2).*var_theta);
    else

      m_s(clstr_subrg(c2))        = m_s(clstr_subrg(c2)) - num_pts_orig(2);
      n_s(doc_id,clstr_subrg(c2)) = n_s(doc_id,clstr_subrg(c2)) - ...
                                    num_pts_orig(2);
      m_c(c2)             = m_c(num_clstrs);
      n_c(:,c2)           = n_c(:,num_clstrs);
      clstr_subrg (c2)    = clstr_subrg(num_clstrs);
      sum_pts_all (c2)    = sum_pts_all(num_clstrs);
      mean_pts_cls(c2)    = mean_pts_cls(num_clstrs);
      var_pts_cls(c2)     = var_pts_cls(num_clstrs);

      c{doc_id}(c{doc_id} ==  c2) = c1;
      % Repoint all words that belonged to the last cluster to its new
      % cluster value
      for d = 1:num_docs
        c{d}(c{d} ==  num_clstrs) = c2;
      end;

      num_clstrs   = num_clstrs - 1;
      m_c          = m_c(1:num_clstrs);
      n_c          = n_c(:,1:num_clstrs);
      clstr_subrg  = clstr_subrg(1:num_clstrs);
      sum_pts_all  = sum_pts_all(1:num_clstrs);
      mean_pts_cls = mean_pts_cls(1:num_clstrs);
      var_pts_cls  = var_pts_cls(1:num_clstrs);

    %sum(m_c)
    end;
  end;
