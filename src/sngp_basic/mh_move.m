% Metropolis-Hastings step to move a "pure" cluster from 1 subreg to another


function [clstr_subrg,m_s,n_s,mv_stat] = mh_move(num_docs,num_clstrs, ...
               doc_subreg, measure, clstr_subrg, v_s, c, m_c, m_s, n_s, ...
               mv_stat)

  mv_stat(2) = mv_stat(2) + 1;
  % Pick a document at random
  doc_id = 1 + floor(num_docs*rand*(1-1e-15));

  % Pick a word in this doc at random
  i  = (1+floor((length(c{doc_id}))*rand*(1-1e-15)));
  cl = c{doc_id}(i);
  old_sub = clstr_subrg(cl);

  indx = find(c{doc_id} == cl);

  % Is this picked cluster clean (ie no points from other documents?)
  if(length(indx) ~= m_c(cl))
%    disp('Reject: not a pure');
    return;
  end;

  % New subregion where c1 is moved to
  old_sr_id = find(doc_subreg{doc_id} == old_sub);
  sr_id = floor(1 + (length(doc_subreg{doc_id}) - 1) * rand*(1-1e-15));
  if sr_id >= old_sr_id
    sr_id = sr_id + 1;
  end;
  new_sub = doc_subreg{doc_id}(sr_id);
  msr     = measure([old_sub, new_sub]);
  v_s_c   = v_s([old_sub, new_sub]);

  denom(1) = sum(log(msr(1) + (m_s(old_sub) - m_c(cl)):(m_s(old_sub)-1)));
  denom(2) = sum(log(msr(2) + (m_s(new_sub)):(m_s(new_sub) + m_c(cl) - 1)));

  p_orig =  log(msr) + m_c(cl) .* log(v_s_c) - denom;
  p_fin  = p_orig(2);
  p_orig = p_orig(1);

  acc = exp(min(0 ,p_fin - p_orig));

  if(new_sub == old_sub)
    error "wtf??";
  end;
  if (rand > acc)
%    disp('Reject!');
  else
%    disp('Accept!');
%    [old_sub, new_sub]
%    disp('*****Accept! Moving from: to*****');

    m_s(old_sub) = m_s(old_sub) - m_c(cl);
    n_s(doc_id,old_sub) = n_s(doc_id,old_sub) - m_c(cl);
    m_s(new_sub) = m_s(new_sub) + m_c(cl);
    n_s(doc_id,new_sub) = n_s(doc_id,new_sub) + m_c(cl);

    clstr_subrg(cl) = new_sub;
    mv_stat(1) = mv_stat(1) + 1;
  end;
