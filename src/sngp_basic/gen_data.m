num_docs = 4;

% num_subreg = 6;
% doc_subreg = {
%               [1,2,4];
%               [4,5,2,3,1];
%               [5,6,2,3,1];
%               [6,3,1];
%              };
%  measure(1:num_subreg) = [.4,.1,.1,.05,.05,.05];

 num_subreg = 5;
 doc_subreg = {
               [1,2];
               [2,3];
               [3,4];
               [4,5];
              };
 measure(1:num_subreg) = .3;

var_x = 1.0;
var_theta = 10;

size_l = 4.*[100, 100,100,100];

v_s = gamrnd(measure,1);

n_c = zeros(num_docs,1);
m_c = [0];
n_s = zeros(num_docs,1);
m_s = [0];
num_clstrs = 0;
clstr_subrg = [];

for doc = [1:num_docs]
  for i = [1:size_l(doc)]

    % Pick a subregion at random
    v = v_s(doc_subreg{doc});
    v = v ./ sum(v);
    s = doc_subreg{doc}(sampleDiscrete(v,1));   % Pick the subregion
    
    % Only components in current subgroup
    p_vec = m_c .* (clstr_subrg == s);

    tail = ismember([1:num_subreg],doc_subreg{doc});
    tail = tail .* measure([1:num_subreg]);

    p_vec = [p_vec, tail];

    sel_clst = sampleDiscrete(p_vec./sum(p_vec),1);

    % If the selected cluster is a new cluster

    if(sel_clst > num_clstrs)
      num_clstrs = num_clstrs + 1;
      sel_clstr_corr = num_clstrs;
      sel_subreg     = sel_clst - (num_clstrs - 1);

      n_c(doc,sel_clstr_corr) = 1;
      m_c(sel_clstr_corr)     = 1;
      theta(sel_clstr_corr)   = var_theta*randn(1,1);

      clstr_subrg(sel_clstr_corr) = sel_subreg;
      clstr_mean(sel_clstr_corr)  = 0;
    else
      sel_clstr_corr = sel_clst;
      sel_subreg     = clstr_subrg(sel_clst);

      n_c(doc,sel_clstr_corr) = n_c(doc,sel_clstr_corr) + 1;
      clstr_mean(sel_clstr_corr) = clstr_mean(sel_clstr_corr) * m_c(sel_clstr_corr);
      m_c(sel_clstr_corr)     = m_c(sel_clstr_corr) + 1;
    end;

    c{doc}(i) = sel_clstr_corr;
    X{doc}(i) = theta(sel_clstr_corr) + var_x*randn(1,1);
    clstr_mean(sel_clstr_corr) = (clstr_mean(sel_clstr_corr) + X{doc}(i)) / ...
                                  m_c(sel_clstr_corr);
  end;
end;
