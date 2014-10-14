function clstrs(theta,m_c,n_c,th,num_clstrs)

clf
clr = colormap(lines);

% What type of bar plot to display
stacked = 0;

num_bns = 100;

mn  = min([theta min(th{num_clstrs,1})]);
mx  = max([theta max(th{num_clstrs,1})]);

theta = [mn theta mx];

[num_docs, tmp] = size(n_c);
if(stacked == 0)
  subplot(num_docs+1,1,1);
  hold on;

  count   = [zeros(num_docs,1)  n_c  zeros(num_docs,1)];
  for i = [1:num_docs]
    bar([theta],count(i,:),1,'Facecolor',clr(i,:))
  end;

  for i = [1:num_docs-1]
    bar([theta],count(i,:)  .* (count(i,:) & count(i+1,:)),1,'Facecolor',clr(num_docs+i,:))
    bar([theta],count(i+1,:).* (count(i,:) & count(i+1,:)),1,'Facecolor',clr(num_docs+i,:))
  end;

else
  count   = [0 m_c 0];
end;

[t_hist(1,:),bn]   = theta_hist(theta,count,num_bns);
for i = [1:num_docs]
  [t_hist(i+1,:),bn] = theta_hist(th{num_clstrs,1},th{num_clstrs,3}(:,i,:),num_bns);
end;

if stacked == 0
  for i = 1:num_docs
    subplot(num_docs+1,1,1+i);
    bar(bn,t_hist(i+1,:),.5,'Facecolor', clr(i,:))
  end;

else
  colormap('default');
  t_hist = t_hist./repmat(sum(t_hist,2),1,num_bns);
  bar(bn,t_hist',.2,'stacked')
  legend({num2str([0:num_docs]')});
end;
