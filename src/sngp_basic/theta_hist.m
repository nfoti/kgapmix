function [t_h, bin] = theta_hist(th,counts,num_bins)


%%%%%%%%%%%%%%%%%%%
%function [t_h, bin, doc_no] = theta_hist(post_samples,num_clstr,num_bins, doc_no)

% count  = 2     % count = 1 => counts across documents     (mc)
%                % count = 2 => counts in document "doc_no" (nc)

% th = post_samples{num_clstr,1}(:);

% if(count == 1)
%   mc = post_samples{num_clstr,2}(:);
% else
%   mc = post_samples{num_clstr,3}(:,doc_no,:);
%   mc = mc(:);
% end;
%%%%%%%%%%%%%
  
  mc = counts(:);
  th = th(:);

  th_min = min(th);
  th_max = max(th);

  bin_wdth = (th_max - th_min)/(num_bins + 1);

  for i = [0:num_bins-1]
    bin(i+1) = sum(th_min + i*bin_wdth);
    j = find(th >= (bin(i+1)) & th <= (bin(i+1) + bin_wdth));
    t_h(1,i+1) = sum(mc(j));
  end;
  bin = bin + bin_wdth/2;

%  bar(bin,t_h,'g')

%xy = 5;
%bar(theta,m_c,'b')
%hist([X{1},X{2}],100)
%r(th{17,1}(xy,:),th{17,2}(xy,:),5,'r')
