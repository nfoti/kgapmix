function h = gen_pred_waterfall_slice(fn,Y,X,times,xstr,ystr)

load(fn);

ntime = max(XXind);

N = size(Y);

h = zeros(1,ntime+1);

% Plot the data points on the z==0 plane
h(1) = plot3(X,Y,zeros(N),'ro');
hold on;

% Surface just doesn't look good yet
% [X Y] = meshgrid(times,predvals);
% Z = zeros(size(X));
% for i = 1:ntime
%   Z(:,i) = pl_slice(XXind==i);
% end
% h(2) = surf(X,Y,Z);

% Now plot the densities
for i = 1:ntime
%  h(i+1) = plot3(times(i).*ones(size(predvals)), predvals, ...
%                 pl_slice(XXind==i)', 'LineWidth',2);
  tt = [times(i).*ones(size(predvals)); predvals; pl_slice(XXind==i)'];
  tf = [1:size(tt,2) size(tt,2):-1:1];
  patch('vertices', tt', 'faces', tf, 'facecolor', 'none', 'edgecolor', 'b', 'edgealpha', 0.5, 'linewidth', 2);
end

set(gca,'FontSize',18);
xlabel(xstr);
ylabel(ystr);
zlabel('Density');

%uistack(h(1),'bottom');