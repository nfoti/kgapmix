% Figure for synthetic data
%
% Contains i.   plot of the data
%          ii.  hist of log(L)
%          iii. trace plot of number of used cluster for each sampler

addpath ../util;

figure;

% Generate the fake data used
load nips_seeds.mat;
s = RandStream('mt19937ar','Seed',seeds(1));
RandStream.setGlobalStream(s);

gen_fake_data_box;

subplot(1,3,1);
set(gca,'FontSize',18);

for i = 1:numel(Ycell)
  plot((i*.1).*ones(size(Ycell{i}')), Ycell{i}', 'o');
  hold on;
end
set(gca,'XLim',[0 1.1]);
set(gca,'XTick',(0:.2:1.1));
set(gca,'XTickLabel',num2str((0:.2:1.1)'));
xlabel('Time');
ylabel('Observation');


load synth_box_slice/synth_box_slice_1.mat;
KKslice = KK;
clearvars -except Ls KKslice;

load synth_box_sngp/synth_box_sngp_1.mat;
KKsngp = KK;
clearvars -except Ls KKslice KKsngp;

load synth_box_finite/synth_box_finite_1.mat;
KKfinite = KK;
clearvars -except Ls KKslice KKsngp KKfinite;

subplot(1,3,2);
set(gca,'FontSize',18);
hslice = plot(KKslice, 'b', 'LineWidth', 2);
hold on;
hsngp = plot(KKsngp, 'r', 'LineWidth', 2);
hfinite = plot(KKfinite, 'g', 'LineWidth', 2);
set(gca,'XLim', [0,numel(KKslice)]);
xlabel('Iteration');
ylabel('# used clusters');

legend([hslice, hsngp, hfinite], 'Slice', 'SNGP', 'Finite');

subplot(1,3,3);
set(gca,'FontSize',18);
hist(log(Ls),25);
xlabel('log(L)');
ylabel('Frequency');

rmpath ../util;