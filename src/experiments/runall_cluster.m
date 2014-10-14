
sdir = regexp(pwd, '/', 'split');
if ~strcmp(sdir{end},'experiments')
  fprintf('Run this script in experiments/ directory\n');
  return;
end

cluster = parcluster('anthill');

qsubargs = '-cwd -V -l h_rt=23:59:59 -l virtual_free=2G'
set(cluster, 'IndependentSubmitFcn', {@independentSubmitFcn, qsubargs});

njobs = 12;
jobs = cell(1,njobs);

%fprintf('RUNNING CMB\n');
%batch(cluster, 'run_cmb_box_finite');
batch(cluster, 'run_cmb_box_slice');
%batch(cluster, 'run_cmb_box_sngp');
%batch(cluster, 'run_cmb_se_finite');
batch(cluster, 'run_cmb_se_slice');

fprintf('RUNNING MOTOR\n');
%batch(cluster, 'run_motor_box_finite');
batch(cluster, 'run_motor_box_slice');
%batch(cluster, 'run_motor_se_finite');
batch(cluster, 'run_motor_se_slice');

fprintf('RUNNING SYNTH\n');
%batch(cluster, 'run_synth_box_finite');
batch(cluster, 'run_synth_box_slice');
%batch(cluster, 'run_synth_box_sngp');

fprintf('Remember to clean up the jobs files when theyre all done');
