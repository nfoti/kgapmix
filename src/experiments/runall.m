
sdir = regexp(pwd, '/', 'split');
if ~strcmp(sdir{end},'experiments')
  fprintf('Run this script in experiments/ directory\n');
  return;
end

fprintf('RUNNING CMB\n');
run_cmb_box_finite;
run_cmb_box_slice;
run_cmb_box_sngp;
run_cmb_se_finite;
run_cmb_se_slice;

fprintf('RUNNING MOTOR\n');
run_motor_box_finite;
run_motor_box_slice;
run_motor_se_finite;
run_motor_se_slice;

fprintf('RUNNING SYNTH\n');
run_synth_box_finite;
run_synth_box_slice;
run_synth_box_sngp;
