

load ~/work/data/cmb/cmb_MMTT.mat;

MM = double(MM(1:600));
TT = TT(1:600);

MM = MM-min(MM);
MM = MM./max(MM);
TT = (TT-mean(TT))./std(TT);

times = 0.1*(.5:.5:9)';
