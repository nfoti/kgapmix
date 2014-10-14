

means = [-5 0 5];
taus = [1/.2 1 1/.2];
p = [.4 .2 .4];

Z = discreternd(N, p);
Y = randn(1,N)./sqrt(taus(Z)) + means(Z);