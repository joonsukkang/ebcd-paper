%%%%%%%%%%%%%%%%%%% using code downloaded from 
% http://www.montefiore.ulg.ac.be/~journee/GPower.zip

clear all;

rng('default');
rng(1);

m=2;
mu=[1, 1/2];
n = 50;
p = 500;

X = readmatrix(strcat('../data/sim1.csv'));
seed_max = size(X,1)/n;


matr = [];
tic
for gg = (0.01:0.01:0.99)

    gamma = gg*ones(1,m);
    matrgg = zeros(m*seed_max, p);

    for seed = 1:seed_max

        nrow_from  = (seed - 1) * n +1;
        nrow_to = seed * n;
        A = X(nrow_from:nrow_to,:);
        
        Z1 = GPower(A,gamma,m,'l1',1,mu);     % Block sparse PCA, l1 penalty
        L = Z1 * diag(mu);

        matrgg(((seed-1)*m + 1): (seed*m), :) = transpose(L);

    end

    matr = [matr; matrgg];

end
toc/seed_max

writematrix(matr, strcat('../output/gpower_sim1.csv'));     
exit();