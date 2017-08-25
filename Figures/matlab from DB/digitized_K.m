K_d=load('digitized_K.dat','-ascii');
lnK=log(K_d(:,1));
var(lnK)
var(exp(lnK))