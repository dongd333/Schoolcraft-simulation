clear all
top=load('Layer_elevation.dat');
CT=load('CT.csv');
nlay=size(top,1)-1;  %44
z=top(:,2);
for i=2:nlay+1
    z(i-1)=0.5*(z(i-1)+z(i));
end
z(nlay+1)=[];
% Ignore first 5 layers which have CT=0 
z=z(6:nlay); nlay=nlay-5

for i=1:nlay
    for j=1:nlay
        lag(i,j)=z(i)-z(j);
    end
end
lag=abs(lag);

nlags=41; tolfact=0.26
sep=linspace(0,10,nlags);

locs=size(CT,1);
%locs=1000
cov=0*sep;num=cov;meanX=0;
for k=1:locs
    X=CT(k,10:48);
    meanX=meanX+mean(X);
    X=X-mean(X);
    X=X'*X;
    for l=1:nlags
        use=( lag<sep(l)+tolfact & lag>sep(l)-tolfact );
%        use=( lag<sep(l)*(1+tolfact) & lag>sep(l)/(1+tolfact) );
        cov(l)=cov(l)+sum(sum(X(use)));
        num(l)=num(l)+sum(sum(use));
    end
end
cov=cov./(num-1);
figure(1);
plot(sep,cov,'o');
axis([0 10 -0.0001 0.0001])
axis square
meanX/locs