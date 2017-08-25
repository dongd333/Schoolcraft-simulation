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

minx=min(CT(:,3));
maxx=max(CT(:,3));
miny=min(CT(:,4));
maxy=max(CT(:,4));

nlags=25; tolfact=.5
sep=linspace(0,nlags-1,nlags);
cov=zeros(nlags,39);
num=zeros(nlags,39);
% each layer first 
for lay=1:39
    meanX=0;
  work=CT(:,3:5);
  work(:,3)=CT(:,lay+9);
  apple=find(work(:,3)>0);
  work = work(apple,:);
  nloc=size(apple,1);
  lag=zeros(nloc);
  
  for i=1:nloc
%      posnow=[work(i,1) work(i,2)];
      blah1=work(:,1)-work(i,1);
      blah2=work(:,2)-work(i,2);
      lag(i,:)=(blah1.^2+blah2.^2).^0.5;
  end
  
    X=work(:,3);
    meanX=meanX+mean(X);
    X=X-mean(X);
    X=X*X';

    for l=1:nlags
        use=( lag<sep(l)+tolfact & lag>sep(l)-tolfact );
%        use=( lag<sep(l)*(1+tolfact) & lag>sep(l)/(1+tolfact) );
        cov(l,lay)=cov(l,lay)+sum(sum(X(use)));
        num(l,lay)=num(l,lay)+sum(sum(use));
    end
  cov(:,lay)=cov(:,lay)./(num(:,lay)-1);
  %cov(1,lay)=var(work(:,3));
  figure(1);
  plot(sep,cov(:,lay));
  %axis([0 10 -10^-6 10^-6])
  axis square
  hold on
  figure(2);
  plot(sep,cov(:,lay)/cov(1,lay));
  %axis([0 10 -10^-6 10^-6])
  axis square
  hold on
  
end
% get thickness-weighted average
dz=abs(diff(top));
dz=dz(6:44,2);
avecov=0*cov(:,1); avecorr=avecov;
for kk=1:nlags 
   avecov(kk)=cov(kk,:)*dz./sum(dz);   
   avecorr(kk)=avecov(kk)/avecov(1);
end
figure(1),(plot(sep,avecov,'k'));
figure(2),(plot(sep,avecorr,'k'));