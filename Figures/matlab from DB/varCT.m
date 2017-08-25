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
numd=size(CT,1);
for lay=1:39
  blahv( (lay-1)*numd + 1:lay*numd, 1 )=CT (:,lay+9);
end
var(blahv)
mean(blahv)