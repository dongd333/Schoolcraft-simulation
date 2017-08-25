% This file reads data from csv file to matrix
clear all;
M = csvread('Initial Conc.csv',1,0);

size(M(1,:));
% Total line of the matrix
line=size(M,1);
% Total column of the matrix
mc=size(M,2);

% The first column is Column of the model
Column=zeros(line,1);
% The second column is Row of the model
Row=zeros(line,1);
% The third column is x axis
X=zeros(line,1);
% The fourth column is y axis
Y=zeros(line,1);


Column(:)=M(:,1);
Row(:)=M(:,2);
X(:)=M(:,3);
Y(:)=M(:,4);

% Keep unique value of column and row
C=unique(Column);
R=unique(Row);

% Maximum number of column and row
max(C);
max(R);
x0=zeros(max(C),1);
y0=zeros(max(R),1);

% The rest columns are initial concentrations, layer 1 to 44 (mc-4)
% Because the dimension Matlab defines 3d matrix is
% [size(y),size(x),size(z)], the order of x and y should be changed.
IC=zeros(max(R),max(C),mc-4);

for i=1:line
    % Keep unique value of x and y
    x0(Column(i))=X(i);
    y0(Row(i))=Y(i);
    for j=1:mc-4
        IC(Row(i),Column(i),j)=M(i,j+4);
    end
end
%IC=log10(1000*IC);
IC=1000*IC;
IC(IC<-10)=-10;
% Find the maximum and minimum values of IC for color display
%[m1,id1]=max(IC);
maxC=max(max(max(IC)))
minC=min(min(min(IC)))

%IC1(:,:)=IC(:,:,5);   % A test for surf plot
z0=zeros(44,1);
%z0(:,1) = linspace(43,0,44);

% Read in Z elevation value; dlmread does not work for txt file
z1=dlmread('Z.dat');
%for i=1:size(z1,1)-1
%    z0(i)=(z1(i)+z1(i+1))/2;
%end
z0(:)=z1(1:44);

% Check the size for plot
size(x0);
size(y0);
size(z0);
size(IC);
%Surf(x0,y0,IC1);
%caxis([0,75]);

%figure(2)
%hs=slice(IC,x0,y0,z0);          % This plot the 4d surface as a block
%set(hs,'EdgeColor','none');

% Define the cutting surface
sx=[40,97.7];
sy=[30,56.2];
sz=0.5;
hs=slice(x0,y0,z0,IC,sx,sy,sz);
colormap(jet)

set(hs,'EdgeColor','none');     % This removes the edge color
axis equal
caxis([0 maxC]);                % Define the range of color display
% opengl software                 % This used to enhance the display of bar
shading flat
h2 = colorbar;
set(h2, 'ylim', [0 maxC]);        % Define the range of bar to show
ylabel(h2, 'Concentration ($ppb$)','Interpreter','LaTex');                % Label the color bar
%hold on;
xlabel('$x$ (m)','Interpreter','LaTex');
ylabel('$y$ (m)','Interpreter','LaTex');
zlabel('$z$ (m)','Interpreter','LaTex');
light('Position',[-1 0 0],'Style','local')
%axis([0 97.75 0 56.25 0 27.5]);
print('-painters','-dpdf','blah.pdf')
