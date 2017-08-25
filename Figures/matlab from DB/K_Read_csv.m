% This file reads data from csv file to matrix
clear all;
M = csvread('K.csv',1,0);

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

% The rest columns are K value, layer 1 to layer 44 (mc-4)
% Because the dimension Matlab defines 3d matrix is
% [size(y),size(x),size(z)], the order of x and y should be changed.
K=zeros(max(R),max(C),mc-4);

for i=1:line
    % Keep unique value of x and y
    x0(Column(i))=X(i);
    y0(Row(i))=Y(i);
    for j=1:mc-4
        K(Row(i),Column(i),j)=M(i,j+4);
    end
end

K=log10(K);
% Find the maximum and minimum values of K for color display
%[m1,id1]=max(K);
max(max(max(K)));
min(min(min(K)));

%k1(:,:)=K(:,:,5);   % A test for surf plot
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
size(K);
%Surf(x0,y0,k1);
%caxis([0,75]);

%figure(2)
%hs=slice(K,x0,y0,z0);          % This plot the 4d surface as a block
%set(hs,'EdgeColor','none');
set(gcf,'papersize',[11 8.5])
% Define the cutting surface
sx=[40,97.7];
sy=[30,56.2];
sz=[0.5];
hs=slice(x0,y0,z0,K,sx,sy,sz);

%cm=flipud(gray);
%cm=bone;
cm=gray;
colormap(cm);

set(hs,'EdgeColor','none');     % This removes the edge color
caxis([0,1.9]);                % Define the range of color display
%opengl software                 % This used to enhance the display of bar
h2 = colorbar;
%h2.Label.String='K (m/day)';
%h2.Label.Interpreter='tex';

ylabel(h2, '$log_{10}(K)$ (m/day)','Interpreter','LaTex');     % Label the color bar

%set(h2, 'ylim', [0 1.9],'position',[.95 .2 .05 .7]);        % Define the range of bar to show
           
hold on;
xlabel('$x$ (m)','Interpreter','LaTex');
ylabel('$y$ (m)','Interpreter','LaTex');
zlabel('$z$ (m)','Interpreter','LaTex');
axis equal
box on
print('-painters','-dpdf','K_map_DB.pdf')




%axis([0 97.75 0 56.25 0 27.5]);
