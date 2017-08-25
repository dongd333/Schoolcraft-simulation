%% --------------------------------------------------------------
clear all
% Read in extracted RW3D files
nf = 50;            % Number of files
ngrid = 40;        % ngrid defined in rt3d input
spe = 14;         % Number of species
nwell = 25;         % Number of wells

c0_Nit=42;          % Initial nitrate concentration
tcutoff = 163.80;  % Cutoff time

t = linspace(0.105555,tcutoff,ngrid);

% Read the factor for the numbers to be divided by
factor(:,:) = csvread('Factor.csv');
% Read in the assembled breakthrough curve from RW3D
for f = 1:nf
    name = sprintf('C_Nitrate_%03d.csv', f+0);
    % S has the dimension of (ngrid, nwell, nfile)
    S(:,:,f) = csvread(name);
    
end
% S_ave = mean(S,3) takes the average of the nf files
% S(:,1,1) gets the data for the first well at different time


%S0 = mean(S,3)./factor./c0_Nit;  % Average c/c0 to be plotted

% The 8th row is day 29.51
S(8,:,:);

D_30(:,:) = S(8,:,:);
% The 18th row is day 71.51
D_72(:,:) = S(18,:,:);
% The 30th row is day 121.91
D_122(:,:) = S(30,:,:);

% Calculate the moving (running) average
% The function only works for new matlab. For old version, just use mean
for n=2:nf
    m_d30(:,n)=mean(D_30(:,1:n),2);
    m_d72(:,n)=mean(D_72(:,1:n),2);
    m_d122(:,n)=mean(D_122(:,1:n),2);
end

% The average of the first value is the first value
m_d30(:,1)=D_30(:,1);
m_d72(:,1)=D_72(:,1);
m_d122(:,1)=D_122(:,1);
% Calculate the difference from the final average
% matrix minus a vector need the function diag()*ones
m_d30_p = m_d30 - diag(m_d30(:,nf))*ones(size(m_d30));
m_d72_p = m_d72 - diag(m_d72(:,nf))*ones(size(m_d72));
m_d122_p = m_d122 - diag(m_d122(:,nf))*ones(size(m_d122));

% The average values
m_d30(:,nf)
m_d72(:,nf)
m_d122(:,nf)

runs = linspace(1,nf,nf);
col =jet(nwell);

subplot(3,1,1); 
for i = 1:nwell
    plot(runs, m_d30_p(i,:)/m_d30(i,nf),'color', col(i,:));
    hold on;
end
ylim([-0.2 0.2])
%ylabel('Day 30') 
subplot(3,1,2); 
for i = 1:nwell
    plot(runs, m_d72_p(i,:)/m_d72(i,nf),'color', col(i,:));
    hold on;
end
%ylabel('Day 72'); 
ylim([-0.2 0.2])

ylabel('$(c_i - c_{ave})/c_{ave}$','interpreter','latex','FontName','Times New Roman','FontSize',12)
subplot(3,1,3); 
for i = 1:nwell
    plot(runs, m_d122_p(i,:)/m_d122(i,nf),'color', col(i,:));
    hold on;
end
xlabel('Realizations','interpreter','latex','FontName','Times New Roman','FontSize',12)

samexaxis('abc','xmt','on','ytac','join','yld',1)

%% --------------------------------------------------------------