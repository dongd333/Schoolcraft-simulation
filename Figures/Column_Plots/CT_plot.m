% In order to plot the legend clearly, load all the data and plot together


%% CT
figure(1)
%% Measurements
CT2 = csvread('CT_Day 2.csv',1,0);
h1(1)=plot(CT2(:,1),CT2(:,2)/1000,'s','markerfacecolor','b','markersize',9);
hold on
CT26 = csvread('CT_Day 26.csv',1,0);
h1(2)=plot(CT26(:,1),CT26(:,2)/1000,'rd','markerfacecolor','r','markersize',9);
hold on

%% Plot PT results
omega = 100;                 % Column length [cm]
ncount = 101;                % Counting interval for particle numbers as a function of length 

CT_d2_m = zeros(ncount-1,1); 
CT_d2_sd = zeros(ncount-1,1);
CT_d26_m = zeros(ncount-1,1); 
CT_d26_sd = zeros(ncount-1,1); 

CT_d2_m2 = zeros(ncount-1,1); 

numx0 = linspace(0,omega,ncount);
numxp = (numx0(1:end-1)+numx0(2:end))/2;
size(numxp);

% Load the results from simulations
load M_CT26_1000.dat
load S_01_13_CT26_600.dat
CT_d26_m(:,1) = M_CT26_1000;
CT_d26_sd(:,1) = S_01_13_CT26_600;

load M_CT2_1001.dat
load S_01_13_CT2_600.dat
CT_d2_m(:,1) = M_CT2_1001;
CT_d2_sd(:,1) = S_01_13_CT2_600;


% Plot the average
errorbar(numxp(1:3:end),CT_d2_m(1:3:end),CT_d2_sd(1:3:end),'+','LineWidth',1);
hold on
h1(7)=plot(numxp(1:3:end),CT_d2_m(1:3:end),'-','LineWidth',2);
hold on
errorbar(numxp(1:3:end),CT_d26_m(1:3:end),CT_d26_sd(1:3:end),'rx','LineWidth',1);
hold on
h1(8)=plot(numxp(1:3:end),CT_d26_m(1:3:end),'r-','LineWidth',2);
hold on

%% Legend
xlabel('Length ($cm$)','Interpreter','LaTex')
ylabel('Concentration ($mg/L$)','Interpreter','LaTex')

axis([0 100 0 0.12]);
% Legend at axes 1
ah1 = gca;
CL1=legend(ah1,h1(1:2),'Measurement Day 2','Measurement Day 26',1,'Location','southoutside','Orientation','horizontal');
set(CL1,'color','none');
% Second column
% Axes handle 2 (unvisible, only for place the second legend)
ah2=axes('position',get(gca,'position'), 'visible','off');
% Fourth column
% Axes handle 2 (unvisible, only for place the second legend)
ah4=axes('position',get(gca,'position'), 'visible','off');
% Legend at axes 2
CL4=legend(ah4,h1(7:8),'PT Simulation Day 2','PT Simulation Day 26',4,'Location','southoutside','Orientation','horizontal');
set(CL4,'color','none');
