% In order to plot the legend clearly, load all the data and plot together


%% Nitrate
figure(1)
%% Measurements
Nitrate2 = csvread('Nitrate Day 2.csv',0,0);
h1(1)=plot(Nitrate2(:,1),Nitrate2(:,2),'^','markerfacecolor','b','markersize',9);
hold on
Nitrate26 = csvread('Nitrate Day 26.csv',0,0);
h1(2)=plot(Nitrate26(:,1),Nitrate26(:,2),'ro','markerfacecolor','r','markersize',9);
hold on

%% Plot PT results
omega = 100;                 % Column length [cm]
ncount = 100;                % Counting interval for particle numbers as a funNitrateion of length 
bin =1;                     % New bin size, input is 1 cm
nbin = ncount/bin;          % Number of array with new bin size

Nitrate_d2_m = zeros(nbin,1); 
Nitrate_d2_sd = zeros(nbin,1);
Nitrate_d26_m = zeros(nbin,1); 
Nitrate_d26_sd = zeros(nbin,1);

numxb0 = linspace(bin/2,omega-bin/2,nbin);
size(numxb0)

% Load the results from simulations

load S_13_N2_1000.dat
Nitrate_d2_sd(:,1) = S_13_N2_1000;

load S_13_N26_1000.dat
Nitrate_d26_sd(:,1) = S_13_N26_1000;

load M_13_N2_1001.dat
Nitrate_d2_m(:,1) = M_13_N2_1001;

load M_13_N26_1001.dat
Nitrate_d26_m(:,1) = M_13_N26_1001;

% Plot the average
errorbar(numxb0(1:3:end),Nitrate_d2_m(1:3:end),Nitrate_d2_sd(1:3:end)/2,'+','LineWidth',1);
hold on
h1(7)=plot(numxb0(1:3:end),Nitrate_d2_m(1:3:end),'-+','LineWidth',2);
hold on
errorbar(numxb0(1:3:end),Nitrate_d26_m(1:3:end),Nitrate_d26_sd(1:3:end)/2,'rx','LineWidth',1);
hold on
h1(8)=plot(numxb0(1:3:end),Nitrate_d26_m(1:3:end),'r-x','LineWidth',2);
hold on


%% Legend
xlabel('Length ($cm$)','Interpreter','LaTex')
ylabel('Concentration ($mg/L$)','Interpreter','LaTex')

axis([0 100 0 30]);
% Legend at axes 1
ah1 = gca;
CL1=legend(ah1,h1(1:2),'Measurement Day 2','Measurement Day 26',1,'Location','southoutside','Orientation','horizontal');
set(CL1,'color','none');
%legend boxoff 
% Second column
% Axes handle 2 (unvisible, only for place the second legend)
ah2=axes('position',get(gca,'position'), 'visible','off');
CL2=legend(ah2,h1(7:8),'PT Simulation Day 2','PT Simulation Day 26',4,'Location','southoutside','Orientation','horizontal');
set(CL2,'color','none');
%legend boxoff 

