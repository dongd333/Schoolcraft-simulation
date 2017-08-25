function [] = plot_nitrate_smooth()
close all; clear; clc;
%-------------------------------------------------------
c0_Nit=42;

tiny = 1.0e-06;

xmin=0; xmax=163.905555; ymin=-0.1; ymax=2.0; 

xt=[0 75 150];yt=[0 0.5 1.0 1.5 2.0];
tlx=0.025; tly=0.025;
pos1=get(0,'screensize');
pos=[pos1(1)+10 pos1(2)+40 pos1(3)-10 pos1(4)-120];

xtext=17;ytext=1.1;
tcutoff = 163.80;  % Cutoff time

nw   = 5; %number of obs. wells
nz   = 5; %number of screen for each well

m=5;  % Rows
n=5;  % Columns

delx = (pos1(3) - pos1(1)-20)/nz;
dely = (pos1(4) - pos1(2)-160)/nw;
%% load workspace containing measured data:
Nitrate = load ('measured_time_conc_nitrate_ant_from_tracer_start.mat');

for p=1:nw
    for q=1:nz
        Nitrate_meas{p,q}=Nitrate.(sprintf('c_%d_%d_m', p+8,q*10+25));
        t_Nitrate_m{p,q}=Nitrate.(sprintf('t_%d_%d_m', p+8,q*10+25));
    end
end
% Nitrate_meas reads the measured Nitrate concentration saved in the mat file
% t_Nitrate_m reads the time associated with measured Nitrate concentrations

%% RT3D Simulation
file2='RT3D005.OBS';  %Tracer Simulation, Bromide

%load the RT3D OBS file:
[g]=read_obs_file('RT3D004.OBS',26,0);   %Nitrate

%% ------------------------------------------------------------------------
% RW3D results
nf = 5;            % Number of files
ngrid = 40;        % ngrid defined in rw3d input
spe = 14;         % Number of species

t = linspace(0.105555,tcutoff,ngrid);

% Read in the assembled breakthrough curve from RW3D
    S1(:,:) = csvread('C_Nitrate_S1.csv');
    SD1(:,:) = csvread('C_Nitrate_std_S1.csv');
    

%% -------------------------------------------------------------------------
h1=figure(1); 
for ii=1:nw
    for jj=1:nz
        kk = (ii-1)*nz + jj;
        subplot(m,n,kk);
        
        % Plot measurements
        t_Nitrate_m{ii,jj} = t_Nitrate_m{ii,jj}(t_Nitrate_m{ii,jj}>0);
        Nitrate_meas{ii,jj} = Nitrate_meas{ii,jj}(t_Nitrate_m{ii,jj}>0);
        
        find(t_Nitrate_m{ii,jj} >= tcutoff);
        ind4= min(find(t_Nitrate_m{ii,jj} >= tcutoff));
        
        plot(t_Nitrate_m{ii,jj}(1:ind4),Nitrate_meas{ii,jj}(1:ind4)./c0_Nit,'bo','markersize',6,'markerfacecolor','c');
        hold on;
        
        % Check and correct the negative values in time
        %for rr=1:ind4
        %    if t_Nitrate_m{ii,jj}(rr) < 0
        %        t_Nitrate_m{ii,jj}(rr) = 0;
        %    end
        %end
        
        % Restrict the data within the same range as simulation
        t_Nitrate_m{ii,jj}(ind4) = tcutoff;
        Nitrate_meas{ii,jj}(ind4) = 0;
        
        % Plot RW3D outputs
        
        plot(t(:), S1(:,kk), 'g-v','markersize',3);   % Plot non-zeros
        hold on;
        errorbar(t(:), S1(:,kk),SD1(:,kk),'g+','LineWidth',1);
        hold on
            
        % Plot RT3D outputs
        ind3=min(find(g(2,:)' >= tcutoff));
        plot(g(2,1:ind3)',g(kk+2,1:ind3)'./c0_Nit,'b-');
        hold on;
        
        % Calculate the root mean square errors
        % Firstly, the closest data with the same time of measurement is
        % obtained. Using "linear" method, can be changed to linear or
        % other
        
        % Re-initialize all array for every loop
        Nitrate_m_v =[];
        t_m_v = [];
        Nitrate_ws = [];
        t_ws = [];
        
        Nitrate_m_v = Nitrate_meas{ii,jj}(1:ind4)./c0_Nit;
        t_m_v = t_Nitrate_m{ii,jj}(1:ind4);
        n_w = ngrid;
        
        t_m_v(isnan(Nitrate_m_v)) = [];
        Nitrate_m_v(isnan(Nitrate_m_v)) = [];
        
        Nitrate_ws(2:n_w+1) = S1(:,kk);
        Nitrate_ws(1) = 1;               % Add two data to the array (top and end)
        %Nitrate_ws(n_w+2) =0;
        t_ws(2:n_w+1) = t(:);
        t_ws(1) = 0;
        
        %if any(t_ws)
        %    t_ws(n_w+2) = 163.905555;
            % Interpolated array on the measured time
            ip_w = interp1(t_ws,Nitrate_ws,t_m_v);
            % Root mean square error
            rms_w=sqrt(sum((ip_w-Nitrate_m_v).^2)/length(Nitrate_m_v));
        %else
        %    rms_w = sqrt(sum(Nitrate_m_v.^2)/length(Nitrate_m_v));
        %end
        %t_ws
        %ip_w
        %t_m_v
        %Nitrate_m_v
        
        % RT3D outputs       
        n_t = ind3;
        Nitrate_ts = [];
        t_ts = [];
        Nitrate_ts(2:n_t+1) = g(kk+2,1:ind3)'./c0_Nit;
        Nitrate_ts(1) = 1;               % Add two data to the array (top and end)
        t_ts(2:n_t+1) = g(2,1:ind3)';
        t_ts(1) = t_Nitrate_m{ii,jj}(1);
        
        ip_t = interp1(t_ts,Nitrate_ts,t_m_v);
        %t_ts
        %t_m_v
        %size(Nitrate_ts)
        %size(t_m_v)
        %length(Nitrate_m_v)
        % Root mean square error for rt3d
        rms_t=sqrt(sum((ip_t-Nitrate_m_v).^2)/length(Nitrate_m_v));
        

        titletext=strcat(num2str(ii+8),'-', num2str(jj*10+25));     % Well starts from 9 and screen from 35
        axis([xmin xmax ymin ymax]);
        %text(xtext,ytext,titletext);
        
        set(gca,'fontsize',14,'fontname','times','FontWeight','normal','LineWidth',1.1,'TickLength',[tlx,tly]);
        set(gca,'xtick',xt,'ytick',yt);
        set(gca, 'xlim', [xmin xmax], 'xtick',[-50 0 50 100 150 200]);
        
        %Fix digit numbers
        rms_w=round(rms_w*1e2)/1e2;
        rms_t=round(rms_t*1e2)/1e2;
        hold on;
        text(xtext+100,ytext+0.5,{strcat('RW3D: ', num2str(rms_w)),strcat('RT3D: ',num2str(rms_t))}, 'FontSize', 8);

    end
end
 
jointfig(h1,m,n);
set(h1,'position', pos);
set(gcf,'PaperPositionMode','auto');

hold on
AxesH = axes('Parent', h1, ...
  'Units', 'normalized', ...
  'Position', [0, 0, 1, 1], ...
  'Visible', 'off', ...
  'XLim', [0, 1], ...
  'YLim', [0, 1], ...
  'NextPlot', 'add');

text(0.5,0.11, 'Time (Days)','Interpreter','LaTex','fontname','times','FontSize',14);
text(0.12,0.5,'$C/C_0$','Interpreter','LaTex','rotation', 90, 'fontname','times','FontSize',14) % y-axis label

end

        
%--------------------------------------------

%pubmode('on');

%-------------------------------END OF PLOT_CONC2--------------------------------
function y = nansum(x)
nans = isnan(x);
i = find(nans);
x(i) = zeros(size(i));
% Protect against an entire column of NaNs
y = sum(x);
i = find(all(nans));
y(i) = i + NaN;
end



%----------------------------------------------------------------------------------
function y = nanabs(x)
%abs ignoring nan's; written by Phani; 16 April 2001
if nargin < 1 
    return
end
if isnan(x)
   y=0;
else
    y=abs(x);
end
end
%------------------------------------------------------------------
function jointfig(hfig,no_row,no_col)
 %   jointfig(hfig,no_row,no_col)
%	--> joint subplots without any space between them
%   hfig : figure handler, if none, keyin gcf instead
%   no_row    : No. of row subplots
%   no_col    : No. of column subplots % All the movement of subplots should be done in unit of points
figure(hfig), hsubplot = get(hfig,'Children');
% convert the position unit from pixel into points : should be restored)
set(hfig,'unit','point');
 % BEWARE! hsubplot has different order from the original subplot sequence
% for instance,
%
%  -----------------------         -----------------------
%  |     1    |     2     |        |     4    |     3     |
%  |----------+-----------|        |----------+-----------|
%  |     3    |     4     |        |     2    |     1     |
%  -----------------------         -----------------------
%       subplot(22i)                  get(gcf,'Children')
%
% THEREFORE, transpose hsubplot into the one in original subplot sequence, like this..
 hsubplot = hsubplot(length(hsubplot):-1:1);
 no_subplot1 = length(hsubplot);
no_space  = no_row*no_col;
no_delta = no_space - no_subplot1;
 % in case of the odd number of subplots
if no_delta,
	for i = 1:no_delta
		addsubplot = subplot(no_row,no_col,no_subplot1+i);
		hsubplot = [hsubplot; addsubplot];
	end
end
 no_subplot = length(hsubplot);
 % Default position of figure in a window in point coord
 for i=1:no_subplot,
	set(hsubplot(i),'unit','point'),
	tmp_ylab_pos = get(get(hsubplot(i),'ylabel'),'position');
	ylab_pos(i) = tmp_ylab_pos(1);
end
 new_ylab_pos = min(ylab_pos);  
 coner1 = get(hsubplot(1),'position');
 coner2 = get(hsubplot(length(hsubplot)),'position');
 % position of lowest-left coner
inix = coner1(1);
iniy = coner2(2)*1.13;
 % axis line width
alinewidth = get(hsubplot(1),'linewidth');
 % total lengths
total_xlength = (coner2(1) + coner2(3) - coner1(1)) + (no_col-1) * alinewidth;
total_ylength = (coner1(2) + coner1(4) - coner2(2)) + (no_row-1) * alinewidth;
 % width of each subplot
 delx = 1.0 * total_xlength / no_col;  
 dely = 0.97 * total_ylength / no_row;  
 % height of each subplot dely = 0.97 * total_ylength / no_row;  
 %index_loop = 0;                   % total subplots index
%for index_row = 1:no_row,         % loop for row index
%    for index_col = 1:no_col          % loop for column index
%        index_loop = index_loop+1;
%
%        startx = inix + (index_col - 1) * delx;
%        starty = iniy + (no_row - index_row) * dely;
%
%%.......It's a bug in MATLAB
%        if alinewidth < 1.0
%           set(hsubplot(index_loop),'position',...
%              [ startx - 0.5 * alinewidth * (index_col-1), ...
%                starty - 0.5 * alinewidth * (index_row-1), delx ,dely]);
%%              [startx-1.0*alinewidth*(index_col-1), starty+1.5*alinewidth*(index_row-1), delx ,dely]);
%        else
%           set(hsubplot(index_loop),'position',[startx, starty, delx ,dely]);
%        end
%
%        subplot(hsubplot(index_loop));
 index_loop = no_subplot+1;              % total subplots index (reverse order)
for index_row = no_row:-1:1,             % loop for row index
    for index_col = no_col:-1:1          % loop for column index
        index_loop = index_loop - 1;
         startx = inix + (index_col - 1) * delx;
        starty = iniy + (no_row - index_row) * dely;
        POSITION = [startx, starty, delx ,dely];
 %.......It's kind of bug of MATLAB
        if alinewidth < 1.0
           POSITION =  [ startx - 0.5 * alinewidth * (index_col-1), ...
                         starty + 0.9 * alinewidth * (index_row-1), delx ,dely];
%          POSITION =  [startx-1.0*alinewidth*(index_col-1), starty+1.5*alinewidth*(index_row-1), delx ,dely]);
        end
         set(hsubplot(index_loop),'position',POSITION);
		
        subplot(hsubplot(index_loop));
         iscale = size(get(gca,'yscale'),2);  % 3:log, 6:linear
         % remove xlabels & xticklabels of subplots located in upper rows other than lowest row
         if index_row ~= no_row,
        	if ~(no_delta & index_row == (no_row - 1) & index_col == no_col),
            	set(get(gca,'xlabel'),'String',[])
            	set(gca,'xticklabel',[]);  %remove xticklabel
            end
        end
         % remove ylabels & yticklabels of subplots located in right columns other than leftmost column
        if index_col ~= 1,
           set(get(gca,'ylabel'),'String',[])
           set(gca,'yticklabel',[]);  %remove yticklabel
        end
         % remove first yticklabel of subplots located in lower rows other than highest row, linear yscale only
		% .... only linear scale
        if index_row ~= 1 & iscale == 6
           	a = get(gca,'ytick'); b = get(gca,'ylim');
			if a(length(a)) == b(length(b)), 
           		a = a(1:length(a)-1); 
           		set(gca,'ytick',a); 
           	end
        end
         % remove first xticklabel of subplots located in left columns other than rightmost column
		% .... only linear scale
		
		if ~no_delta,
  	      	if index_col ~= no_col & iscale == 6
       	       	a = get(gca,'xtick'); b = get(gca,'xlim');
           	   	if a(length(a)) == b(length(b)), 
          			a = a(1:length(a)-1); 
           			set(gca,'xtick',a); 
           		end
        	end
        else
	   	    if index_col == no_col & index_row == no_row - 1 & iscale == 6,
       	       	a = get(gca,'xtick'); 
          		a = a(2:length(a)); 
           		set(gca,'xtick',a); 
           	end	
        end	
    end
end
 % get back to initial unit
set(hfig,'unit','default')
for i=1:no_subplot,	set(hsubplot(i),'unit','default'),end
 % delete dummy subplots
if no_delta, for i = 1:no_delta, delete(hsubplot(no_subplot1+i)); 
    end
end
end
%-------------------------------------------------------------------------------------------------------------------
function [a] = read_obs_file(filnam, nobs,iplot)
%--------------------------------------------------------------------
% Usage: [A] = read_obs_file(filnam,nobs,iplot);
% This function reads an MT3D/RT3D ".OBS" file (eg., RT3D001.OBS) and 
% creates an array A(NOBS,NPTS) for further processing inside MATLAB.
% NOBS is the  number of observational points.
% It also optionally plots the data if IPLOT=1. The first argument 
% "filnam" is the name of the OBS file. Note that the two vectors
% A(1,i) and (A(2,i) contain a counter and the time (units depend
% on the simulation) respectively. The time series start from A(3,i).
% To make sure that the data is read correctly into MATLAB:
% 1. Type out the vector A(1,:). It should contain consecutive integers.
% 2. max(A(2,:)) should be equal to the total simulation time.
% 3. Use iplot=1; The plot can be used to see if the data is read properly.
% (This is the only intended use of this feature. More useful/customized plots 
% can be generated from inside MATLAB once the data is read successfully.)
% TROUBLESHOOTING: 
% If the data is not read in properly, the usual cause is that there are
% multiple datasets in the file. This happens, for example, when the
% MT3D/RT3D execution is terminated for some reason, leaving an ".OBS"
% file behind. When the execution is restarted, data is appended to the
% old file, if it exists. This utility makes use of the FSCANF 
% function which reads  ALL the data in the file. This will cause problems 
% if multiple datasets are interpreted as a single dataset.
% Please send corrections, if any,  to: 
% Phani, Geological Sciences Dept., Michigan State University, East Lansing.
% e-mail: phani@msu.edu
%---------------------------------------------------------------------
%tic;         % Get the clock running to compute time taken.
kmax=nobs+2; % Increment by 2 because we also have time and a counter;
nobs3=3*nobs;
fid1 = fopen(filnam,'r');
%disp('Reading data...');
%pause(1);
[dummy]=fscanf(fid1,'%c', 71); % Read the title into a dummy variable. We don't need it;
[loc] = fscanf(fid1,'%3f',nobs3); %Read the locations (I,J,K) of the points. This is not returned to the calling program, but it can be, if there is interest.
c = fscanf(fid1,'%12f'); %Now read the data (Note, '%12f' is needed to take care of an occasional negative sign in the otherwise '%11f' data. This is very important! Zero is sometimes computed as a small negative quantity.)
fclose(fid1);
%disp('Completed Reading...');
%pause(1);
nxmax=length(c);
nx=fix(nxmax/kmax); %use 'fix' in place of 'floor'; Also, 'fix' here will take care of datafiles that have their "tails"  uneven.
nxmax=nx*kmax;
a = reshape(c,kmax,nx); 

if (iplot==1) 
%   fprintf(1,'%s', 'Plotting Data...');
%   pause(1);
 x=a(2,:);
for k=3:kmax;
   y=a(k,:);
   hold on;
   h = plot(x,y);
   set(h, 'LineWidth',1.25);
end
hold off;
h=xlabel('Time');
set(h,'FontName','times');
set(h, 'FontSize',14);
h=ylabel('Concentration');
set(h,'FontName','times');
set(h, 'FontSize',14);
%axis square;
h=gca; %Get the handle to the current axis.
set(h,'TickLength',[0.02, 0.02]); %Change attributes using the handle.
set(h,'FontName','times');
set(h,'FontSize',14);
set(h, 'LineWidth', 1.5);
box on;
h=figure(1);
set(h, 'color','w');
end
end


%disp('All Done.'); % We are all done, except for the time estimate below.
%tim = toc;  % "toc" gives us the total time reckoned from the moment we called "tic" earlier. Nice feature!
%tim_min=tim/60.0; % seconds to  minutes
%fprintf(1, 'Total Elapsed Time in Minutes = %6.2f\n', tim_min);
%-------------------------------------End of READ_OBS_FILE-----------------------------