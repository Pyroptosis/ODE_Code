% File to run the pyroptosis model through ode file pyroptosis_ode

%% Clear any previous data
clear
close all

%% Choose filename to save data to (optional)
%filename='Filename.mat';
%% Set up to run ODE
M=eye(20);

% Define options for first run (until inflammasome formed) and second run
% (until cell rupture)
options1 = odeset('Mass',M,'MassSingular','yes', 'RelTol',1e-4,'AbsTol',1e-4,'Events', @switchrate);
options2 = odeset('Mass',M,'MassSingular','yes', 'RelTol',1e-4,'AbsTol',1e-4,'Events', @membrabeRupture);

%% Set initial conditions
IC_NFkBc=0.75; % Initial NF-kB cytoplasmic relative  concentration
IC_NFkBn=0.25; % Initial NF-kB nuclear relative  concentration
IC_ASCf=1; % Initial ASC (free) relative  concentration
IC_PC1=1; % Initial pro-caspase-1 relative  concentration
IC_GSDMD=1; % Initial GSDMD relative  concentration
IC_IL18=1; % Initial pro-IL-18 relative  concentration
IC_TR=0; % Initial TR, drug, relative  concentration
IC_V=1; % Inititial relative volume of cell

% Set initial condition vector
y0=[IC_NFkBc IC_NFkBn 0 0 0 IC_ASCf 0 IC_PC1 0 IC_GSDMD 0 0 0 0 IC_IL18 0 0 IC_TR 0 IC_V];

%% Run ODE until 1st switch event - inflammasome base formation
% Set timespan
tspan1 = [0 60*5];

% Set initial value of step function F
F=1;

% Run ODE
[t1,y1] = ode15s(@(t1,y1) pyroptosis_ode(y1,F), tspan1, y0, options1);

%% Run ODE until 2nd switch event - cell rupture

% Only occurs if 1st switch event occured before final simulation time
if t1(end)<(60*5)
    
% Set timespan to start at final time from 1st run of ODE   
tspan2 = [t1(end) 60*5];

% Set updated step function value
F=0;
% Note we can actually 
% Use the last data from y1 as initial condition
y0 = y1(end,:);

% Run ODE
[t2,y2] = ode15s(@(t2,y2) pyroptosis_ode(y2,F), tspan2, y0, options2);

% Merge the 1st and 2nd run solutions (if both occur)
t=[t1;t2];
y=[y1;y2];
else
    % Set 1st run solution as final solution (if 1st event doesn't occur in timeframe)
    t=t1;
    y=y1;
end


%% Plot the results as seperate subplots
% Set figure defaults (style for plotting)
figure('DefaultLegendFontSize',22,'DefaultLegendFontSizeMode','manual', 'DefaultAxesFontSize', 20,'DefaultLineLineWidth', 4,'Units','normalized','Position',[0 0 1 1])

At=60*3.5;
% Plot NF-kB dynamics
subplot(2,4,1)
set(gca, 'ColorOrder',[0 0.7 0; 0 0.7 0],'NextPlot', 'replacechildren');
plot(t,y(:,1),':',t,y(:,2))
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[NF-\kappaB_{c}]','[NF-\kappaB_{n}]','Location','northeast')
xlim([0 60*5])
ylim([0 2])
box on

% Plot NLRP3 dynamics
subplot(2,4,2)
set(gca, 'ColorOrder',[0 0 0.8; 0 0 0.8; 0 0 0.8], 'NextPlot', 'replacechildren');
hold on
plot(t,y(:,3),':',t,y(:,4),'-.',t,y(:,5))
plot(1*ones(1,60*5),'b--','LineWidth',2);
hold off
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[NLRP3_{i}]','[NLRP3_{a}]','[NLRP3_{o}]','n','Location','northeast')
xlim([0 60*5])
ylim([0 2])
box on

% Plot ASC dynamics
subplot(2,4,3)
set(gca, 'ColorOrder',[1 0.9 0; 1 0.9 0; ],'NextPlot', 'replacechildren');
plot(t,y(:,6),':',t,y(:,7))
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[ASC_{f}]','[ASC_{b}]');
xlim([0 60*5])
ylim([0 2])
box on
 
% Plot caspase1 dynamics
subplot(2,4,4)
set(gca, 'ColorOrder',[ 0.9 0 0; 0.9 0 0],'NextPlot', 'replacechildren');
plot(t,y(:,8),':',t,y(:,9))
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[pro-C1]','[C1]');
xlim([0 60*5])
ylim([0 2])
box on
 
% Plot GSDMD dynamics
subplot(2,4,5)
set(gca, 'ColorOrder',[0.3 0.1 0.6; 0.3 0.1 0.6;],'NextPlot', 'replacechildren');
plot(t,y(:,10),':',t,y(:,11))
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[GSDMD]','[GSDMD-N]','Location','northeast');
xlim([0 60*5])
ylim([0 2])
box on
 
% Plot Il-1b dynamics
subplot(2,4,6)
set(gca, 'ColorOrder',[0.3 0.4 0.1; 0.3 0.4 0.1; 0.3 0.4 0.1;],'NextPlot', 'replacechildren');
plot(t,y(:,12),':',t,y(:,13),'-.',t,y(:,14))
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[pro-IL-1\beta]','[IL-1\beta_{c}]','[IL-1\beta_{e}]','Location','northeast')
xlim([0 60*5])
ylim([0 2])
box on
 
% Plot Il-18 dynamics
subplot(2,4,7)
set(gca, 'ColorOrder',[1 0.5 0.2; 1 0.5 0.2; 1 0.5 0.2],'NextPlot', 'replacechildren');
plot(t,y(:,15),':',t,y(:,16),'-.',t,y(:,17))
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[pro-IL-18]','[IL-18_{c}]','[IL-18_{e}]','Location','northeast')
xlim([0 60*5])
ylim([0 2])
box on
 
% Plot Volume
subplot(2,4,8)
set(gca, 'ColorOrder',[0 0 0],'NextPlot', 'replacechildren');
hold on
plot(t,y(:,20));
plot(1.5*ones(1,60*5),'k--','LineWidth',2);
%xline(At,'r--','LineWidth',3);
hold off
xlabel('Time (minutes)')
ylabel('Volume')
legend('V','V_{c}','Location','northeast')
xlim([0 60*5])
ylim([0 2])
box on

% Save figure (optional)
% savefig('Figure.fig')
%% Save data to file (optional)
% save(filename)
%% Define switch event functions

% Event 1: inflammasome base formation once NLRP3o=y(5)=1
function [value, isterminal, direction] = switchrate(tspan, y)
            value      = (y(5) - 1);  
            isterminal = 1;   % Stop the integration
            direction  = 1;   % Can be reached from above
end
  
% Event 2: cell rupture once V=y(2)=terminal_vol
function [value, isterminal, direction] = membrabeRupture(tspan, y)
    terminal_vol = 1.5; % volume value when membrane completely ruptures. 
    value      = (y(20) - terminal_vol);  
    isterminal = 1;   % Stops the integration.  
    direction  = 0; % Value can be reached from above or below 
end
