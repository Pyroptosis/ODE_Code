%% Run Single Cell Pyroptosis Model
% This is the main file that solves the system of ODEs (formulated in the
% pyroptosis_ODE.m) and creates data as well as single-dose plots.

%% Clear previous data and close figures
clear
close all

%% Generate data
%Specify data filename
data_filename = 'Filename.mat';
%Specify drug dose
drug_dose =0;

%% Set up ODE options
% Set up options to use prior to inflammasome base formation:
options1 = odeset('RelTol',1e-4,'AbsTol',1e-4);

%% Set initial conditions
ic_NFkB_c=0.75; % Initial conc. of cytoplasmic NF-kB
ic_NFkB_n=0.25; % Initial conc. of nuclear NF-kB
ic_ASC_f=1; % Initial conc. of free ASC
ic_pro_c1=1; % Initial conc. of pro-caspase-1
ic_GSDMD=1; % Initial conc. of (uncleaved) GSDMD
ic_IL_18=1;  % Initial conc. of pro-IL-18
ic_V=1; % Inititial cell volume
ic_TR=drug_dose; % Initial conc. of the drug, this is a function input

% Set the initial condition vector
y0=zeros(20,1);
y0(1) = ic_NFkB_c;
y0(2) = ic_NFkB_n;
y0(6) = ic_ASC_f;
y0(8) = ic_pro_c1;
y0(10) = ic_GSDMD;
y0(15) = ic_IL_18;
y0(18) = ic_TR;
y0(20) = ic_V;


tspan1 = [0.01 300];
% Run ODE
[t,y] = ode15s(@(t,y) pyroptosis_ODE_cont(t,y), tspan1, y0, options1);

save(data_filename);

%% Make subplots for one drug dose

% Cut the data once volume>1.5 as cell has burst
for k=1:length(y)
    if y(k,20)<=1.5
        yd(k,:)=y(k,:);
        td(k)=t(k);
    else
        yd(k,:)=NaN;
        
        td(k)=NaN;
    end
end

figure('DefaultLegendFontSize',20,'DefaultLegendFontSizeMode','manual', 'DefaultAxesFontSize', 20,'DefaultLineLineWidth', 4,'Units','normalized','Position',[0 0 1 1])
plotendtime=200;

% Plot NF-kB dynamics
subplot(2,4,1)
set(gca, 'ColorOrder',[0 0.7 0],'NextPlot', 'replacechildren');
plot(td,yd(:,2))
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[NF-\kappaB_{n}]','Location','northeast')
xlim([0 plotendtime])
ylim([0 2])
box on

% Plot NLRP3 dynamics
subplot(2,4,2)
set(gca, 'ColorOrder',[0 0 0.8; 0 0 0.8; 0 0 0.8], 'NextPlot', 'replacechildren');
hold on
plot(td,yd(:,3),':',td,yd(:,4),'-.',td,yd(:,5))
plot(1*ones(1,plotendtime),'b--','LineWidth',2);
hold off
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[NLRP3_{i}]','[NLRP3_{a}]','[NLRP3_{o}]','n','Location','northeast')
xlim([0 plotendtime])
ylim([0 2])
box on

% Plot ASC dynamics
subplot(2,4,3)
set(gca, 'ColorOrder',[1 0.9 0],'NextPlot', 'replacechildren');
plot(td,yd(:,7));
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[ASC_{b}]');
xlim([0 plotendtime])
ylim([0 2])
box on

% Plot caspase1 dynamics
subplot(2,4,4)
set(gca, 'ColorOrder',[ 0.9 0 0],'NextPlot', 'replacechildren');
plot(td,yd(:,9))
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[C1]');
xlim([0 plotendtime])
ylim([0 2])
box on

% Plot GSDMD dynamics
subplot(2,4,5)
set(gca, 'ColorOrder',[ 0.3 0.1 0.6;],'NextPlot', 'replacechildren');
plot(td,yd(:,11))
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[GSDMD-N]','Location','northeast');
xlim([0 plotendtime])
ylim([0 2])
box on

% Plot Il-1b dynamics
subplot(2,4,6)
set(gca, 'ColorOrder',[0.3 0.4 0.1; 0.3 0.4 0.1; 0.3 0.4 0.1],'NextPlot', 'replacechildren');
plot(td,yd(:,12),':',td,yd(:,13),'-.',td,yd(:,14))
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[pro-IL-1\beta]','[IL-1\beta_{c}]','[IL-1\beta_{e}]','Location','northeast')
xlim([0 plotendtime])
ylim([0 2])
box on

% Plot Il-18 dynamics
subplot(2,4,7)
set(gca, 'ColorOrder',[ 1 0.5 0.2; 1 0.5 0.2],'NextPlot', 'replacechildren');

plot(td,yd(:,16),'-.',td,yd(:,17))
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[IL-18_{c}]','[IL-18_{e}]','Location','northeast')
xlim([0 plotendtime])
ylim([0 2])
box on

% Plot Volume
subplot(2,4,8)
set(gca, 'ColorOrder',[0 0 0],'NextPlot', 'replacechildren');
hold on
plot(td,yd(:,20));
plot(1.5*ones(1,plotendtime),'k--','LineWidth',2);
%xline(At,'r--','LineWidth',3);
hold off
xlabel('Time (minutes)')
ylabel('Volume')
legend('V','V_{c}','Location','northeast')
xlim([0 plotendtime])
ylim([0 2])
box on

savefig('Figure.fig')
% End of File %