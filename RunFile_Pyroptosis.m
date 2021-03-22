%% File to run the pyroptosis pathway ODE model
% This is the main file that solves the system of ODEs (formulated in the
% pyroptosis_ODE.m) and creates data as well as single-dose plots.
%% Clear any previous data
clear
close all
%% Section 1: Set initial concentrations 
% Set up vector to store concentrations of each component
y0=zeros(15,1);
% where:
% % y(1) = nuclear NF-kB
% % y(2) = inactive NLRP3 (NLRP3i)
% % y(3) = active NLRP3 (NLRP3a)
% % y(4) = oligomerised/bound NLRP3 (NLRP3o)
% % y(5) = bound ASC (ASCb)
% % y(6) = cleaved caspase-1 (C1)
% % y(7) = cleaved gasdermin N terminal (GSDMD-N)
% % y(8) = Pro-interleukin-1b (Pro-IL-1b) (b=beta)
% % y(9) = cytoplasmic interleukin-1b (IL-1bc)
% % y(10) = external interleukin-1b (IL-1be)
% % y(11) = cytoplasmic interleukin-18 (IL-18c)
% % y(12) = external interleukin-18 (IL-18e)
% % y(13) = drug
% % y(14) = drug-NLRP3a complex
% % y(15) = relative cell volume (V)

%Specify drug dose
drug_dose =0;

% Set initial conditions for components intially > 0
y0(13)=drug_dose; % Initial drug concentration
y0(15) = 1;       % Initial relative cell volume

%% Section 2: Set up ODE solver options

% Set up options to use prior to inflammasome base formation:
options1 = odeset('RelTol',1e-4,'AbsTol',1e-4);

% Define the time simulations should be run over
tspan1 = [0 300];

%% Section 3: Set up NF-kB function
nfkb_0=0.25;        % Initial concentration of NF-kBn
h=0.55;             % maximum heigh elevation of the NF-kBn peak
s=0.8;              % skewness of the NF-kB peak
tau=10;             % time when the Nf-kBn peak occurs

% Set vector to store the NF-kBn relevant parameters
nfkb_vars=[nfkb_0, h, tau, s];

%% Section 4: Run ODE solver

% Run ODE solver 15s from file conserved_pyroptosis_ODE.m
[t,y] = ode15s(@(t,y) conserved_pyroptosis_ODEs(t,y,nfkb_vars), tspan1, y0, options1);

%% Section 5: Plot results (optional)


%Cut the data once volume>1.5 as cell has burst (for plots)
for k=1:length(y)
    if y(k,end)<=1.5
        yd(k,:)=y(k,:);
        td(k)=t(k);
    else
        yd(k,:)=NaN;
        
        td(k)=NaN;
    end
end


figure('DefaultLegendFontSize',20,'DefaultLegendFontSizeMode','manual', 'DefaultAxesFontSize', 20,'DefaultLineLineWidth', 4,'Units','normalized','Position',[0 0 1 1])
t = tiledlayout(2,4,'TileSpacing','none','Padding','none');
plotendtime=300;

% Plot NF-kB dynamics
nexttile
set(gca, 'ColorOrder',[0 0.7 0],'NextPlot', 'replacechildren');
nfkbn=@(t_in) h.* exp(-log(t_in./tau).^2/s) + nfkb_0;
plot(td,nfkbn(td))
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[NF-\kappaB_{n}]','Location','northeast')
xlim([0 plotendtime])
ylim([0 2])
yticks([0 1])
box on


% Plot NLRP3 dynamics
nexttile
if(drug_dose>0)
    set(gca, 'ColorOrder',[0 0 0.8; 0 0 0.8; 0 0 0.8; 0 0 0; 0.07 0.07 0.07], 'NextPlot', 'replacechildren');
    hold on
    plot(td,yd(:,2),':',td,yd(:,3),'-.',td,yd(:,4),td,yd(:,14))
    plot(1*ones(1,plotendtime),'--','LineWidth',2);
    hold off
    xlabel('Time (minutes)')
    ylabel('Concentration (a.u)')
    legend('[NLRP3_{i}]','[NLRP3_{a}]','[NLRP3_{o}]','[NLRP3_{a}\cdotdrug]','Location','northeast')
    xlim([0 plotendtime])
    ylim([0 2])
    yticks([0 1])
    box on
    yyaxis right
    ax=gca;
    ax.YAxis(2).Color = 'k';
    yticks([0.5])
    yticklabels({'n'})
else
    set(gca, 'ColorOrder',[0 0 0.8; 0 0 0.8; 0 0 0.8; 0.07 0.07 0.07], 'NextPlot', 'replacechildren');
    hold on
    plot(td,yd(:,2),':',td,yd(:,3),'-.',td,yd(:,4))
    plot(1*ones(1,plotendtime),'--','LineWidth',2);
    hold off
    xlabel('Time (minutes)')
    ylabel('Concentration (a.u)')
    legend('[NLRP3_{i}]','[NLRP3_{a}]','[NLRP3_{o}]','Location','northeast')
    xlim([0 plotendtime])
    ylim([0 2])
    yticks([0 1])
    box on
    yyaxis right
    ax=gca;
    ax.YAxis(2).Color = 'k';
    yticks([0.5])
    yticklabels({'\it{n}'})
end

% Plot ASCb  dynamics
nexttile
set(gca, 'ColorOrder',[1 0.9 0],'NextPlot', 'replacechildren');
plot(td,yd(:,5));
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[ASC_{b}]');
xlim([0 plotendtime])
ylim([0 2])
yticks([0 1])
box on

% Plot C1 dynamics
nexttile
set(gca, 'ColorOrder',[0.9 0 0 ],'NextPlot', 'replacechildren');
plot(td,yd(:,6));
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[C1]');
xlim([0 plotendtime])
ylim([0 2])
yticks([0 1])
box on

% Plot GSDMD-N dynamics
nexttile
set(gca, 'ColorOrder',[0.3 0.1 0.6 ],'NextPlot', 'replacechildren');
plot(td,yd(:,7));
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[GSDMD-N]');
xlim([0 plotendtime])
ylim([0 2])
yticks([0 1])
box on


% Plot IL-1b dynamics
nexttile
set(gca, 'ColorOrder',[0.3 0.4 0.1; 0.3 0.4 0.1; 0.3 0.4 0.1],'NextPlot', 'replacechildren');
plot(td,yd(:,8),':',td,yd(:,9),'-.',td,yd(:,10));
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[pro-IL-1\beta]','[IL-1\beta_{c}]','[IL-1\beta_{e}]','Location','northeast')
xlim([0 plotendtime])
ylim([0 2])
yticks([0 1])
box on


% Plot IL-18 dynamics
nexttile
set(gca, 'ColorOrder',[ 1 0.5 0.2; 1 0.5 0.2; 1 0.5 0.2],'NextPlot', 'replacechildren');
plot(td,yd(:,11),'-.',td,yd(:,12));
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[IL-18_{c}]','[IL-18_{e}]','Location','northeast')
xlim([0 plotendtime])
ylim([0 2])
yticks([0 1])
box on

% Plot relative cell Volume
nexttile
set(gca, 'ColorOrder',[0 0 0; 0.07 0.07 0.07],'NextPlot', 'replacechildren');
hold on
plot(td,yd(:,end));
plot(1.5*ones(1,plotendtime),'--','LineWidth',2);
hold off
xlabel('Time (minutes)')
ylabel('Volume')
legend('V','Location','northeast')
xlim([0 plotendtime])
ylim([0 2])
%yticks([0 1 1])
box on
yyaxis right
yticks([0.75]);
ax=gca;
ax.YAxis(2).Color = 'k';
yticklabels({'\it{V_c}'})

% End of File %
