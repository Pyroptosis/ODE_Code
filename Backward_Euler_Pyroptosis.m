%% File to run the Backward Euler discretised version of the conserved pyroptosis pathway ODE model

%% Clear any previous data
clear
close all

%% Section 1: Set up vector

dt=1.0;          % time-step of simulations (mins)
T =300;          % final time we want to consider

% Set up vector to store concentrations of each component over time
y=zeros(T,13);
% where:
% % y(:,1) = nuclear NF-kB (NF-kBn) (k=kappa)
% % y(:,2) = inactive NLRP3 (NLRP3i)
% % y(:,3) = active NLRP3 (NLRP3a)
% % y(:,4) = oligomerised/bound NLRP3 (NLRP3o)
% % y(:,5) = bound ASC (ASCb)
% % y(:,6) = cleaved caspase-1 (C1)
% % y(:,7) = cleaved gasdermin N terminal (GSDMD-N)
% % y(:,8) = Pro-interleukin-1b (Pro-IL-1b) (b=beta)
% % y(:,9) = cytoplasmic interleukin-1b (IL-1bc)
% % y(:,10) = external interleukin-1b (IL-1be)
% % y(:,11) = cytoplasmic interleukin-18 (IL-18c)
% % y(:,12) = external interleukin-18 (IL-18e)
% % y(:,13) = relative cell volume (V)


%% Section 2: Set parameter values

% Set signals to be on
S1=1;               % signal 1 (initiation of NF-kB translocation)
S2=1;               % signal 2 (activation of NLRP3)

a=-1;               % sigmoid constant for step-function approximation for NLRP3o-ASC binding
alpha1=0.07;        % transcription rate of NLRP3 mediated by NF-kB
alpha2 = 0.1;       % cleavage rate of GSDMD by C1
alpha3=0.06;        % transcription rate of Pro-IL-1b mediated by NF-kB
alpha4=1;           % cleavage rate of Pro-IL-1b by C1
alpha5=1;           % cleavage rate of Pro-IL-18 by C1
b=2;                % sigmoid constant for step-function approximation for NLRP3o-ASC binding
c=1000;             % sigmoid power constant for step-function approximation for NLRP3o-ASC binding
C1_50=0.3;          % half-max value of C1 required for cleavage actions
delta1=0.002;       % decay rate of NLRP3i and NLRP3a
delta2=0.004;       % decay rate of Pro-IL-1b and IL-1bc
h=0.55;             % maximum heigh elevation of the NF-kBn peak
gammaC1=2;          % Hill function coefficient of C1 for cleavage actions
gammaNF=2;          % Hill function coefficient of NF-kBn for transcription actions
k1=0.7;             % activation rate of NLRP3i to NLRP3a
k2=1;               % oligomerisation rate of NLRP3a to NLRP3o
k3=0.04;            % binding rate of ASC to NLRP3o
k4=0.03;            % binding rate of C1 to ASC
k5=1;               % transport rate of IL-1bc out of cell to IL-1be
k6=1;               % transport rate of IL-18c out of cell to IL-18e
k7=0.2;             % rate at which cell volume increases
NF_50=0.3;          % half-max value of NF-kBn required for transcription actions
s=0.8;              % skewness of the NF-kB peak
tau=10;             % time when the Nf-kBn peak occurs
Vc=1.5;             % critical relative volume at which the cell ruptures

%% Section 3: Set up initial concentrations and conditions
% Set initial conditions for each component
y(1,1)=0.25;        % Initial NF-kBn concentration
y(1,13)=1;          % Initial relative cell volume
y(1,2:12)=0;        % All other components zero initially


%% Section 4: Run simulation
ta=1; % keeps track of time for NF-kB

% For all time-steps
for t=1:T-1

    ta=ta+dt; % keeps track of time for NF-kB
    
    % If cell volume is less than critical value
    if y(t,end)<Vc
       
        y(t+1,1)=h.* exp(-(log((ta)./tau).*log((ta)./tau))/s)+y(1,1);
        
        % Define NF-kBn hill function for transcription of NLRP3i and Pro-IL-1b
        Hill1= ((y(t+1,1)-y(1,1))^(gammaNF))/((NF_50^(gammaNF))+((y(t+1,1)-y(1,1))^(gammaNF)));
        
        % Update NLRP3i
        y(t+1,2)= ((y(t,2)+dt*alpha1*Hill1)/(1+dt*S2*k1+dt*delta1));
        
        % Set up quadratic for backward euler of NLRP3a
        m1=dt*k2; % quadratic a term
        m2=1+dt*delta1; % quadratic b term
        m3=-(y(t,3)+dt*S2*k1*y(t+1,2)); % quadratic c term
        
        % Update NLRP3a (negative root omitted as y>0)
        y(t+1,3)=(-m2+sqrt((m2^2)-4*m1*m3))/(2*m1);
        
        % Update NLRP3o
        y(t+1,4)=y(t,4)+dt*k2*(y(t+1,3)*y(t+1,3));
        
        % Set up function for ASC to NLRP3o binding
        n1=((y(t+1,4)-a)/b)^(-c);
        F1=1/(1+n1);
        
        % Update bound ASC
        y(t+1,5)=(y(t,5)+dt*k3*F1*y(t+1,4))/(1+dt*k3*F1*y(t+1,4));
        
        % Update cleaved caspase 1
        y(t+1,6)=(y(t,6)+dt*k4*y(t+1,5))/(1+dt*k4*y(t+1,5));
        
        % Define caspase 1 cleavage hill function
        Hill2=(y(t+1,6)^gammaC1)/((C1_50^gammaC1)+(y(t+1,6)^gammaC1));
        
        % Update cleaved GSDMD
        y(t+1,7)=(y(t,7)+dt*alpha2*Hill2)/(1+dt*alpha2*Hill2);
        
        % Update Pro-IL-1b
        y(t+1,8)=(y(t,8)+dt*alpha3*Hill1)/(1+delta2*dt+dt*alpha4*Hill2);
        
        % Update cytoplasmic IL-1b
        y(t+1,9)= (y(t,9)+dt*alpha4*Hill2*y(t+1,8))/(1+dt*k5*y(t+1,7)+dt*delta2);
        
        % Update external IL-1b
        y(t+1,10)= y(t,10)+dt*k5*y(t+1,7)*y(t+1,9);
        
        % Update cytoplasmic IL-18
        y(t+1,11)= (y(t,11)+dt*alpha5*Hill2-dt*alpha5*Hill2*y(t,12))/(1+(dt*k6*y(t+1,7))+(dt*alpha5*Hill2+dt*dt*alpha5*k6*Hill2*y(t+1,7)));
        
        % Update external IL-18
        y(t+1,12)= y(t,12)+dt*k6*y(t+1,7)*y(t+1,11);
        
        % Update relative cell volume
        y(t+1,13)= (y(t,13))/(1-dt*k7*y(t+1,7));
        
    else
        % If cell volume is over the critical value, the cell bursts so all
        % concentrations are removed
        y(t+1,:)=NaN;
    end
end


%% Optional: Plot figure of concentrations of each component

figure('DefaultLegendFontSize',22,'DefaultLegendFontSizeMode','manual', 'DefaultAxesFontSize', 20,'DefaultLineLineWidth', 4,'Units','normalized','Position',[0 0 1 1])
t=linspace(1,T*dt,T);
% Plot NF-kB dynamics
subplot(2,4,1)
set(gca, 'ColorOrder',[ 0 0.7 0],'NextPlot', 'replacechildren');
hold on
plot(t,y(:,1))
hold off
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[NF-\kappaB_{n}]','Location','northeast')
xlim([0 150])
ylim([0 2])
box on

% Plot NLRP3 dynamics
subplot(2,4,2)
set(gca, 'ColorOrder',[0 0 0.8; 0 0 0.8; 0 0 0.8], 'NextPlot', 'replacechildren');
hold on
plot(t,y(:,2),':')
plot(t,y(:,3),'-.')
plot(t,y(:,4))
hold off
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[NLRP3_{i}]','[NLRP3_{a}]','[NLRP3_{o}]','Location','northeast')
xlim([0 150])
ylim([0 2])
box on

% Plot ASC dynamics
subplot(2,4,3)
set(gca, 'ColorOrder',[ 1 0.9 0; ],'NextPlot', 'replacechildren');
plot(t,y(:,5))
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[ASC_{b}]');
xlim([0 150])
ylim([0 2])
box on

% Plot caspase1 dynamics
subplot(2,4,4)
set(gca, 'ColorOrder',[  0.9 0 0],'NextPlot', 'replacechildren');
plot(t,y(:,6))
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[C1]');
xlim([0 150])
ylim([0 2])
box on

% Plot GSDMD dynamics
subplot(2,4,5)
set(gca, 'ColorOrder',[ 0.3 0.1 0.6;],'NextPlot', 'replacechildren');
plot(t,y(:,7))
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[GSDMD-N]','Location','northeast');
xlim([0 150])
ylim([0 2])
box on

% Plot Il-1b dynamics
subplot(2,4,6)
set(gca, 'ColorOrder',[0.3 0.4 0.1; 0.3 0.4 0.1; 0.3 0.4 0.1;],'NextPlot', 'replacechildren');
hold on
plot(t,y(:,8),':')
plot(t,y(:,9),'-.')
plot(t,y(:,10))
hold off
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[pro-IL-1\beta]','[IL-1\beta_{c}]','[IL-1\beta_{e}]','Location','northeast')
xlim([0 150])
ylim([0 2])
box on

% Plot Il-18 dynamics
subplot(2,4,7)
set(gca, 'ColorOrder',[ 1 0.5 0.2; 1 0.5 0.2],'NextPlot', 'replacechildren');
hold on
plot(t,y(:,11),'-.')
plot(t,y(:,12))
hold off
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[IL-18_{c}]','[IL-18_{e}]','Location','northeast')
xlim([0 150])
ylim([0 2])
box on

% Plot Volume
subplot(2,4,8)
set(gca, 'ColorOrder',[0 0 0],'NextPlot', 'replacechildren');
hold on
plot(t,y(:,13));
plot(1.5*ones(1,150),'k--','LineWidth',2);
hold off
xlabel('Time')
ylabel('Volume')
legend('V','V_{c}','Location','northeast')
xlim([0 150])
ylim([0 2])
box on
%End of File %


