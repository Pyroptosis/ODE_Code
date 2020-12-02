%% File to run the discretised version of the reduced pyroptosis model
%% Clear any previous data
clear
close all

%% Set parameter values
S1=1; % assume signal 1 is "on" here 
S2=1; % assume signal 2 is "on" here 
% Parameter values:
alpha1=0.025; % NLRP3 transcription rate
alpha2 = 0.08; % cleavage coefficient of GSDMD by caspase
alpha3=0.007; % Pro-ILib transcription rate
alpha4=0.8; % cleavage coefficient of IL-1b by caspase
alpha5=0.8; % cleavage coefficient of IL-18 by caspase
C1_50=0.3; % caspase cleavage half max occupancy
delta1=0.002; % decay rate of NLRP3i and NLRP3o 
delta2=0.004; % decay rate of Pro-Il-1b and Il-1b 
gammaC1=2; % Hill coefficient fro caspase cleavage
gammaNF=2; % Hill function coefficient - transcription
k1a=0.3; % NFkB cytoplasm->nucleus rate 
k1b=k1a/10; % NFkB cytoplasm<-nucleus rate 
k2=0.07; % inactive-->active NLRP3 rate 
k3=0.07; % NLRP3 oligomerisation rate
k4=0.02; % ACS + NLRP3o binding rate
k5=0.04; % Pro-caspase1 -> caspase1 cleavage rate
k6=0.8; % rate at which IL-1b can leave cell
k7=0.8; % rate at which IL-18 can leave cell
k8=0.1; % rate at which cell volume increases once pores formed
k_pTR=0; % rate at which TR binds to NLRP3a
k_mTR=0; % rate at which [TR.NLRP3a] complex dissociates
NF_50=0.3; % NFKB_50 - transcription half max occupancy of TFs on DNA
Vc=1.5;

%%
dt=1; % time step 
T =300/dt; % time scale
%%
y=zeros(T+1,2);

% Set initial conditions
y(1,2)=0.25; %Nf-kB nucleus
y(1,3)=0;% inactive NLRP3
y(1,4)=0; % active NLRP3
y(1,5)=0; % oligomerized NLRP3
y(1,7)=0; %bound ASC
y(1,9)=0; %caspase
y(1,11)=0; % GSDMD
y(1,12)=0; %Pro-IL-1b
y(1,13)=0;% IL-1bc
y(1,14)=0;% IL-1be
y(1,16)=0;% IL-18c
y(1,17)=0;% IL-18e
y(1,18)=1;%volume

t=1;
for t=1:T
    if y(t,5)<=1 % If inflammasome base has not yet been formed
        F=1;
    else
        F=0;
    end
    
    if y(t,18)<Vc % If cell volume has not reached critical level
    %NF-kB nucleus
    y(t+1,2)=(y(t,2)+k1a*F*dt)/(1+(k1b*dt)+(k1a*F*dt));
    % Hill function for transcription
    Hill1=(y(t+1,2)^(gammaNF))/((NF_50^(gammaNF))+y(t+1,2)^(gammaNF));
    % NLRP3i
    y(t+1,3)=(y(t,3)+alpha1*dt*(Hill1))/(1+(k2*dt)+(delta1*dt));
    % NLRP3a
    y(t+1,4)=(y(t,4)+k2*dt*(y(t+1,3)))/(1+(k3*dt)+(delta1*dt));%without drug
    %NLRP3o
    y(t+1,5)=y(t,5)+dt*k3*F*y(t+1,4);
    % ASC bound
    y(t+1,7)=(y(t,7)+dt*k4*(1-F)*y(t+1,5))/(1+dt*k4*(1-F)*y(t+1,5));
    % Caspase1
    y(t+1,9)=(y(t,9)+dt*k5*y(t+1,7))/(1+dt*k5*y(t+1,7));
    % Hill2
    Hill2=(y(t+1,9)^(gammaC1))/((C1_50^(gammaC1))+y(t+1,9)^(gammaC1));
    % GSDMD
    y(t+1,11)=(y(t,11)+dt*alpha2*Hill2)/(1+dt*alpha2*Hill2);
    % ProIL1b
    y(t+1,12)=(y(t,12)+dt*alpha3*Hill1)/(1+dt*alpha4*Hill2+delta2*dt); % Pro-IL1b
    % IL1b cytoplasm
    y(t+1,13)=(y(t,13)+alpha4*dt*Hill2*y(t+1,12))/(1+delta2*dt+k6*dt*y(t+1,11)); %IL-1b
    % IL1b external
    y(t+1,14)=y(t,14)+dt*k6*y(t+1,11)*y(t+1,13); % IL-1b c
    
    % IL-18 cytoplasm
    y(t+1,16)= (y(t,16)+alpha5*dt*Hill2*(1-y(t,17)))/((1+alpha5*dt*Hill2)*(1+k7*dt*y(t+1,11)));
 
    % IL-18 external
    y(t+1,17)=y(t,17)+dt*k7*y(t+1,11)*y(t+1,16);
    % Volume
    y(t+1,18)=y(t,18)/(1-k8*dt*y(t+1,11));
    else
        y(t+1,:)=NaN;
    end
end
%% Plotting results (optional)
figure('DefaultLegendFontSize',22,'DefaultLegendFontSizeMode','manual', 'DefaultAxesFontSize', 20,'DefaultLineLineWidth', 4,'Units','normalized','Position',[0 0 1 1])

t=linspace(0,300,T+1)
% Plot NF-kB dynamics
subplot(2,4,1)
set(gca, 'ColorOrder',[ 0 0.7 0],'NextPlot', 'replacechildren');
plot(t,y(:,2))
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[NF-\kappaB_{n}]','Location','northeast')
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
set(gca, 'ColorOrder',[ 1 0.9 0; ],'NextPlot', 'replacechildren');
plot(t,y(:,7))
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[ASC_{b}]');
xlim([0 60*5])
ylim([0 2])
box on
 
% Plot caspase1 dynamics
subplot(2,4,4)
set(gca, 'ColorOrder',[  0.9 0 0],'NextPlot', 'replacechildren');
plot(t,y(:,9))
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[C1]');
xlim([0 60*5])
ylim([0 2])
box on
 
% Plot GSDMD dynamics
subplot(2,4,5)
set(gca, 'ColorOrder',[ 0.3 0.1 0.6;],'NextPlot', 'replacechildren');
plot(t,y(:,11))
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[GSDMD-N]','Location','northeast');
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
set(gca, 'ColorOrder',[ 1 0.5 0.2; 1 0.5 0.2],'NextPlot', 'replacechildren');
plot(t,y(:,16),'-.',t,y(:,17))
%xline(At,'r--','LineWidth',3);
xlabel('Time (minutes)')
ylabel('Concentration (a.u)')
legend('[IL-18_{c}]','[IL-18_{e}]','Location','northeast')
xlim([0 60*5])
ylim([0 2])
box on
 
% Plot Volume
subplot(2,4,8)
set(gca, 'ColorOrder',[0 0 0],'NextPlot', 'replacechildren');
hold on
plot(t,y(:,18));
plot(1.5*ones(1,60*5),'k--','LineWidth',2);
%xline(At,'r--','LineWidth',3);
hold off
xlabel('Time (minutes)')
ylabel('Volume')
legend('V','V_{c}','Location','northeast')
xlim([0 60*5])
ylim([0 2])
box on
% End of File %