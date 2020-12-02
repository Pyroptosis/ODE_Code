function yp = pyroptosis_ODE_cont(t,y)
% ODE function file
% Model variables
% y(1) = [NFkB_c] : conc. of NFkB in cytoplasm
% y(2) = [NFkB_n] : conc. of NFkB in nucleus
% y(3) = [NLRP3_i] : conc. of inactive NLRP3
% y(4) = [NLRP3_a] : conc. of active NLRP3
% y(5) = [NLRP3_o] : conc. of oligomerised NLRP3 (forming the inflammasome base)
% y(6) = [ASC_f] : conc. of free ASC
% y(7) = [ASC_b] : conc. of inflammasome-bound ASC 
% y(8) = [pro-C1] : conc. of pro-caspase-1
% y(9) = [C1] : conc. of inflammasome-bound (i.e. active) caspase-1
% y(10) = [GSDMD] : conc. of (uncleaved) GSDMD 
% y(11) = [GSDMD_N] : conc. of GSDMD with N-terminal domain cleaved
% y(12) = [pro-IL-1b] : conc. of pro-interleukin-1b
% y(13) = [IL-1b_c] : conc. of cleaved/mature interleukin-1b in the cytoplasm
% y(13) = [IL-1b_e] : conc. of cleaved/mature interleukin-1b in the extracellular domain
% y(14) = [pro-IL-18] : conc. of pro-interleukin-18
% y(15) = [IL-18_c] : conc. of cleaved/mature interleukin-18 in the cytoplasm
% y(16) = [IL-18_e] : conc. of cleaved/mature interleukin-18 in the extracellular domain
% y(18) = [TR] : unbound drug (tranilast)
% y(19) = [TR * NLRP3_a] : complex of drug and active NLRP3_o
% y(20) = V : cell volume 

S1=1; % signal 1 - here turned on 
S2=1; % signal 2 - here turned on

% Parameter values:
alpha1=0.06; % NLRP3 transcription rate%0.025
alpha2 = 0.05; % cleavage coefficient of GSDMD by caspase 0.08
alpha3=0.06; % Pro-ILib transcription rate 0.007
alpha4=0.5; % cleavage coefficient of IL-1b by caspase 0.8
alpha5=0.5; % cleavage coefficient of IL-18 by caspase
C1_50=0.3; % caspase cleavage half max occupancy
delta1=0.002; % decay rate of NLRP3i and NLRP3o 
delta2=0.004; % decay rate of Pro-Il-1b and Il-1b 
gammaC1=2; % Hill coefficient fro caspase cleavage
gammaNF=2; % Hill function coefficient - transcription
k2=0.1; % inactive-->active NLRP3 rate 0.07
k3=0.2; % NLRP3 oligomerisation rate 0.07
k4=0.04; % ACS + NLRP3o binding rate 0.022
k5=0.03; % Pro-caspase1 -> caspase1 cleavage rate 0.04
k6=0.5; % rate at which IL-1b can leave cell 0.8
k7=0.5; % rate at which IL-18 can leave cell
k8=0.08; % rate at which cell volume increases once pores formed 0.1
k_pTR=0.05; % rate at which TR binds to NLRP3a
k_mTR=0.0005; % rate at which [TR.NLRP3a] complex dissociates
NF_50=0.3; % NFKB_50 - transcription half max occupancy of TFs on DNA
nu=30; % scale parameter one of NF-kB
s=0.9; % scale parameter two of NF-kB
t=t/15; % scale parameter of time of NF-kB

% Sigmoid function for ASC binding parameters
a=-1;
c=100;
b=2;

%System of ODEs
yp=[ exp(-log(t)^2/s) * (s+2*log(t)) * 1/(nu*t^2*s) %(log(t)-2)*log(t)*exp(-(log(t)^2)/stdev*t)/(a*t^2);%1-(exp((-log(t)).^2)/0.2)./x)/1.5 %ODE for [NF-kB_c]
    - exp(-log(t)^2/s) * (s+2*log(t)) * 1/(nu*t^2*s) 
    alpha1* (y(2)-0.25)^(gammaNF) / ( NF_50^(gammaNF) + (y(2)-0.25)^(gammaNF) ) - S2*k2*y(3) - delta1*y(3) %ODE for [NLRP3_i]
    S2*k2*y(3) - k3*y(4) - k_pTR*y(18)*y(4) + k_mTR*y(18)-delta1*y(4) %ODE for [NLRP3_a]
    k3*y(4) %ODE for [NLRP3_o]
    -k4*y(6)* 1/(1+((y(5)-a)/b).^(-c)) %ODE for [ASC_f]
    +k4*y(6)* 1/(1+((y(5)-a)/b).^(-c)) %ODE for [ASC_b]
    -y(7)*k5*y(8) %ODE for [pro-C1]
    +y(7)*k5*y(8) %ODE for [C1]
    -alpha2*( (y(9))^gammaC1/( C1_50^gammaC1 + (y(9))^gammaC1))*y(10) %ODE for [GSDMD]
    +alpha2*( (y(9))^gammaC1/( C1_50^gammaC1 + (y(9))^gammaC1))*y(10) %ODE for [GSDMD-N]
    +alpha3* (y(2)-0.25)^(gammaNF) / ( NF_50^(gammaNF) + (y(2)-0.25)^(gammaNF) ) - alpha4*( y(9)^gammaC1/( C1_50^gammaC1 + y(9)^gammaC1))*y(12)-delta2*y(12) %ODE for [pro-IL-1b]
    + alpha4*( (y(9))^gammaC1/( C1_50^gammaC1 + (y(9))^gammaC1))*y(12)-y(11)*k6*y(13)-delta2*y(13) %ODE for [IL-1b_c]
    +y(11)*k6*y(13) %ODE for [IL-1b_e]
    - alpha5*( (y(9))^gammaC1/( C1_50^gammaC1 + (y(9))^gammaC1))*y(15) %ODE for [pro-IL-18]
    + alpha5*( (y(9))^gammaC1/( C1_50^gammaC1 + (y(9))^gammaC1))*y(15) - y(11)*k7*(y(16)) %ODE for [IL-18_c]
    + y(11)*k7*(y(16)) %ODE for [IL-18_e]
  -k_pTR*y(18)*y(4) + k_mTR*y(19) %18: tr drug
    +k_pTR*y(18)*y(4) - k_mTR*y(19) %20: [tr na] bound tranilast & active NLRP3 
    +k8*y(20)*y(11)]; %ODE for V
end
% End of File %