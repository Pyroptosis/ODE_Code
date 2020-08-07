function yp = pyroptosis_ode(y,F)
% y(1) = NFkB_c = NFkB in cytoplasm
% y(2) = NFkB_n = NFkB_tot - NFkB = NFkB in nucleus
% y(3) = NLRP3_i = inactive NLRP3
% y(4) = NLRP3_a = active NLRP3
% y(5) = NLRP3_o = oligomerised NLRP3
% y(6) = ASC 
% y(7) = ASC bound to [NLRPR_o + ASC] complex
% y(8) = pC1 = pro-caspase-1
% y(9) = C1 = caspase-1
% y(10) = GSDMD (not cleaved)
% y(11) = GSDMD_N
% y(12) = pIL1b = pro-interleukin-1B
% y(13) = IL1b_c = interleukin-1B in cytoplasm
% y(14) = IL1b_ex in the extracellular domain
% y(15) = p1L18 = pro-interleukin-18
% y(16) = 1L18_c = interleukin-18 in cytoplasm
% y(17) = 1L18_ex = interleukin-18 in extracellular domain
% y(18) = drug 
% y(19) = drug bound to NLRP3_a
% y(20) = volume 

S1=1; % assume signal 1 is "on" here 
S2=1; % assume signal 2 is "on" here 
% Parameter values:
alpha1=0.02; % NLRP3 transcription rate
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
k_pTR=0.003; % rate at which TR binds to NLRP3a
k_mTR=k_pTR/100; % rate at which [TR.NLRP3a] complex dissociates
NF_50=0.1; % NFKB_50 - transcription half max occupancy of TFs on DNA



%Set ODEs
yp=[-S1*k1a*y(1)*F + k1b*y(2) %1: NFkB cytoplasm
    +S1*k1a*y(1)*F - k1b*y(2) %2: NFkB nucleus
    alpha1* y(2)^(gammaNF) / ( NF_50^(gammaNF) + y(2)^(gammaNF) ) - S2*k2*y(3)-delta1*y(3) %3: NLRP3 inactive.
    S2*k2*y(3) - k3*y(4)*F - k_pTR*y(18)*y(4) + k_mTR*y(20)-delta1*y(4) %4: NLRP3 active
    k3*y(4)*F %5: NLRP3 oligomerised
    -k4*y(6)*(1-F) %6: ASC unbound
    +k4*y(6)*(1-F) %7: ASC bound
    -y(7)*k5*y(8)%8: procasp1
    +y(7)*k5*y(8)%9: casp1
    -alpha2*( (y(9))^gammaC1/( C1_50^gammaC1 + (y(9))^gammaC1))*y(10) %10: GSDMS
    +alpha2*( (y(9))^gammaC1/( C1_50^gammaC1 + (y(9))^gammaC1))*y(10) %11: GSDMD-N
    +alpha3* y(2)^(gammaNF) / ( NF_50^(gammaNF) + y(2)^(gammaNF) ) - alpha4*( y(9)^gammaC1/( C1_50^gammaC1 + y(9)^gammaC1))*y(12)-delta2*y(12)  %12: pro-IL-1b
    + alpha4*( (y(9))^gammaC1/( C1_50^gammaC1 + (y(9))^gammaC1))*y(12)-y(11)*k6*y(13)-delta2*y(13) %13: IL-1b (cytoplasm)
    +y(11)*k6*y(13) %14: IL-1b (external).
    - alpha5*( (y(9))^gammaC1/( C1_50^gammaC1 + (y(9))^gammaC1))*y(15)  %15: pro-IL-18
    + alpha5*( (y(9))^gammaC1/( C1_50^gammaC1 + (y(9))^gammaC1))*y(15) - y(11)*k7*(y(16))  %16: IL-18 (cytoplasm)
    + y(11)*k7*(y(16))  %17: IL-18 (external)
    -k_pTR*y(18)*y(4) + k_mTR*y(19) %18: tr drug
    +k_pTR*y(18)*y(4) - k_mTR*y(19) %20: [tr na] bound tranilast & active NLRP3 
    +k8*y(20)*y(11) ;%21: cell volume
    
 ];


end