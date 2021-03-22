function yp = conserved_pyroptosis_ODEs(t,y,nfkb_vars)
%% Function to solve the conserved system of ODEs for pyroptosis pathway
% % Run through the file RunFile_Pyroptosis
% Input:
% % t =  time space of simulations
% % y =  vector of variables of system of equations
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

% % nfkb_vars = parameters for the NF-kB peak
% Output: 
% % y(t) =  vector of variables of system of equations over each time-step

%% Set up NF-kB function
% Nuclear NF-kB peak parameters
nfkb0=nfkb_vars(1); h=nfkb_vars(2); tau=nfkb_vars(3); s=nfkb_vars(4);
% Function at each time position
nfkbn=@(t_in) h.* exp(-log(t_in./tau).^2/s)+nfkb0;

%% Parameter settings
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
k_pTR=0.005;        % rate at which TR binds to NLRP3a
k_mTR=0.00005;      % rate at which [TR.NLRP3a] complex dissociates

%% System of ODEs to be solved for each time
yp=[nfkbn(t) % NF-kB
    alpha1* (nfkbn(t)-nfkb0)^(gammaNF) / ( NF_50^(gammaNF) + (nfkbn(t)-nfkb0)^(gammaNF) ) - S2*k1*y(2) - delta1*y(2) % NLRP3i
    S2*k1*y(2) - (k2*y(3)*y(3)) - k_pTR*y(13)*y(3) + k_mTR*y(14)-delta1*y(3)% NLRP3a
    k2*y(3)*y(3) % NLRP3o
    +k3*(1-y(5))* 1/(1+((y(4)-a)/b).^(-c)) % ASCb
    +y(5)*k4*(1-y(6)) % C1
    +alpha2*( (y(6))^gammaC1/( C1_50^gammaC1 + (y(6))^gammaC1))*(1-y(7)) % GSDMD-N
    +alpha3* (nfkbn(t)-nfkb0)^(gammaNF) / ( NF_50^(gammaNF) + (nfkbn(t)-nfkb0)^(gammaNF) ) - alpha4*( y(6)^gammaC1/( C1_50^gammaC1 + y(6)^gammaC1))*y(8)-delta2*y(8) % Pro-IL-1b
    + alpha4*( (y(6))^gammaC1/( C1_50^gammaC1 + (y(6))^gammaC1))*y(8)-y(7)*k5*y(9)-delta2*y(9) % IL-1bc
    +y(7)*k5*y(9) % IL-1be
    + alpha5*( (y(6))^gammaC1/( C1_50^gammaC1 + (y(6))^gammaC1))*(1-y(11)-y(12)) - y(7)*k6*(y(11)) % IL-18c
    + y(7)*k6*(y(11)) % IL-18e
    -k_pTR*y(13)*y(3) + k_mTR*y(14) % Drug
    +k_pTR*y(13)*y(3) - k_mTR*y(14) % Drug bound to NLRP3a
    +k7*y(7)*y(15)]; % Volume V
end
% end of function %