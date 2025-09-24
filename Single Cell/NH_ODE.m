function dydt = NH_ODE(t,y,p)

%% Select modules to use
Na_clamp = 0;

flag_ECC = 1;    % if 0, module clamped
flag_cam = 1;    % if 0, module clamped
flag_CaMKII = 1; % if 0, module clamped
flag_BAR = 1;    % if 0, module clamped

%% Collect params and ICs for each module

ny_ECC = 92; % NH - ECC 92
ny_cam = 15; % 15*3 NH CaM
ny_CaMKII = 6; % 6 NH CaMKII module 137-143
ny_BAR = 41; % 41 NH in BAR module

% Allocate ICs for each moduel
% Ca_j is y(36), Ca_sl is y(37), Ca_cytosol is y(38)
% y(54) -> y(59) are state transitions for myofilament model
% y(60) -> y(65) are state transitions for mode 1 junctional LCCs
% y(66) -> y(71) are state transitoins for mode 2 junctional LCCs
% y(72) -> y(77) are state transitions for mode 1 sarcolemmal LCCs
% y(78) -> y(83) are state transitions for mode 2 sarcolemmal LCCs

y_ecc = y(1:92); 
y_camDyad = y(93:107);
y_camSL = y(108:122);
y_camCyt = y(123:137);
y_CaMKII = y(138:143);
y_BAR = y(144:184);

CaMKII_inhit            = 1;
CaMKII_double           = 1;


cycleLength         = p(1);
AF_index            = p(2);
prot_index          = p(3);
Ligtot              = p(4);
Currents_record     = p(5);
Output_current      = p(6);
gender_flag           = p(7);
Female_GNa_ODE        = p(8);
Female_Gk1_ODE        = p(9);
Female_Gkur_ODE       = p(10);
Female_Gkach_ODE      = p(11);
Female_Ito_ODE        = p(12);
Female_CSQ_ODE        = p(13);
Female_NaK_ODE        = p(14);
Female_SK_ODE         = p(15);

gender_specific_array = [gender_flag Female_GNa_ODE Female_Gk1_ODE Female_Gkur_ODE Female_Gkach_ODE ... 
      Female_Ito_ODE Female_CSQ_ODE Female_NaK_ODE Female_SK_ODE];


% Used for population building
INa_Scale               = p(16);
INaL_Scale              = p(17);
INab_Scale              = p(18);
INaK_Scale              = p(19);
Itof_Scale              = p(20);
IKur_Scale              = p(21);
IK2p_Scale              = p(22);
IKr_Scale               = p(23);
IKs_Scale               = p(24);
IK1_Scale               = p(25);
IKp_Scale               = p(26);
IKach_Scale             = p(27);
ISK_Scale               = p(28);
ICaL_Scale              = p(29);
ICab_Scale              = p(30);
ICap_Scale              = p(31); % IS this VPMCA
INCX_Scale              = p(32);
IClCa_Scale             = p(33);
IClb_Scale              = p(34);
Jrel_Scale              = p(35); % VRyR?
Jserca_Scale            = p(36);
Cleft_Buffer_Scale      = p(37);
Cytosol_Buffer_Scale    = p(38);
L_type_k1o_Scale        = p(39);
L_type_Scale_camkii     = p(40);
PLB_kmf_Scale           = p(41);
RyR_EC50_Scale          = p(42);
Ca_SR_leak_Scale        = p(43);
CSQ_B_max_Scale         = p(44);
KiCa_Scale              = p(45);
KoCa_Scale              = p(46);
Kim_Scale               = p(47);
Kom_Scale               = p(48);



if prot_index == 100
    DTE_flag = p(49);
    Prot_interval = p(50);
else
    DTE_flag = 0;
    Prot_interval = 0;
end
%% Parameters

% Input parameter for stimulation protocols
% prot_input_par = 400;

% VClamp_V_Step =0;
LCCtotBA    = 0.025;           % [uM] - [umol/L cytosol]
RyRtotBA    = 0.135;           % [uM] - [umol/L cytosol]
PLBtotBA    = 38;          % [uM] - [umol/L cytosol]
TnItotBA    = 70;              % [uM] - [umol/L cytosol]
IKstotBA    = 0.025;           % [uM] - [umol/L cytosol]
ICFTRtotBA  = 0.025;         % [uM] - [umol/L cytosol]
PP1_PLBtot  = 0.89;          % [uM] - [umol/L cytosol]
PLMtotBA    = 48;              % [uM] - [umol/L cytosol] as in Yang & Saucerman (mouse) model
MyototBA    = 70;              % [uM] - [umol/L cytosol] as TnI
IKrtotBA    = 0.025;           % [uM] - [umol/L cytosol] as IKs
IKurtotBA   = 0.025;          % [uM] - [umol/L cytosol] as IKr
INatotBA    = 0.025;           % [uM] - [umol/L cytosol] as IKr
IClCatotBA  = 0.025;         % [uM] - [umol/L cytosol] as ICFTR
ItototBA    = 0.025;           % [uM] - [umol/L cytosol] as IKr
IK1totBA    = 0.025;           % [uM] - [umol/L cytosol] as IKr
%CaM to Dyad
CaMtotDyad = 418;           % [uM]       % CaM (not used)
BtotDyad = 1.54 / 8.293e-4;   % [uM]
CaMKIItotDyad = 120;        % [uM]
CaNtotDyad = 3.6;           % [uM]
PP1totDyad = 96.5;          % [uM] % PP1

%CaM to SL
CaMtotSL = 5.65;            % [uM]       % CaM (not used)
BtotSL = 24.2;              % [uM]
CaMKIItotSL = 0.099516;      % [uM]
CaNtotSL = 3e-3;            % [uM]
PP1totSL = 0.57;            % [uM] % PP1

%CaM to Cyt
CaMtotCyt = 5.65;           % [uM]       % CaM (not used)
BtotCyt = 24.2;             % [uM]
CaMKIItotCyt = 0.099516; % [uM]
CaNtotCyt = 3e-3;           % [uM]
PP1totCyt = 0.57;           % [uM] % PP1
% CKIIOE = 0; % not being used

% Parameters for CaMKII module
LCCtotDyad = 28.26;       % [uM] - Total Dyadic [LCC] - (umol/l dyad)
LCCtotSL = 0.0846;          % [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
RyRtot = 382.6;             % [uM] - Total RyR (in Dyad)
PP1_dyad = 95.7;            % [uM] - Total dyadic [PP1] % why different from PP1totDyad?
PP1_SL = 0.57;              % [uM] - Total Subsarcolemmal [PP1]
PP2A_dyad = 95.76;          % [uM] - Total dyadic PP2A
OA = 0;                     % [uM] - PP1/PP2A inhibitor Okadaic Acid
PLBtot = 38;                % [uM] - Total [PLB] in cytosolic units
NaVtot = 30; % also used for Itof and IKur



K = 135; % [mM]
Mg = 1;  % [mM]
%% Distribute parameters by module

%% CaM module
CaDyad = y(36)*1e3; % from ECC model, *** Converting from [mM] to [uM] ***
compart_dyad = 2;
Btot = 0;
% ** NOTE: Btotdyad being sent to the dyad camODEfile is set to zero, but is used below for transfer between SL and dyad
pCaMDyad = [K,Mg,CaMtotDyad,Btot,CaMKIItotDyad,CaNtotDyad,PP1totDyad,CaDyad,cycleLength,compart_dyad]; %added by NH
% pCaMDyad = [K, Mg, CaMtotDyad, 0, CaMKIItotDyad, CaNtotDyad, PP1totDyad, CaDyad, cycleLength, compart_dyad]; %old one
CaSL = y(37)*1e3; % from ECC model, *** Converting from [mM] to [uM] ***
compartSL = 1;
pCaMSL = [K, Mg, CaMtotSL, BtotSL, CaMKIItotSL, CaNtotSL, PP1totSL, CaSL, cycleLength, compartSL];
CaCyt = y(38)*1e3; % from ECC model, *** Converting from [mM] to [uM] ***
compartCyt = 0;
pCaMCyt = [K, Mg, CaMtotCyt, BtotCyt, CaMKIItotCyt, CaNtotCyt, PP1totCyt, CaCyt, cycleLength, compartCyt];

%% CaMKII phosphorylation module 
CaMKIIact_Dyad = CaMKIItotDyad.*(y(100)+y(101)+y(102)+y(103)); % Multiply total by fraction
CaMKIIact_SL = CaMKIItotSL.*(y(115)+y(116)+y(117)+y(118));
%PP1_PLB_avail = y(ny_ECC+3*ny_cam+ny_CaMKII+22)./PP1_PLBtot + .0091;  % Active PP1 near PLB / total PP1 conc + basal value
PP1_PLB_avail = 1 - y(167)/PP1_PLBtot + 0.081698; % NEW ODEs - Mouse model % PP1_PLB_avail should be 100% (1) with NO ISO (Derived). With ISO -14% in Mouse model


pCaMKII = [CaMKIIact_Dyad,LCCtotDyad,RyRtot,PP1_dyad,PP2A_dyad,OA,PLBtot,...
           NaVtot,CaMKIIact_SL,LCCtotSL,PP1_SL,PP1_PLB_avail];

LCC_CKdyadp = (y(139)./LCCtotDyad)*L_type_Scale_camkii; % fractional CaMKII-dependent LCC dyad phosphorylation
RyR_CKp = y(141)./RyRtot;         % fractional CaMKII-dependent RyR phosphorylation
PLB_CKp = y(142)./PLBtot;         % fractional CaMKII-dependent PLB phosphorylation
NaV_CKp = y(140)./NaVtot;         % fractional CaMKII-dependent NaV phosphorylation

% Need to change before populations
Itof_CKp = NaV_CKp;
IK1_CKp = NaV_CKp;
Gapjunct_CKp = NaV_CKp;
Ikur_CKp = NaV_CKp;

%% BAR (PKA phosphorylation) module

pBAR = [Ligtot,LCCtotBA,RyRtotBA,PLBtotBA,TnItotBA,IKstotBA,ICFTRtotBA,...
    PP1_PLBtot,PLMtotBA,MyototBA,IKrtotBA,IKurtotBA,INatotBA,IClCatotBA,...
    ItototBA,IK1totBA,AF_index];

LCCa_PKAp = y(171)./LCCtotBA;
LCCb_PKAp = y(172)./LCCtotBA;
PLB_PKAn = (PLBtotBA - y(169))./PLBtotBA; % non-phosphorylated PLB targets
RyR_PKAp = y(173)./RyRtotBA;
TnI_PKAp = y(174)./TnItotBA;
IKs_PKAp = y(175)./IKstotBA;
PLM_PKAp = y(170)./PLMtotBA;
Myo_PKAp = y(180)./MyototBA;
IKr_PKAp = y(176)./IKrtotBA;
IKur_PKAp = y(177)./IKurtotBA;
INa_PKAp = y(181)./INatotBA;
IClCa_PKAp = y(179)./IClCatotBA;
Ito_PKAp = y(182)./ItototBA;
IK1_PKAp = y(183)./IK1totBA;

% gender_specific_array = [gender_flag Female_GNa_ODE Female_Gk1_ODE Female_Gkur_ODE Female_Gkach_ODE ... 
%       Female_Ito_ODE Female_CSQ_ODE Female_NaK_ODE Female_SK_ODE];


pECC = [cycleLength,AF_index,prot_index, CaMKII_inhit, CaMKII_double, INa_Scale, INaL_Scale, INab_Scale... %1-8
    INaK_Scale Itof_Scale IKur_Scale IK2p_Scale IKr_Scale IKs_Scale IK1_Scale IKp_Scale IKach_Scale... %9-17
    ISK_Scale ICaL_Scale ICab_Scale ICap_Scale INCX_Scale IClCa_Scale IClb_Scale Jrel_Scale Jserca_Scale... % 18-26
    Cleft_Buffer_Scale Cytosol_Buffer_Scale L_type_k1o_Scale PLB_kmf_Scale... % 27-30
    RyR_EC50_Scale Ca_SR_leak_Scale CSQ_B_max_Scale KiCa_Scale KoCa_Scale Kim_Scale Kom_Scale... % 31-37
    Na_clamp Currents_record Output_current...% 38-40
    gender_specific_array, DTE_flag, Prot_interval]; %-NH updated version

pPhos = [NaV_CKp INa_PKAp LCCb_PKAp LCC_CKdyadp RyR_CKp RyR_PKAp PLB_CKp PLB_PKAn PLM_PKAp TnI_PKAp ...  %1 - 10
    IKs_PKAp IKr_PKAp Ikur_CKp IKur_PKAp IClCa_PKAp Myo_PKAp Itof_CKp Ito_PKAp IK1_PKAp IK1_CKp Gapjunct_CKp LCCa_PKAp]; %11-22
%% Solve dydt in each module

if flag_ECC==0,
    dydt_ecc = zeros(1,length(y_ecc))';
else
    dydt_ecc = NH_ECC(t,y_ecc,pECC,pPhos);
end

if flag_cam==0,
    dydt_camDyad = zeros(1,length(y_camDyad))';
	dydt_camSL = zeros(1,length(y_camSL))';
	dydt_camCyt = zeros(1,length(y_camCyt))';
    JCaDyad = 0; JCaSL = 0; JCaCyt = 0;
else
    [dydt_camDyad,JCaDyad]  = NH_cam(t,y_camDyad,pCaMDyad);
    [dydt_camSL,JCaSL]  = NH_cam(t,y_camSL,pCaMSL);
    [dydt_camCyt,JCaCyt] = NH_cam(t,y_camCyt,pCaMCyt);
end

if flag_CaMKII==0,
    dydt_CaMKIIDyad = zeros(1,length(y_CaMKII))';
else
    dydt_CaMKIIDyad = NH_camkii(t,y_CaMKII,pCaMKII);
end

if flag_BAR==0,
    dydt_BAR = zeros(1,length(y_BAR))';
else
    dydt_BAR = NH_BetaAR(t,y_BAR,pBAR);
end


% incorporate Ca buffering from CaM, convert JCaCyt from uM/msec to mM/msec
dydt_ecc(36) = dydt_ecc(36) + 1e-3*JCaDyad;
dydt_ecc(37) = dydt_ecc(37) + 1e-3*JCaSL;
dydt_ecc(38) = dydt_ecc(38) + 1e-3*JCaCyt; 

% Cell geometry % from ECC module
cellLength = 100;     % cell length [um]
cellRadius = 10.25;   % cell radius [um]
Vcell = pi*cellRadius^2*cellLength*1e-15;    % [L]
Vmyo = 0.65*Vcell; VSL = 0.02*Vcell; Vdyad = 1*0.0539*.01*Vcell; 
% incorporate CaM diffusion between compartments
% kDyadSL = 3.6363e-16;	% [L/msec]
kSLmyo = 8.587e-15;     % [L/msec]
k0Boff = 0.0014;        % [s^-1] 
k0Bon = k0Boff/0.2;     % [uM^-1 s^-1] kon = koff/Kd
k2Boff = k0Boff/100;    % [s^-1] 
k2Bon = k0Bon;          % [uM^-1 s^-1]
% k4Boff = k2Boff;        % [s^-1]
k4Bon = k0Bon;          % [uM^-1 s^-1]
CaMtotDyad = sum(y_camDyad(1:6))+CaMKIItotDyad*sum(y_camDyad(7:10))+sum(y_camDyad(13:15));
Bdyad = BtotDyad - CaMtotDyad; % [uM dyad]
J_cam_dyadSL = 1e-3*(k0Boff*y_camDyad(1) - k0Bon*Bdyad*y_camSL(1)); % [uM/msec dyad]
J_ca2cam_dyadSL = 1e-3*(k2Boff*y_camDyad(2) - k2Bon*Bdyad*y_camSL(2)); % [uM/msec dyad]
J_ca4cam_dyadSL = 1e-3*(k2Boff*y_camDyad(3) - k4Bon*Bdyad*y_camSL(3)); % [uM/msec dyad]
J_cam_SLmyo = kSLmyo*(y_camSL(1)-y_camCyt(1)); % [umol/msec]
J_ca2cam_SLmyo = kSLmyo*(y_camSL(2)-y_camCyt(2)); % [umol/msec]
J_ca4cam_SLmyo = kSLmyo*(y_camSL(3)-y_camCyt(3)); % [umol/msec]

% adding comments: diffusion of CaM between 3 compartments.
dydt_camDyad(1) = dydt_camDyad(1) - J_cam_dyadSL;
dydt_camDyad(2) = dydt_camDyad(2) - J_ca2cam_dyadSL;
dydt_camDyad(3) = dydt_camDyad(3) - J_ca4cam_dyadSL;
dydt_camSL(1) = dydt_camSL(1) + J_cam_dyadSL*Vdyad/VSL - J_cam_SLmyo/VSL;
dydt_camSL(2) = dydt_camSL(2) + J_ca2cam_dyadSL*Vdyad/VSL - J_ca2cam_SLmyo/VSL;
dydt_camSL(3) = dydt_camSL(3) + J_ca4cam_dyadSL*Vdyad/VSL - J_ca4cam_SLmyo/VSL;
dydt_camCyt(1) = dydt_camCyt(1) + J_cam_SLmyo/Vmyo;
dydt_camCyt(2) = dydt_camCyt(2) + J_ca2cam_SLmyo/Vmyo;
dydt_camCyt(3) = dydt_camCyt(3) + J_ca4cam_SLmyo/Vmyo;
%% Collect all dydt terms

dydt = [dydt_ecc; dydt_camDyad; dydt_camSL; dydt_camCyt; dydt_CaMKIIDyad; dydt_BAR];