function [t, y,output] = NH_single_cell(duration, frequencey, AF, gender_flag,output_currents, ISO_param,Plot_Flag,beat_analysis)

%% Parameters for external modules

if gender_flag == 0 && AF == 0
    load yfin_nSR_1Hz_0_ISO_male.mat;
elseif gender_flag == 1 && AF == 0
    load yfin_nSR_1Hz_0_ISO_female.mat;
elseif gender_flag == 0 && AF == 1
    load yfin_cAF_1Hz_0_ISO_male.mat;
elseif gender_flag == 1 && AF == 1
    load yfin_cAF_1Hz_0_ISO_female.mat;
end

y0n = yfinal;

AF_index = AF;
freq = frequencey;
cycleLength = 1e3/freq;
ISO = ISO_param;
Currents_record = output_currents;

%Prot_index
%protocol(1) = 'Paced';
%protocol(4) = 'DAD';

prot_index = 1;

%% Collect all parameters


% 1 = female; 0 = male (baseline)

Female_diff = sex_diff;

if gender_flag == 1

    Female_GNa        = Female_diff(1,1);
    Female_Gk1        = Female_diff(1,2);
    Female_Gkur       = Female_diff(1,3);
    Female_Gkach      = Female_diff(1,4);
    Female_Ito        = Female_diff(1,5);
    Female_CSQ        = Female_diff(1,6);
    Female_NaK        = Female_diff(1,7);
    Female_SK         = Female_diff(1,8);

else

    Female_GNa        = 0;
    Female_Gk1        = 0;
    Female_Gkur       = 0;
    Female_Gkach      = 0;
    Female_Ito        = 0;
    Female_CSQ        = 0;
    Female_NaK        = 0;
    Female_SK         = 0;

end

gender_param = [gender_flag Female_GNa Female_Gk1 Female_Gkur Female_Gkach...
    Female_Ito Female_CSQ Female_NaK Female_SK]; %


%% Establish and define globals

if Currents_record  == 1

    global tStep tArray ICa_store Ito_store INa_store IK1_store
    global Jserca_store IKs_store Jleak_store Incx_store
    global INaK_store Ikr_store INabk_store
    global Lmyo_store Fmyo_store Vmax_store
    global IKp_store Iclca_store Iclbk_store Ipca_store Icabk_store Cai_store
    global I_app_store I_NaLstore I_kurstore I_k2pstore I_kistore I_kachstore
    global I_skstore I_ClCFTRstore I_pcastore I_J_SRCarelstore

    tStep = 1; tArray = zeros(1,1e6); ICa_store = zeros(1,1e6); Ito_store = zeros(1,1e6);
    INa_store = zeros(1,1e6); IK1_store = zeros(1,1e6);
    Jserca_store = zeros(1,1e6); IKs_store = zeros(1,1e6);
    Jleak_store = zeros(1e6,2); Incx_store = zeros(1,1e6);
    INaK_store = zeros(1,1e6); Ikr_store = zeros(1,1e6); INabk_store = zeros(1,1e6);
    Fmyo_store = zeros(1,1e6); Lmyo_store = zeros(1,1e6); Vmax_store = zeros(1,1e6);
    IKp_store = zeros(1,1e6); Iclca_store = zeros(1,1e6); Iclbk_store = zeros(1,1e6);
    Ipca_store = zeros(1,1e6); Icabk_store = zeros(1,1e6); Cai_store = zeros(1,1e6); I_J_SRCarelstore = zeros(1,1e6);
    I_app_store = zeros(1,1e6); I_NaLstore = zeros(1,1e6); I_kurstore = zeros(1,1e6);
    I_k2pstore = zeros(1,1e6); I_kistore = zeros(1,1e6); I_kachstore = zeros(1,1e6);
    I_skstore = zeros(1,1e6); I_ClCFTRstore = zeros(1,1e6); I_pcastore = zeros(1,1e6);

end


%% Run single simulation
prot_rate = freq;
period = 1000/prot_rate;
num_beats = floor(duration/period);   % add +1 here for extra beat
duration = (num_beats*period);

tic
par_SA = ones(1,33);
p = [cycleLength, AF_index, prot_index, ISO, Currents_record, output_currents,...
    gender_param, par_SA];
tspan = [0 duration]; % [ms]
options = odeset('RelTol',1e-6,'MaxStep',1);
[t,y] = ode15s(@NH_ODE,tspan,y0n,options,p);
yfinal = y(end,:);

toc

if Plot_Flag == 1
    % Nicer multi-panel figure
    figure('Color','w','Position',[100 100 600 800]); % white background, larger size

    % --- Subplot 1 ---
    subplot(3,1,1)
    plot(t, y(:,34), 'b', 'LineWidth', 1.5)
    ylabel('[Na]_i','FontSize',12,'FontWeight','bold')
    grid on

    % --- Subplot 2 ---
    subplot(3,1,2)
    plot(t, y(:,39), 'r', 'LineWidth', 1.5)
    ylabel('Vm','FontSize',12,'FontWeight','bold')
    grid on

    % --- Subplot 3 ---
    subplot(3,1,3)
    plot(t, y(:,38)*1000000, 'k', 'LineWidth', 1.5)
    ylabel('[Ca]_i (nM)','FontSize',12,'FontWeight','bold')
    xlabel('Time (ms)','FontSize',12,'FontWeight','bold')
    grid on

    % Adjust subplot spacing
    set(findall(gcf,'-property','FontName'),'FontName','Arial')
    sgtitle('Simulation Outputs','FontSize',16,'FontWeight','bold')

end



%% Rename outputs
if Currents_record == 1 || output_currents == 1

    tA = tArray(1:tStep); dVm = Vmax_store(1:tStep);
    ICa = ICa_store(1:tStep); Ito = Ito_store(1:tStep);
    INa = INa_store(1:tStep); IK1 = IK1_store(1:tStep);
    IKs = IKs_store(1:tStep); IKr = Ikr_store(1:tStep);
    IKp = IKp_store(1:tStep);
    IClCa = Iclca_store(1:tStep); IClbk = Iclbk_store(1:tStep);
    INCX = Incx_store(1:tStep); INaK = INaK_store(1:tStep);
    INabk = INabk_store(1:tStep); IPMCA = Ipca_store(1:tStep);
    ICabk = Icabk_store(1:tStep); Cai_tA = Cai_store(1:tStep);
    Jserca = Jserca_store(1:tStep); Jleak = Jleak_store(1:tStep,:);
    I_app = I_app_store(1:tStep); I_NaL = I_NaLstore(1:tStep); I_Kur = I_kurstore(1:tStep);
    I_K2P = I_k2pstore(1:tStep); I_Ki = I_kistore(1:tStep); I_Kach = I_kachstore(1:tStep);
    I_SK = I_skstore(1:tStep); I_ClFTR = I_ClCFTRstore(1:tStep); I_pca = I_pcastore(1:tStep);
    I_J_SRCarel = I_J_SRCarelstore(1:tStep);
end


%% Current Plot
if Currents_record == 1
    figure(2)
    set(gcf,'color','w')
    subplot(4,2,1); hold on, plot(tA,INa); title('INa');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])
    subplot(4,2,3); hold on, plot(tA,ICa); title('ICa');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])
    subplot(4,2,5); hold on, plot(tA,INCX); title('INCX');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])
    subplot(4,2,7); hold on, plot(tA,INaK); title('INaK');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])
    subplot(4,2,2); hold on, plot(tA,IKr); title('IKr');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])
    subplot(4,2,4); hold on, plot(tA,IKs); title('IKs');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])
    subplot(4,2,6); hold on, plot(tA,Ito); title('Ito');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])
    subplot(4,2,8); hold on, plot(tA,IK1); title('IK1');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)

    figure(3)
    hold on,
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,1);hold on, plot(tA,IKp); title('IKp');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,2);hold on, plot(tA,IClCa); title('IClCa');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,3);hold on, plot(tA,IClbk); title('IClbk');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,4);hold on, plot(tA,INabk); title('INabk');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,5);hold on, plot(tA,IPMCA); title('IPMCA');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,6);hold on, plot(tA,ICabk); title('ICabk');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,7);hold on, plot(tA,Cai_tA); title('Cai_tA');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,8);hold on, plot(tA,Jserca); title('Jserca');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,9);hold on, plot(tA,Jleak); title('Jleak');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)

    figure(4)
    hold on,
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,1);hold on, plot(tA,I_app); title('Iapp');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,2);hold on, plot(tA,I_NaL); title('INaL');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,3);hold on, plot(tA,I_Kur); title('IKur');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,4);hold on, plot(tA,I_K2P); title('IK2P');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,5);hold on, plot(tA,I_Ki); title('IKi');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,6);hold on, plot(tA,I_Kach); title('IKach');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,7);hold on, plot(tA,I_SK); title('ISK');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,8);hold on, plot(tA,I_ClFTR); title('IClFTR');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,9);hold on, plot(tA,I_pca); title('Ipca');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)

    figure(5)
    hold on,
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(3,2,1);hold on, plot(t,y(:,39)); title('Vm');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(3,2,3);hold on, plot(t,y(:,34)); title('[Na]_i');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])
    subplot(3,2,2);hold on, plot(tA,I_Ki); title('IK1');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])
    subplot(3,2,4);hold on, plot(tA,I_Kach); title('IKach');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])
    subplot(3,2,6);hold on, plot(tA,(I_Kach+I_Ki)); title('IKach + IK1');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])

    figure(6)
    hold on,
    plot(tA,ICa); title('ICa');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])

end

if beat_analysis == 1
    time = t; % (ms)
    Vm = y(:,39); % (mV)
    Ca = y(:,38); % (mM)
    Na = y(:,34); % (mM)
    dVm_array = (y(2:end,39)-y(1:end-1,39))./(t(2:end)-t(1:end-1));
    dVm = [dVm_array(1); dVm_array];
    AP_index = 2;
    period = 1000/freq;

    outputs = function_beat_analysis(time,Vm,Ca,Na,dVm,period,AP_index);
    dVm_max = outputs(1);
    Vm_max = outputs(2);
    RMP = outputs(3);
    AP_amp = outputs(4);
    APD90 = outputs(5);
    APD70 = outputs(6);
    APD50 = outputs(7);
    APD30 = outputs(8);
    Ca_max = outputs(9)*1000000;
    Ca_min = outputs(10)*1000000;
    CaT_amp = outputs(11)*1000000;
    CaT_rise = outputs(12);
    CaT_decay_50 = outputs(13);
    CaT_decay_63 = outputs(14);
    Na_min = outputs(15);
    VPLT = outputs(16);
    APD20 = outputs(17);
    %     disp(APD90)


    output = [dVm_max Vm_max -RMP AP_amp APD90 APD70 APD50 APD30 Ca_max...
        Ca_min CaT_amp CaT_rise CaT_decay_50 CaT_decay_63 Na_min VPLT freq APD20];
else

    output = [];
end



Vm = y(:,39);
Ca = y(:,38);

duration_extend = duration-period*5; duration_last = find (t == duration(end));
duration_extend_t = find(t > duration_extend);
duration_extend_time = duration_extend_t(1);
new_time = t(duration_extend_time:duration_last);
new_Vm = Vm(duration_extend_time:duration_last);
new_Ca = Ca(duration_extend_time:duration_last);


time = t; % (ms)
Vm = y(:,39); % (mV)
Ca = y(:,38); % (mM)
AP_index = 2;
figures = 1;
outputs = function_beat_analysis_2017_alternans(time,Vm,Ca,period,AP_index, figures); %only use for high frequency
APD90_1 = outputs(1);
APD90_2 = outputs(2);
Calcium_amp_1 = outputs(3);
Calcium_amp_2 = outputs(4);
APD_difference = APD90_1-APD90_2;
Calcium_difference = Calcium_amp_2 - Calcium_amp_1;
cell_setup = {};
cell_setup{1} = APD_difference;       %disp(APD_difference)
cell_setup{2} = Calcium_difference;   %disp(Calcium_difference)
cell_setup{3} = prot_rate;

