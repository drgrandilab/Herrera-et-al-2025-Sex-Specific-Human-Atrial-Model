function [output] = A0_Obtain_intial_conditions(duration, frequencey, AF, ISO_param, gender, first, second,Plot_Flag)

%% Parameters for external modules

if gender == 0 && AF == 0
    load yfin_nSR_1Hz_0_ISO_male.mat;
elseif gender == 1 && AF == 0
    load yfin_nSR_1Hz_0_ISO_female.mat;
elseif gender == 0 && AF == 1
    load yfin_cAF_1Hz_0_ISO_male.mat;
elseif gender == 1 && AF == 1
    load yfin_cAF_1Hz_0_ISO_female.mat;
end

load parameter_matrix_1000_all_parameters_addional_600_cells.mat % sigma 0.1
all_parameters = all_parameters;
[N_trials N_par] = size(all_parameters);
y0n = yfinal;
beat_analysis = 1;
%% Param
AF_index = AF;
freq = frequencey;
cycleLength = 1e3/freq;
ISO = ISO_param;
Currents_record = 1;
gender_flag = gender;
output_currents = 1;


%Prot_index
%protocol(1) = 'Paced';
%protocol(4) = 'DAD';

prot_index = 1;

%% Collect all parameters

% 1 = female; 0 = male (baseline)

Female_diff = sex_diff;

if gender_flag == 1
    Female_GNa        = Female_diff(1,1); % Ambrosi 2013; added
    Female_Gk1        = Female_diff(1,2); % Ambrosi 2013; added
    Female_Gkur       = Female_diff(1,3); % Ambrosi 2013; added
    Female_Gkach      = Female_diff(1,4); % Ambrosi 2013; added
    Female_Ito        = Female_diff(1,5); % Ambrosi 2013; added
    Female_CSQ        = Female_diff(1,6); % Madsen 2021 ; added
    Female_NaK        = Female_diff(1,7); %
    Female_SK         = Female_diff(1,8); %
elseif gender_flag == 0
    Female_GNa        = 1;
    Female_Gk1        = 1;
    Female_Gkur       = 1;
    Female_Gkach      = 1;
    Female_Ito        = 1;
    Female_CSQ        = 1;
    Female_NaK        = 1;
    Female_SK         = 1;

end

gender_param = [gender_flag Female_GNa Female_Gk1 Female_Gkur Female_Gkach...
    Female_Ito Female_CSQ Female_NaK Female_SK]; %


%% Establish and define globals

if Currents_record  == 1 || output_currents == 1

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

for ii = first:second
    %     y0n = all_ICs(ii,:);
    par_SA = all_parameters(ii,:); % 25 parameters
    disp(ii)
    p = [cycleLength, AF_index, prot_index, ISO, Currents_record, output_currents,...
        gender_param, par_SA];
    tspan = [0 duration]; % [ms]
    options = odeset('RelTol',1e-5,'MaxStep',1);
    [t,y] = ode15s(@NH_ODE,tspan,y0n,options,p);
    all_ICs(ii,:) = y(end,:)';
    period = 1000/prot_rate;

    if Plot_Flag == 1
        % Create the figure only on first iteration
        if ii == first
            fig = figure('Color','w','Position',[100 100 600 800], 'Name','Population traces');
            ax1 = subplot(3,1,1); hold(ax1,'on'); grid(ax1,'on');
            ylabel(ax1,'[Na]_i','FontSize',12,'FontWeight','bold');

            ax2 = subplot(3,1,2); hold(ax2,'on'); grid(ax2,'on');
            ylabel(ax2,'Vm','FontSize',12,'FontWeight','bold');

            ax3 = subplot(3,1,3); hold(ax3,'on'); grid(ax3,'on');
            ylabel(ax3,'[Ca]_i (nM)','FontSize',12,'FontWeight','bold');
            xlabel(ax3,'Time (ms)','FontSize',12,'FontWeight','bold');

            sgtitle('Simulation Outputs','FontSize',16,'FontWeight','bold');
            set(findall(fig,'-property','FontName'),'FontName','Arial');
        end

        % Plot into existing axes
        plot(ax1, t, y(:,34), 'LineWidth', 0.8);
        plot(ax2, t, y(:,39), 'LineWidth', 0.8);
        plot(ax3, t, y(:,38)*1e6, 'LineWidth', 0.8);

        % Force a UI refresh so you see progress as it runs
        drawnow limitrate
        % (Optional) pause(0.01)  % tiny pause if you still don’t see updates
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

        output(ii,:) = [dVm_max Vm_max RMP AP_amp APD90 APD70 APD50 APD30 Ca_max...
            Ca_min CaT_amp CaT_rise CaT_decay_50 CaT_decay_63 Na_min VPLT APD20];
    end
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

if output_currents == 1
    total_store_values = [dVm; ICa; Ito; INa; IKs; IKr; IKp; IClCa; IClbk; INCX; INaK;...
        INabk; IPMCA; ICabk; Cai_tA; Jserca; I_app; I_NaL; I_Kur; I_K2P;...
        I_Ki; I_Kach; I_SK; I_ClFTR; I_pca; IK1; I_J_SRCarel; tA]';

    total_store_values(:,29) = Jleak(:,1);
    total_store_values(:,30) = Jleak(:,2);
else
    total_store_values = 0;
end
%% Current Plot
if Currents_record == 1
    figure(2)
    set(gcf,'color','w')
    subplot(4,2,1); plot(tA,INa); title('INa');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(4,2,3); plot(tA,ICa); title('ICa');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(4,2,5); plot(tA,INCX); title('INCX');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(4,2,7); plot(tA,INaK); title('INaK');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(4,2,2); plot(tA,IKr); title('IKr');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(4,2,4); plot(tA,IKs); title('IKs');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(4,2,6); plot(tA,Ito); title('Ito');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(4,2,8); plot(tA,IK1); title('IK1');

    figure(3)
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,1); plot(tA,IKp); title('IKp');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,2); plot(tA,IClCa); title('IClCa');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,3); plot(tA,IClbk); title('IClbk');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,4); plot(tA,INabk); title('INabk');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,5); plot(tA,IPMCA); title('IPMCA');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,6); plot(tA,ICabk); title('ICabk');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,7); plot(tA,Cai_tA); title('Cai_tA');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,8); plot(tA,Jserca); title('Jserca');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,9); plot(tA,Jleak); title('Jleak');
    set(gca,'box','off','tickdir','out','fontsize',12)

    figure(4)
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,1); plot(tA,I_app); title('Iapp');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,2); plot(tA,I_NaL); title('INaL');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,3); plot(tA,I_Kur); title('IKur');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,4); plot(tA,I_K2P); title('IK2P');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,5); plot(tA,I_Ki); title('IKi');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,6); plot(tA,I_Kach); title('IKach');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,7); plot(tA,I_SK); title('ISK');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,8); plot(tA,I_ClFTR); title('IClFTR');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,9); plot(tA,I_pca); title('Ipca');
    set(gca,'box','off','tickdir','out','fontsize',12)

end


