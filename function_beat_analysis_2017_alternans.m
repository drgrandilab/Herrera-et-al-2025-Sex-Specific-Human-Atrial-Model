function outputs = function_beat_analysis_2017_alternans(time,Vm,Ca,period,AP_index, figures)
%% Check
%% Check
try
    if AP_index == 1,
        t1 = 0; t2 = t1+period-1.5; t3 = t1+2*period-1.5;
    else
        tt = time(end) - 4*period; t0 = time(end) - 3*period; t1 = time(end)- 2*period; t2 = t1+period-1.5; t3 = t1+2*period-1.5;
    end

    %% Alternans analysis
    t0_roi = find(time>t0); t0_index = t0_roi(1); Vm0 = Vm(t0_index);
    t1_roi = find(time>t1); t1_index = t1_roi(1); Vm1 = Vm(t1_index);
    t2_roi = find(time>t2); t2_index = t2_roi(1); Vm2 = Vm(t2_index);
    t3_roi = find(time>t3); t3_index = t3_roi(1); Vm3 = Vm(t3_index);

    max_1st = max(Vm(t0_index:t1_index)); index_max_1st = find (Vm == max_1st);
    max_2nd = max(Vm(t1_index:t2_index)); index_max_2nd = find (Vm == max_2nd);
    Vm_min_1 = min(Vm(index_max_1st:t1_index)); index_Vm_min = find (Vm == Vm_min_1);
    Vm_min_2 = min(Vm(index_max_2nd:t2_index)); index_Vm_min_2 = find (Vm == Vm_min_2);

    Vm_1_dvdt_setup = Vm(t0_index: t1_index); time_1 = time(t0_index: t1_index); dv_dt_1 = (diff(Vm_1_dvdt_setup))./(diff(time_1)); max_dv_dt_1 = max(dv_dt_1);
    time_dvdt = find (dv_dt_1 == max_dv_dt_1); index_time = time_1(time_dvdt); vm_dvdt_time = find(time == index_time); vm_dvdt_int = Vm(vm_dvdt_time);


    Vm_2_dvdt_setup = Vm(t1_index:t2_index); time_2 = time(t1_index: t2_index); dv_dt_2 = (diff(Vm_2_dvdt_setup))./(diff(time_2)); max_dv_dt_2 = max(dv_dt_2);
    time_dvdt_2 = find (dv_dt_2 == max_dv_dt_2); index_time_2 = time_2(time_dvdt_2); vm_dvdt_time_2 = find(time == index_time_2); vm_dvdt_int_2 = Vm(vm_dvdt_time_2);


    time_min_2 = time(index_Vm_min_2);
    % time_min_1 = time(index_Vm_min);
    % time_peak_2 = time(index_max_2nd);
    % time_peak_1 = time(index_max_1st);


    time_45 = Vm(index_max_1st:index_Vm_min); time_op = time(index_max_1st:index_Vm_min); time_45_point = find(Vm(time_45>-55)); %changed to -60 2/4/2023
    time_45_point_index = time_45_point(end); time_at_45 = time_op(time_45_point_index); vm_45_time = find(time == time_at_45); vm_45_int = Vm(vm_45_time);
    time_45_2 = Vm(index_max_2nd:index_Vm_min_2); time_op_2 = time(index_max_2nd:index_Vm_min_2); time_45_point_2 = find(Vm(time_45_2>-55)); %changed to -60 2/4/2023
    time_45_point_index_2 = time_45_point_2(end); time_at_45_2 = time_op_2(time_45_point_index_2); vm_45_time_2 = find(time == time_at_45_2); vm_45_int_2 = Vm(vm_45_time_2);
    vm_int = (Vm(t0_index));
    vm_int_2 = (Vm(t1_index));



    APD90_1 = (time_at_45-t0);
    APD90_2 = (time_at_45_2-t1);

    if APD90_1<APD90_2
        APD90_1 = (time_at_45_2-t1);
        APD90_2 = (time_at_45-t0);
    else
        APD90_1 = APD90_1;
        APD90_2= APD90_2;
    end

    APD_difference = APD90_1-APD90_2;
    APD_cutoff = 5;
    if APD_difference <APD_cutoff   % READ TO CHANGE THIS VALUE
        alternans_APD = false;
    else
        alternans_APD = true;
    end


    if figures == 1
        figure(1)
        subplot(3,1,1),
        hold on
        plot (time_at_45, vm_45_int, 'bx')
        plot (time_at_45_2, vm_45_int_2, 'rx')
        plot (t0, vm_int, 'bx')
        plot (t1, vm_int_2, 'rx')
    end



    %% Calcium ampltitude analysis
    [Ca_max, index_ca] = max(Ca(t0_index:t1_index)); index_ca = find (Ca == Ca_max); % peak CaT
    [Ca_max_2, index_ca_2] = max(Ca(t1_index:t2_index)); index_ca_2 = find (Ca == Ca_max_2); % peak CaT

    Peak_calcium = Ca(t1_index:t2_index); Peak_calcium_index = findpeaks(Peak_calcium);
    Peak_calcium_ca = max(Peak_calcium_index); Index_calcium_max = find(Ca == Peak_calcium_ca); time_peak_ca = time(Index_calcium_max);


    Ca_min_new = min(Ca(index_ca:t1_index)); index_ca_min = find(Ca == Ca_min_new); time_ca_min = time(index_ca_min);
    Ca_min_0 = min(Ca(t0_index:t1_index)); index_ca_min_0 = find (Ca == Ca_min_0); % diast Ca
    % Ca_min = min(Ca(t1_index:t2_index)); index_ca_min = find (Ca == Ca_min); % diast Ca

    time_min_0 = time(index_ca_min_0);
    time_peak = time(index_ca);
    time_peak_2 = time(index_ca_2);
    time_min = time (index_ca_min);


    Ca_duration_1 = (time_min-time_peak)*.80;
    Ca_duration_2 = (time_min_2-time_peak_2)*0.80;


    time_index_min = time(index_ca_min);
    time_index = time(index_ca);
    time_index_2 = time(index_ca_2);


    CaT_amp = Ca_max-Ca_min_0; %red
    CaT_amp_2 = Peak_calcium_ca-Ca_min_new; %blue
    %
    if CaT_amp_2<CaT_amp
        CaT_amp = (Ca_max_2-Ca_min_new);
        CaT_amp_2 = (Ca_max-Ca_min_0);
    else
        CaT_amp = CaT_amp;
        CaT_amp_2 = CaT_amp_2;
    end


    if figures == 1
        figure(1)
        subplot(3,1,2),
        hold on, plot (time_index, Ca_max*1000000, 'ro')
        plot (time_min_0,Ca_min_0*1000000, 'ro')
        plot (time_peak_ca, Peak_calcium_ca*1000000, 'bo')
        plot (time_ca_min, Ca_min_new*1000000, 'bo')
    end
    value_check = 1;
    outputs = [APD90_1 APD90_2 CaT_amp CaT_amp_2 Ca_duration_1 Ca_duration_2 Ca_min_new alternans_APD  APD_difference value_check];


catch


    outputs = [0 0 0 0 0 0 0 0 0 0];


end


