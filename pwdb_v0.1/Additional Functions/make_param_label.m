function [label, units, abbr, graph_title, graph_title_no_units] = make_param_label(curr_param)


switch curr_param
    case {'HR', 'hr'}
        label = 'HR: Heart rate [bpm]';
        units = 'bpm';
        abbr = 'HR';
        graph_title = 'Heart Rate [bpm]';
    case {'SV', 'sv'}
        label = 'SV: Stroke volume [ml]';
        units = 'ml';
        abbr = 'SV';
        graph_title = 'Stroke Volume [ml]';
    case 'CO'
        label = 'CO: Cardiac output [l/min]';
        units = 'l/min';
        abbr = 'CO';
        graph_title = 'Cardiac Output [l/min]';
    case {'LVET', 'lvet'}
        label = 'LVET: Left ventricular ejection time [ms]';
        units = 'ms';
        abbr = 'LVET';
        graph_title = 'Left Ventricular Ejection Time [ms]';
    case 'dPdt'
        label = 'dP/dt: Maximum aortic value [mmHg/s]';
        units = 'mmHg/s';
        abbr = 'dp/dt';
        graph_title = 'dP/dt [mmHg/s]';
    case {'t_PF', 't_pf', 'PFT', 'pft'}
        label = 'PFT: Peak flow time [ms]';
        units = 'ms';
        abbr = 'PFT';
        graph_title = 'Peak Flow Time [ms]';
    case {'reg_vol', 'RFV', 'rfv'}
        label = 'Reverse flow volume [ml]';
        units = 'ml';
        abbr = 'RFV';
        graph_title = 'Reverse Flow Volume [ml]';
    case 'SBP_a'
        label = 'Aortic pressure [mmHg]: SBP, systolic';
        units = 'mmHg';
        abbr = 'SBP_a';
        graph_title = 'Aortic SBP [mmHg]';
    case 'DBP_a'
        label = '\\hspace{4.1cm} DBP, diastolic';
        units = 'mmHg';
        abbr = 'DBP_a';
        graph_title = 'Aortic Disatolic Blood Pressure [mmHg]';
    case 'dbp'
        label = 'dbp';
        units = 'mmHg';
        abbr = 'DBP';
        graph_title = 'Diastolic Blood Pressure [mmHg]';  % 'DBP [mmHg]'
    case {'MBP_a', 'mbp', 'MBP'}
        label = '\\hspace{4.1cm} MAP, mean';
        units = 'mmHg';
        abbr = 'MAP_a';
        graph_title = 'Mean Arterial Pressure [mmHg]'; %'Aortic MBP [mmHg]';
    case 'PP_a'
        label = '\\hspace{4.1cm} PP, pulse pressure';
        units = 'mmHg';
        abbr = 'PP_a';
        graph_title = 'Aortic PP [mmHg]';
    case 'SBP_b'
        label = 'Brachial pressure [mmHg]: SBP';
        units = 'mmHg';
        abbr = 'SBP_b';
        graph_title = 'Brachial SBP [mmHg]';
    case 'DBP_b'
        label = '\\hspace{4.5cm} DBP';
        units = 'mmHg';
        abbr = 'DBP_b';
        graph_title = 'Brachial DBP [mmHg]';
    case 'MBP_b'
        label = '\\hspace{4.5cm} MBP';
        units = 'mmHg';
        abbr = 'MBP_b';
        graph_title = 'Brachial MBP [mmHg]';
    case 'PP_b'
        label = '\\hspace{4.5cm} PP';
        units = 'mmHg';
        abbr = 'PP_b';
        graph_title = 'Brachial PP [mmHg]';
    case 'SBP_f'
        label = 'Digital pressure [mmHg]: SBP';
        units = 'mmHg';
        abbr = 'SBP_b';
        graph_title = 'Digital SBP [mmHg]';
    case 'DBP_f'
        label = '\\hspace{4.5cm} DBP';
        units = 'mmHg';
        abbr = 'DBP_b';
        graph_title = 'Digital DBP [mmHg]';
    case 'MBP_f'
        label = '\\hspace{4.5cm} MBP';
        units = 'mmHg';
        abbr = 'MBP_b';
        graph_title = 'Digital MBP [mmHg]';
    case 'PP_f'
        label = '\\hspace{4.5cm} PP';
        units = 'mmHg';
        abbr = 'PP_b';
        graph_title = 'Digital PP [mmHg]';
    case 'PP_amp'
        label = 'Pulse Pressure Amplification';
        units = 'ratio';
        abbr = 'PP_{amp}';
        graph_title = 'Pulse Pressure Amplification';
    case {'AP', 'AP_c', 'AP_a'}
        label = 'Augmentation Pressure [mmHg]';
        units = 'mmHg';
        abbr = 'AP';
        graph_title = 'Augmentation Pressure [mmHg]';
    case {'AIx', 'AI_c', 'AI_a'}
        label = 'Augmentation Index [\\%%]';
        units = '%%';
        abbr = 'AIx';
        graph_title = 'Augmentation Index [%]';
    case 'PWV_a'
        label = 'Pulse wave velocity [m/s]: aortic';
        units = 'm/s';
        abbr = 'PWV_{a}';
        graph_title = 'Aortic PWV [m/s]';
    case {'PWV_cf', 'pwv_cf', 'expected_pwv_aorta'}
        label = '\\hspace{4.43cm} carotid-femoral';
        units = 'm/s';
        abbr = 'PWV_{cf}';
        graph_title = 'Carotid-Femoral PWV [m/s]';
    case 'pwv'
        label = 'Pulse wave velocity [m/s]';
        units = 'm/s';
        abbr = 'PWV';
        graph_title = 'Pulse wave velocity [m/s]';
    case {'PWV_br', 'pwv_br', 'expected_pwv_arm'}
        label = '\\hspace{4.43cm} brachial-radial';
        units = 'm/s';
        abbr = 'PWV_{br}';
        graph_title = 'Brachial-Radial PWV [m/s]';
    case {'PWV_fa', 'pwv_fa', 'expected_pwv_leg'}
        label = '\\hspace{4.43cm} femoral-ankle';
        units = 'm/s';
        abbr = 'PWV_{fa}';
        graph_title = 'Femoral-Ankle PWV [m/s]';
    case {'dia_asc_a', 'dia_asc'}
        label = 'Diameter [mm]: ascending aorta';
        units = 'mm';
        abbr = 'dia_{asca}';
        graph_title = 'Asc. Aorta Dia [mm]';
    case {'dia_desc_thor_a', 'dia_desc_thor'}
        label = '\\hspace{2.7cm} descending thoracic aorta';
        units = 'mm';
        abbr = 'dia_{dta}';
        graph_title = 'Desc. Thor. Aorta Dia [mm]';
    case {'dia_abd_a', 'dia_abd'}
        label = '\\hspace{2.7cm} abdominal aorta';
        units = 'mm';
        abbr = 'dia_{abda}';
        graph_title = 'Abd. Aorta Dia [mm]';
    case {'dia_car', 'dia_carotid'}
        label = '\\hspace{2.7cm} carotid';
        units = 'mm';
        abbr = 'dia_{car}';
        graph_title = 'Carotid Dia [mm]';
    case 'len_prox_a'
        label = 'Length of proximal aorta [mm]';	
        units = 'mm';	
        abbr = 'Len';
        graph_title = 'Prox. Aortic Length [mm]';
    case 'svr'
        label = 'Systemic vascular resistance [10$^6$ Pa s m$^{-3}$]';	
        units = '10^6 Pa s / m3';	
        abbr = 'SVR';
        graph_title = 'Systemic Vascular Resistance';
    case {'Tr', 'Tr_a', 'Tr_c'}
        label = 'Time to Reflected Wave [ms]';	
        units = 'ms';	
        abbr = 'Tr';
        graph_title = 'Time to Reflected Wave [ms]';	
    case 'ht'
        label = 'Height [m]';	
        units = 'm';	
        abbr = 'Ht';
        graph_title = 'Height [m]';	
    case 'MBP_drop_finger'
        label = 'Pressure drop to finger [mmHg]';
        units = 'mmHg';
        abbr = 'drop fin';
        graph_title = 'Pressure drop to finger [mmHg]';
    case 'MBP_drop_ankle'
        label = 'Pressure drop to ankle [mmHg]';
        units = 'mmHg';
        abbr = 'drop ankle';
        graph_title = 'Pressure drop to ankle [mmHg]';
    case 'dia'
        label = 'Arterial diameter';	
        units = 'm';	
        abbr = 'Dia';
        graph_title = 'Arterial diameter [m]';
    case 'len'
        label = 'Proximal aortic length [m]';	
        units = 'm';	
        abbr = 'Len';
        graph_title = 'Proximal Aortic Length [m]';
    case {'rho','density'}
        label = 'Density';	
        units = 'kg /m^3';	
        abbr = 'rho';
        graph_title = 'Density [kg /m^3]';
    case {'pvc', 'PVC'}
        label = 'Peripheral vasc comp.';	
        units = 'scaling factor';	
        abbr = 'PVC';
        graph_title = 'Peripheral Vasc Comp. []';
    case {'pvc_adj'}
        label = 'Peripheral vasc compliance [10^9 m^3/Pa]';	
        units = 'm3/Pa';	
        abbr = 'PVC';
        graph_title = 'Peripheral vasc comp. [10^9 m^3/Pa]';
    case 'p_out'
        label = 'Outflow pressure [mmHg]';
        units = 'mmHg';
        abbr = 'p_out';
        graph_title = 'Outflow pressure [mmHg]';
    case 'p_drop'
        label = 'Assumed pressure drop [mmHg]';
        units = 'mmHg';
        abbr = 'p_drop';
        graph_title = 'Assumed pressure drop [mmHg]';
    case 'time_step'
        label = 'Simulation time step [s]';
        units = 's';
        abbr = 'time_step';
        graph_title = 'Simulation time step [s]';
    case {'mu', 'viscosity'}
        label = 'Viscosity [Pa s]';
        units = 'Pa s';
        abbr = 'mu';
        graph_title = 'Viscosity [Pa s]';
    case 'alpha'
        label = 'Alpha [-]';
        units = '-';
        abbr = '\alpha';
        graph_title = 'Alpha [-]';
    case 'age'
        label = 'Age [years]';
        units = 'years';
        abbr = 'age';
        graph_title = 'Age [years]';
    case 'visco_elastic_log'
        label = 'Elastic (0) or visco-elastic (1)';
        units = '';
        abbr = 'viscLog';
        graph_title = 'Elastic vs Visco-Elastic';
    case 'AGI_mod'
        label = 'Modified Ageing Index [au]';
        units = 'au';
        abbr = 'agi_mod';
        graph_title = 'Modified Ageing Index';
    case 'gamma'
        label = 'Wall viscosity';
        units = 'au';
        abbr = 'gamma';
        graph_title = 'Wall Viscosity';
    case 'RI'
        label = 'Reflection Index [au]';
        units = 'au';
        abbr = 'ri';
        graph_title = 'Reflection Index';
    case 'SI'
        label = 'Stiffness Index [m/s]';
        units = 'm/s';
        abbr = 'ri';
        graph_title = 'Stiffness Index';
    case 'IAD'
        label = 'Inter-arm SBP Difference';
        units = 'mmHg';
        abbr = 'IAD';
        graph_title = 'Inter-arm SBP Difference';
    case 'IADabs'
        label = 'abs Inter-arm SBP Difference';
        units = 'mmHg';
        abbr = 'IADabs';
        graph_title = 'abs Inter-arm SBP Difference';
    case 'pvr'
        label = 'Peripheral Vascular Resistance';
        units = 'Pa s/m^3';
        abbr = 'PVR';
        graph_title = label;
    case {'b0', 'gamma_b0'}
        label = 'Wall Viscosity b0 Parameter';
        units = 'g/s';
        abbr = 'b0';
        graph_title = label;
    case {'b1', 'gamma_b1'}
        label = 'Wall Viscosity b1 Parameter';
        units = 'g cm/s';
        abbr = 'b1';
        graph_title = label;
    case 'k1'
        label = 'Stiffness k1 Parameter';
        units = 'g/s^2/cm';
        abbr = 'k1';
        graph_title = label;
    case 'k2'
        label = 'Stiffness k2 Parameter';
        units = '/cm';
        abbr = 'k2';
        graph_title = label;
    case 'k3'
        label = 'Stiffness k3 Parameter';
        units = 'g/s^2/cm';
        abbr = 'k3';
        graph_title = label;
    case 'reflect_log'
        label = 'Reflections';
        units = '';
        abbr = 'Reflections';
        graph_title = label;
    case 'SBP_diff'
        label = 'Difference in aortic and brachial SBP [mmHg]';
        units = 'mmHg';
        abbr = 'SBP_diff';
        graph_title = 'Difference in aortic and brachial SBP [mmHg]';
    case 'SMBP_a'
        label = 'Systolic - Mean Aortic BP [mmHg]';
        units = 'mmHg';
        abbr = 'SMBP';
        graph_title = 'Systolic - Mean Aortic BP [mmHg]';
    case 'tau'
        label = 'Time constant [s]';
        units = 's';
        abbr = 'tau';
        graph_title = 'Time constant [s]';
    case 'P1t_a'
        label = 'Time of aortic P1';
        units = 's';
        abbr = 'P1t';
        graph_title = 'Time of aortic P1 [s]';
    case 'P1t_c'
        label = 'Time of carotid P1';
        units = 's';
        abbr = 'P1t';
        graph_title = 'Time of carotid P1 [s]';
    case 'P1_a'
        label = 'Magnitude of aortic P1';
        units = 'mmHg';
        abbr = 'P1';
        graph_title = 'Magnitude of aortic P1 [mmHg]';
    case 'P1_c'
        label = 'Magnitude of carotid P1';
        units = 'mmHg';
        abbr = 'P1';
        graph_title = 'Magnitude of carotid P1 [mmHg]';
    case 'P2t_a'
        label = 'Time of aortic P2';
        units = 's';
        abbr = 'P2t';
        graph_title = 'Time of aortic P2 [s]';
    case 'P2t_c'
        label = 'Time of carotid P2';
        units = 's';
        abbr = 'P2t';
        graph_title = 'Time of carotid P2 [s]';
    case 'P2_a'
        label = 'Magnitude of aortic P2';
        units = 'mmHg';
        abbr = 'P2';
        graph_title = 'Magnitude of aortic P2 [mmHg]';
    case 'P2_c'
        label = 'Magnitude of carotid P2';
        units = 'mmHg';
        abbr = 'P2';
        graph_title = 'Magnitude of carotid P2 [mmHg]';
    case 'trial'
        label = 'Trial Parameter';
        units = '';
        abbr = 'trial';
        graph_title = 'Trial Parameter []';
    case 'trial2'
        label = 'Trial Parameter 2';
        units = '';
        abbr = 'trial';
        graph_title = 'Trial Parameter 2 []';
    case 'trial3'
        label = 'Trial Parameter 3';
        units = '';
        abbr = 'trial';
        graph_title = 'Trial Parameter 3 []';
end

% create graph title without units
graph_title_no_units = graph_title;
temp = strfind(graph_title_no_units, ' [');
if ~isempty(temp)
    graph_title_no_units = graph_title_no_units(1:temp-1);
end
graph_title_no_units = strrep(graph_title_no_units, 'MBP', 'Mean Arterial Pressure');
graph_title_no_units = strrep(graph_title_no_units, 'SBP', 'Systolic Blood Pressure');
graph_title_no_units = strrep(graph_title_no_units, 'DBP', 'Diastolic Blood Pressure');
graph_title_no_units = strrep(graph_title_no_units, 'PP', 'Pulse Pressure');

end
