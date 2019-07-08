function co_case_study(pwdb_no)
% CO_CASE_STUDY generates the plots reported in the case study on
% estimating cardiac output from blood pressure waveforms
%
%               co_case_study
%
%   Inputs:     - the 'pwdb_data' file produced by 'export_pwdb.m'.
%
%   Outputs:    - plots illustrating the results of the case study
%           
%   Accompanying Article:
%       This code is provided to facilitate reproduction of the Pulse Wave
%       Database described in: 
%           Charlton P.H. et al. Modelling arterial pulse waves in healthy
%           ageing: a database for in silico evaluation of haemodynamics
%           and pulse wave indices, ~~ under review ~~  
%           DOI: ~~ tbc ~~
%       Further information on the pwdb Pulse Wave Database is provided in
%       this article and at the project website:
%           https://peterhcharlton.github.io/pwdb
%
%   Author: Peter H. Charlton
%
%   Licence:
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%  Copyright (C) 2019  King's College London
%
% v.1.0  Peter H. Charlton, King's College London

fprintf('\n --- Running Cardiac Output (CO) Case Study ---')

% Setup paths with current simulation paths
PATHS = setup_paths_for_post_processing(pwdb_no);

% Cardiac output algorithm settings
up.settings = cardiac_output_settings;
display_settings(up);

% Setup universal parameters
up = setup_up(up);

% Collate data
data = load_data(PATHS, up);

% Find cardiac output
data = find_cardiac_output(data, up);

% Make results plots
make_results_plots(data, up, PATHS);

% Make additional plots
if up.make_plots
    make_additional_plots(co_data, up)
end

% % Display results
% display_results(co_data, up);

end

function settings = cardiac_output_settings

% This specifies the settings used for cardiac output estimation

%% %%%%%%%%%%%%%%%%% Anatomical Site for ABP measurement %%%%%%%%%%%%%%%%%
settings.artery.options = { ...
    'Femoral', ...      % 1
    'Radial', ...       % 2
    'AorticRoot', ...   % 3
    'Brachial', ...     % 4
    };  
settings.artery.choice = 2;

end

function display_settings(up)

% This displays the settings that the user has chosen

different_settings = fieldnames(up.settings);
desired_len = 10;

message = '\n     ---------\n     Settings:\n     ---------';
for s = 1 : length(different_settings)
    curr_setting_name = strrep(different_settings{s}, '_', ' ');
    eval(['curr_setting_choice = up.settings.' different_settings{s} '.options{up.settings.' different_settings{s} '.choice};']);
    spaces = repmat(' ', [1, desired_len-length(curr_setting_name)]);
    message = [message, '\n     ', curr_setting_name, ':', spaces, curr_setting_choice];
    clear curr_setting*
end
clear s

message = [message, '\n'];

fprintf(message)


end

function up = setup_up(up)

close all
% identify the artery of interest
chosen_option = identify_chosen_option(up, 'artery');

% physiological parameters which are varied in the database
up.params = {'pwv', 'dia', 'mbp', 'lvet', 'hr', 'sv'};

% CO algorithms
up.algorithm_types =        {...   
%                                    'area under curve', ...       % 1  (pub 5)
%                                    'systolic area', ...          % 2  (pub 6?)
                                   'rms', ...                    % 3  (pub 10)
%                                    'rms_pulsatile', ...          % 4  (pub 11)
%                                    'integrate compliance', ...   % 5  (pub 8)
%                                    'constant compliance', ...    % 6  (pub 9)
%                                    'esvb', ...                   % 7
%                                    'erlanger_hooker', ...        % 8
                                   'liljestrand', ...            % 9 
%                                    'herd',...                    % 10
%                                    'harley',...                  % 11
%                                    'kouchoukos',...              % 12
%                                    'wesseling',...               % 13
%                                    'erlanger',...                % 14
%                                    'wesseling2',...              % 15
                                   };

%% - make additional plots
up.make_plots = 0;

%% - set current directory
cd(fileparts(which(mfilename)));

end

function data = load_data(PATHS, up)

% identify the artery of interest
chosen_option = identify_chosen_option(up, 'artery');

fprintf(['\n - Loading data for ' chosen_option ' artery'])

% Load all data
load(PATHS.exported_data_mat_pwdb_data);
orig_data = data; 
% Identify relevant data
rel_subjs = data.plausibility.plausibility_log';
clear data

% - set up groups
rel_cols = strcmp(orig_data.config.variations.param_names, 'pwv') | ... 
    strcmp(orig_data.config.variations.param_names, 'mbp') | ... 
    strcmp(orig_data.config.variations.param_names, 'dia');
% Group 1: changes in these variables, and no change in other variables
data.group1 = sum(abs(orig_data.config.variations.params(:,~rel_cols))')==0 & orig_data.config.baseline_sim_for_age' == 0;
data.group1 = data.group1(:);
% Group 2: changes in other variables, but no change in specified variables
data.group2 = ~data.group1 & orig_data.config.baseline_sim_for_age == 0;

% - set up groups
rel_cols = strcmp(orig_data.config.variations.param_names, 'pwv');
% Group 1: changes in these variables, and no change in other variables
data.group1 = sum(abs(orig_data.config.variations.params(:,~rel_cols))')==0 & orig_data.config.baseline_sim_for_age' == 0;
data.group1 = data.group1(:);

rel_cols = strcmp(orig_data.config.variations.param_names, 'dia');
% Group 2: changes in these variables, and no change in other variables
data.group2 = sum(abs(orig_data.config.variations.params(:,~rel_cols))')==0 & orig_data.config.baseline_sim_for_age' == 0;
data.group2 = data.group2(:);

rel_cols = strcmp(orig_data.config.variations.param_names, 'mbp');
% Group 3: changes in these variables, and no change in other variables
data.group3 = sum(abs(orig_data.config.variations.params(:,~rel_cols))')==0 & orig_data.config.baseline_sim_for_age' == 0;
data.group3 = data.group3(:);

rel_cols = strcmp(orig_data.config.variations.param_names, 'lvet');
% Group 4: changes in these variables, and no change in other variables
data.group4 = sum(abs(orig_data.config.variations.params(:,~rel_cols))')==0 & orig_data.config.baseline_sim_for_age' == 0;
data.group4 = data.group4(:);

rel_cols = strcmp(orig_data.config.variations.param_names, 'hr');
% Group 5: changes in these variables, and no change in other variables
data.group5 = sum(abs(orig_data.config.variations.params(:,~rel_cols))')==0 & orig_data.config.baseline_sim_for_age' == 0;
data.group5 = data.group5(:);

rel_cols = strcmp(orig_data.config.variations.param_names, 'sv');
% Group 6: changes in these variables, and no change in other variables
data.group6 = sum(abs(orig_data.config.variations.params(:,~rel_cols))')==0 & orig_data.config.baseline_sim_for_age' == 0;
data.group6 = data.group6(:);

% Group 7: all except baseline
data.group7 = orig_data.config.baseline_sim_for_age == 0;

% - extract data
eval(['data.p = orig_data.waves.P_' chosen_option '(:);']);
data.age = orig_data.config.age;
data.ref_log = orig_data.config.baseline_sim_for_age;
data.end_sys_samp_no = orig_data.waves.fs*extractfield(orig_data.haemods, 'LVET')/1000; data.end_sys_samp_no = data.end_sys_samp_no(:);
data.settings.fs = orig_data.waves.fs;

% - extract CO
data.ref_co = extractfield(orig_data.haemods, 'CO');
data.ref_co = data.ref_co(:);

% - identify calibration subject
data.cal_subj = nan(length(data.age),1);
for s = 1 : length(data.age)
    curr_age = data.age(s);
    data.cal_subj(s,1) = find(data.ref_log & data.age == curr_age);
end

% Extract param values
for param_no = 1 : length(up.params)
    curr_param = up.params{param_no};
    eval(['data.params.' curr_param ' = orig_data.config.' curr_param '_SD(rel_subjs);'])
end

% Remove physiologically implausible subjs from analysis
fields = fieldnames(data);
for field_no = 1 : length(fields)
    eval(['curr_field = data.' fields{field_no} ';']);
    % skip if this doesn't need adjusting
    if isstruct(curr_field)
        continue
    end
    % otherwise remove those subjects which were implausible
    curr_field = curr_field(rel_subjs);
    eval(['data.' fields{field_no} ' = curr_field;']);
end

end

function data = find_cardiac_output(data, up)

% algorithm to use
fprintf(['\n - Estimating cardiac output'])

co_data.sv_est = nan(length(data.age),length(up.algorithm_types));


% cycle through each subject
for subj_no = 1 : length(data.age)
    
    % cycle through each algorithm
    for alg_no = 1 : length(up.algorithm_types)
        
        curr_alg = up.algorithm_types{alg_no};
        
        % calculate uncalibrated sv
        switch curr_alg
            case 'area under curve'   % (charlton 2014, eqn 5)
                data.sv_est(subj_no, alg_no) = sum(data.p{subj_no})/data.settings.fs;
                
            case 'systolic area'    % (charlton 2014, eqn 6?)
                rel_t = 1:data.end_sys_samp_no(subj_no);
                data.sv_est(subj_no, alg_no) = sum(data.p{subj_no}(rel_t))/data.settings.fs;
                
            case 'rms'    % (charlton 2014, eqn 10)
                data.sv_est(subj_no, alg_no) = sqrt(mean(data.p{subj_no}.^2));
                
            case 'rms_pulsatile'    % (charlton 2014, eqn 11)
                data.sv_est(subj_no, alg_no) = sqrt(mean((data.p{subj_no}-mean(data.p{subj_no})).^2));
                
            case 'integrate compliance'       % (charlton 2014, eqn 8)
                start_dia_el = data.end_sys_samp_no(subj_no);
                compliance = calc_compliance(data.p{subj_no}, data.settings.fs);
                data.sv_est(subj_no, alg_no) = -1 * compliance * ...
                    (data.p{subj_no}(end) - data.p{subj_no}(start_dia_el) ) ...
                    * (1 + ( ...
                    (sum(data.p{subj_no}(1:start_dia_el))/data.settings.fs)/ ...   % assume Pout is 0
                    (sum(data.p{subj_no}(start_dia_el:end))/data.settings.fs) ...
                    ) );
                
            case 'constant compliance'    % (charlton 2014, eqn 9)
                compliance = calc_compliance(data.p{subj_no}, data.settings.fs);
                p_out = 0;
                [p_s, sys_peak_el] = max(data.p{subj_no});
                t_s = sys_peak_el/data.settings.fs;
                t_d = length(data.p{subj_no})/data.settings.fs;
                p_d = data.p{subj_no}(end);
                p_0 = data.p{subj_no}(1);
                tau = -1*(t_d-t_s)/log(p_d/p_s); % use D1
                % assuming P(t_d) = P(t_0)
                data.sv_est(subj_no, alg_no) = compliance*( (p_d-p_0) + ( (1/tau)*sum(data.p{subj_no}-p_out)/data.settings.fs ) );
                
            case 'esvb'   %?\cite{Papaioannou2012}
                a1 = 0.4;   % ?\cite{Reymond2009}
                b1 = 5;
                pmaxc = 20;
                pwidth = 30;
                C = a1 + (b1./ (1 + ( (data.p{subj_no}-pmaxc)/pwidth ).^2 ) );
                T = (length(data.p{subj_no})-1)/data.settings.fs;
                PP = max(data.p{subj_no}) - min(data.p{subj_no});
                data.sv_est(subj_no, alg_no) = mean(C*PP./T);
                
            case 'erlanger_hooker'   %?\cite{Papaioannou2012}
                PP = max(data.p{subj_no}) - min(data.p{subj_no});
                T = (length(data.p{subj_no})-1)/data.settings.fs;
                data.sv_est(subj_no, alg_no) = PP/T;
                
            case 'liljestrand'   %?\cite{Papaioannou2012}
                PP = max(data.p{subj_no}) - min(data.p{subj_no});
                T = (length(data.p{subj_no})-1)/data.settings.fs;
                p_d = data.p{subj_no}(end);
                p_s = max(data.p{subj_no});
                data.sv_est(subj_no, alg_no) = PP/(T*(p_s+p_d));
                
            case 'herd'   %?\cite{Papaioannou2012}
                p_d = data.p{subj_no}(end);
                p_m = mean(data.p{subj_no});
                T = (length(data.p{subj_no})-1)/data.settings.fs;
                data.sv_est(subj_no, alg_no) = (p_m - p_d)/T;
                
            case 'harley'   %?\cite{Papaioannou2012}
                PP = max(data.p{subj_no}) - min(data.p{subj_no});
                start_dia_el = data.end_sys_samp_no(subj_no);
                Ts = start_dia_el/data.settings.fs;
                T = (length(data.p{subj_no})-1)/data.settings.fs;
                data.sv_est(subj_no, alg_no) = PP*Ts/T;
                
            case 'kouchoukos'   %?\cite{Papaioannou2012}
                Psa = sum(data.p{subj_no}(rel_t))/data.settings.fs;
                start_dia_el = data.end_sys_samp_no(subj_no);
                Ts = start_dia_el/data.settings.fs;
                T = (length(data.p{subj_no})-1)/data.settings.fs;
                Td = T - Ts;
                data.sv_est(subj_no, alg_no) = (Psa/T)*(1 + (Ts/Td));
                
            case 'wesseling'   %?\cite{Papaioannou2012}
                Psa = sum(data.p{subj_no}(rel_t))/data.settings.fs;
                T = (length(data.p{subj_no})-1)/data.settings.fs;
                data.sv_est(subj_no, alg_no) = Psa/T;
                
            case 'erlanger'   %?\cite{Papaioannou2012}
                PP = max(data.p{subj_no}) - min(data.p{subj_no});
                data.sv_est(subj_no, alg_no) = PP;
                
            case 'wesseling2'   %?\cite{Papaioannou2012}
                Psa = sum(data.p{subj_no}(rel_t))/data.settings.fs;
                T = (length(data.p{subj_no})-1)/data.settings.fs;
                HR = 60/T;
                Pm = mean(data.p{subj_no});
                data.sv_est(subj_no, alg_no) = (Psa/T)*(163+HR-0.48*Pm);
        end
        
        % calibrate sv
        cal_const = data.sv_est(data.cal_subj(subj_no), alg_no)./data.ref_co(subj_no);
        data.co_cal(subj_no, alg_no) = data.sv_est(subj_no, alg_no)/cal_const;
        
    end
    
end

end

function compliance = calc_compliance(p, fs)

% using C4

dt = 1/fs;
p_out = 0;
compliance = 1./(sum(p-p_out)*dt);

% using C2

p_sys = max(p);
p_dia = p(end);
compliance = 1/(p_sys+p_dia);

% using C3
p_sys = max(p);
p_dia = p(end);
compliance = 1/(p_sys+(2*p_dia));



end

function chosen_option = identify_chosen_option(up, rel_setting)

eval(['rel_options = up.settings.' rel_setting '.options;']);
eval(['rel_choice = up.settings.' rel_setting '.choice;']);
chosen_option = rel_options{rel_choice};

end

function display_results(co_data, up)

d = co_data;

for subj_no = 1 : length(d)
    
    % identify the parameter of interest
    chosen_option = identify_chosen_option(up, 'physiological_parameter');
    eval(['param_data = co_data(subj_no).' chosen_option '.v;'])
    param_vals = unique(param_data);
    baseline_val = param_data(co_data(subj_no).baseline_el);
    baseline_el = find(param_vals == 0);
    init_color = 0.6;
    for s = 1 : length(param_vals)
        curr_param_val = param_vals(s);
        if curr_param_val == 0
            colors(s,:) = init_color*[0,0,1];
        elseif curr_param_val < 0
            const = 0.3/sum(param_vals<0);
            colors(s,:) = (init_color + (s-baseline_el)*const)*[0,0,1];
        else
            const = 0.4/sum(param_vals>0);
            colors(s,:) = (init_color + (s-baseline_el)*const)*[0,0,1];
        end
    end
    
    % - make figure
    screensize = get( 0, 'Screensize' );
    figure('Position', [125, 125, 900, 500])
    subplot(1,3,1:2)
    
    ftsize = 14;
    
    % Plot baseline value
    %plot(co_data.ref_co(co_data.baseline_el), co_data.co_est(co_data.baseline_el), 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','b','MarkerFaceColor','b'),
    %hold on,
    
    % Plot other values
    for s = 1 : length(co_data(subj_no).ref_co)
        color_number = find(param_vals == param_data(s));
        if s == find(co_data(subj_no).baseline_el)
            label_h(color_number) = plot(co_data(subj_no).ref_co(s), co_data(subj_no).co_est(s), 'o', 'color', colors(color_number,:), ...
                'MarkerSize', 10, 'MarkerEdgeColor', colors(color_number,:), 'MarkerFaceColor', colors(color_number,:));
            hold on
        else
            label_h(color_number) = plot(co_data(subj_no).ref_co(s), co_data(subj_no).co_est(s), 'o', 'color', colors(color_number,:), ...
                'MarkerSize', 6, 'MarkerEdgeColor', colors(color_number,:), 'MarkerFaceColor', colors(color_number,:));
            hold on
        end
    end
    
    % set limits and plot line of identity
    const = 0.2;
    lims = [min([co_data(subj_no).ref_co; co_data(subj_no).co_est]), max([co_data(subj_no).ref_co; co_data(subj_no).co_est])];
    lims = [floor(lims(1)-(const*range(lims))), ceil(lims(2)+(const*range(lims)))];
    xlim(lims), ylim(lims), axis equal, plot(lims, lims, 'k'), xlim(lims), ylim(lims)
    
    % Find individual errors
    mape = nan(size(param_vals));
    for s = 1 : length(mape)
        rel_els = find(param_data == param_vals(s));
        mape(s) = round(1000*mean(abs((co_data(subj_no).co_est(rel_els) - co_data(subj_no).ref_co(rel_els))./co_data(subj_no).ref_co(rel_els))))/10;
    end
    
    
    % make legend
    for s = 1 : length(param_vals)
        labels{s} = [chosen_option ' = ' num2str(param_vals(s)) ', MAPE = ' num2str(mape(s))];
    end
    legend(label_h', labels, 'Position', [0.65,0.5, 0.3, 0.2])
    
    % Add axis labels
    xlabel('Reference CO [l/min]', 'FontSize', ftsize)
    ylabel('Estimated CO [l/min]', 'FontSize', ftsize)
    set(gca, 'FontSize', ftsize, 'XTick', ceil(lims(1)):lims(end), 'YTick', ceil(lims(1)):lims(end))
    
    % Find error
    mape = round(10*100*mean(abs(co_data(subj_no).co_est - co_data(subj_no).ref_co)./co_data(subj_no).ref_co))/10;
    title(['Varying ' strrep(chosen_option, '_', ' ') '. MAPE = ' num2str(mape) ' %'])
    
    % Display results
    
    message = ['\n\nResults:', '\n--------'];
    message = [message, '\n     ', 'Mean absolute percentage error: ', num2str(mape), ' %%'];
    message = [message, '\n'];
    
    fprintf(message)
    
end

end

function make_additional_plots(co_data, up)

%% - Plot of ABP pulse

% extract data
data.v = co_data.p{find(co_data.baseline_el)};
data.t = (0:length(data.v)-1)./data.settings.fs;

% make figure
figure('Position', [100,100,500,400])
ftsize = 14;

% plot areas
rel_els = data.t <= up.analysis.duration_systole/1000;
h = area(data.t(rel_els), data.v(rel_els)); hold on
h.FaceColor = 0.7*ones(1,3);
h.LineStyle = 'none';
dim = [.23 .3 .1 .1];
str = {'Systolic','Area'};
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'FontSize', ftsize);
rel_els = data.t > up.analysis.duration_systole/1000;
h = area(data.t(rel_els), data.v(rel_els)); hold on
h.FaceColor = 0.5*ones(1,3);
h.LineStyle = 'none';
dim = [.53 .3 .1 .1];
str = {'Diastolic','Area'};
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'FontSize', ftsize);

% plot pulse
lwidth = 2; 
plot(data.t,data.v, 'b', 'LineWidth', lwidth)
xlabel('Time   [s]', 'FontSize', ftsize)
ylabel('ABP   [Pa]', 'FontSize', ftsize)
set(gca, 'FontSize', ftsize)
ylim([0, 1.1*max(data.v)])
box off

savepath = 'C:\Users\pc13\Dropbox\Work\AppliedMaths_SummerSchool\Classes\Day7\Group Activity\Figures\abp_pulse';
print(gcf,savepath,'-depsc')
close('all')

end

function make_results_plots(data, up, PATHS)

fprintf('\n - Making plots')

% calculate errors
data.error = data.co_cal-data.ref_co;
data.ape = abs(100*(data.error./data.ref_co));

% cycle through each algorithm
for rel_alg = 1:length(data.co_cal(1,:))
    
    % Extract results for this algorithm
    params = fieldnames(data.params);
    rel_ape_data = nan(sum(data.group1),length(params)); %sum(data.group1),length(params));
    for param_no = 1 : length(params)
        curr_param = params{param_no};
        eval(['rel_subjs = data.group' num2str(param_no) ';']);
        rel_ape_data(:, param_no) = data.ape(rel_subjs,rel_alg);
    end
    
    %% make scatter plots
    
    % setup
    paper_size = [500,400];
    figure('Position', [20,20,paper_size])
    ftsize = 22;
    all_color = 0.4*[1,1,1];
    
    % plot
    plot([0,14], [0,14], 'Color', 0.4*[1,1,1]), hold on
    x_data = data.ref_co(data.group7);
    y_data = data.co_cal(data.group7,rel_alg);
    plot(x_data, y_data, 'o', 'Color', all_color, 'MarkerFaceColor', all_color)
    
    x_data = data.ref_co(data.group3);
    y_data = data.co_cal(data.group3,rel_alg);
    plot(x_data, y_data, 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 8)
    
    x_data = data.ref_co(data.group6 | data.group5);
    y_data = data.co_cal(data.group6 | data.group5,rel_alg);
    plot(x_data, y_data, 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 8)
    
    
    % tidy up
    all_data = [data.ref_co(data.group7); data.co_cal(:,rel_alg)];
    lims = [floor(min(all_data)), ceil(max(all_data))];
    xlim(lims), ylim(lims)
    box off
    set(gca, 'FontSize', ftsize)
    xlabel('Reference CO (l/min)', 'FontSize', ftsize)
    ylabel('Estimated CO (l/min)', 'FontSize', ftsize)
    curr_alg = up.algorithm_types{rel_alg};
    switch curr_alg
        case 'liljestrand'
            title_text = 'Pulse Pressure';
        case 'rms'
            title_text = 'Root-Mean-Square';
    end
    dim = [0.15,0.7,0.7,0.3];
    annotation('textbox',dim,'String',[title_text, ' Algorithm'],'LineStyle', 'none', 'FontSize', ftsize, 'HorizontalAlignment', 'center');
    
    % annotate MAPE
    data.mape(rel_alg,1) = mean(data.ape(data.group7, rel_alg));
    dim = [.5 .14 .2 .2];
    str = ['MAPE (all) = ' num2str(data.mape(rel_alg),'%5.1f'), ' %'];
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'FontSize', ftsize-2, 'Color', all_color);
    
    data.mape_mbp(rel_alg,1) = mean(data.ape(data.group3, rel_alg));
    dim = [.5 .08 .2 .2];
    str = ['MAPE (MAP) = ' num2str(data.mape_mbp(rel_alg),'%5.1f'), ' %'];
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'FontSize', ftsize-2, 'Color', 'r');
    
    data.mape_co(rel_alg,1) = mean(data.ape( (data.group6 | data.group5), rel_alg));
    dim = [.5 .02 .2 .2];
    str = ['MAPE (CO) = ' num2str(data.mape_co(rel_alg),'%5.1f'), ' %'];
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'FontSize', ftsize-2, 'Color', 'b');
    
    % save plot
    PrintFigs(gcf, paper_size/70, [PATHS.CaseStudies, 'CO_analysis_corr_', curr_alg])
    
end

end

function PrintFigs(h, paper_size, savepath, close_plot)
set(h,'PaperUnits','inches');
set(h,'PaperSize', [paper_size(1), paper_size(2)]);
set(h,'PaperPosition',[0 0 paper_size(1) paper_size(2)]);
set(gcf,'color','w');
print(h,'-dpdf',savepath)
print(h,'-depsc',savepath)
%print(h,'-dpng',savepath)

% if you want .eps illustrations, then do as follows:
up.eps_figs = 0;
if up.eps_figs
    % you need to download 'export_fig' from:
    % http://uk.mathworks.com/matlabcentral/fileexchange/23629-export-fig
    export_fig_dir_path = 'C:\Documents\Google Drive\Work\Projects\PhD\Github\phd\Tools\Other Scripts\export_fig\altmany-export_fig-76bd7fa\';
    addpath(export_fig_dir_path)
    export_fig(savepath, '-eps')
end

if nargin > 3 && ~close_plot
else
    close all;
end

% save 
fid = fopen([savepath, '.txt'], 'w');
p = mfilename('fullpath');
p = strrep(p, '\', '\\');
fprintf(fid, ['Figures generated by:\n\n ' p '.m \n\n on ' datestr(today)]);
fclose all;

end