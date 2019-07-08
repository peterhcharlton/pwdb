function analyse_pwdb(pwdb_no)
% ANALYSE_PWDB analyses the Pulse Wave Database, producing results figures
% and tables.
%
%               analyse_pwdb
%
%   Inputs:     - 'pwdb_data.mat', which is produced by 'export_pwdb.m'.
%               - [optional] the '...path.mat' files produced by the
%                  same script
%
%   Outputs:    - Figures and Tables containing the results of the analyses.
%           
%   Accompanying Article:
%       This code is provided to facilitate reproduction of the Pulse Wave
%       Database described in: 
%           Charlton P.H. et al. Modelling arterial pulse waves in healthy
%           ageing: in silico evaluation of haemodynamics and
%           pulse wave indices, ~~ under review ~~  
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
% Contributed to by: Marie Willemet, Jordi Alastruey, and Peter H. Charlton
% v.1.0

fprintf('\n --- Analysing PWDB ---')

%% Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTINGS TO CHANGE: This function specifies where to save the outputs   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PATHS = setup_paths_for_post_processing(pwdb_no);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create folders in which to store results
create_folders(PATHS)

% Make table of haemodynamic characteristics
make_haem_characteristics_table(PATHS);

% Make demo PPG analysis plot
make_demo_ppg_analysis_plot(PATHS);

% Compare simulated characteristics with prescribed characteristics
make_simulated_vs_prescribed_figure(PATHS);

% Make wave speed figure
make_wave_speed_figure(PATHS);

% Compare simulated characteristics with those from literature
make_simulated_vs_literature_figure(PATHS);

% Make figure of baseline waves of 25-year old at common measurement sites
make_baseline_waves_figure(PATHS, pwdb_no);

% Make pressure-area plots to investigate hysteresis
make_pressure_area_plots(PATHS);

% Make plot of how gamma changes with arterial diameter
make_gamma_plots(PATHS);

% Make plots of how 25-year old waves vary with changes in initial params
make_initial_parameters_waves_figure(PATHS);

% Make figures of changes in wave shape with parameters
make_parameters_waves_figure(PATHS);

% Make figure
make_wrist_ppg_parameters_figure(PATHS);

% Make figure of PPG estimation
make_ppg_estimation_figure(PATHS);

% Make figure of changes in wave shape with age (comparing to literature)
make_changes_in_waves_age_figure(PATHS);

% Compare simulated characteristics when parameters are varied individually
make_parameters_characteristics_figure(PATHS);

% Perform sensitivity analyses
make_sensitivity_analysis_figures(PATHS);

% Make comparison of waves with Mynard's
comparison_w_mynard_waves_figure(PATHS);

% Make figure showing pulse waves along particular paths
make_pw_along_path_figure(PATHS);

fprintf('\n --- Finished analysing PWDB ---')

end

function create_folders(PATHS)

folders_to_make = {PATHS.Analysis_tables, PATHS.Analysis_figures};
for s = 1 : length(folders_to_make)
    if ~exist(folders_to_make{s}, 'dir')
        mkdir(folders_to_make{s})
    end
end

end

function make_demo_ppg_analysis_plot(PATHS)

fprintf('\n - Making BP and PPG PW analysis demo plots')

% load collated data
load(PATHS.exported_data_mat_pwdb_data)

% Setup figure for Pressure and Flow vel
fig_settings.sig_types = {'P', 'PPG'};
fig_settings.req_sites = {'Radial'};
baseline_sim_no = find(data.config.baseline_sim_for_all);

% cycle through different signals
for sig_type_no = 1 : length(fig_settings.sig_types)
    curr_sig_type = fig_settings.sig_types{sig_type_no};
    
    % cycle through different sites
    for req_site_no = 1 : length(fig_settings.req_sites)
        
        % extract relevant signal at this site
        curr_site = fig_settings.req_sites{req_site_no};
        eval(['rel_sig.v = data.waves.' curr_sig_type '_' curr_site '{baseline_sim_no};'])
        rel_sig.fs = data.waves.fs;
        rel_sig.t = [0:length(rel_sig.v)-1]/rel_sig.fs;
        
        % convert to friendly units if needed
        if sum(strcmp(curr_sig_type, {'P', 'Pe'}))
            rel_sig.v = rel_sig.v/133.33;   % Pa to mmHg
            rel_sig.units = 'mmHg';
        elseif strcmp(curr_sig_type, 'A')
            rel_sig.v = rel_sig.v*1000*1000;   % m^2 to mm^2
            rel_sig.units = 'mm^2';
        elseif strcmp(curr_sig_type, 'PPG')
            rel_sig.units = 'au';
        end
        
        % Make plot
        options.save_folder = PATHS.Analysis_figures;
        options.do_demo_plot = false;
        options.do_plot = 1; options.plot_third_deriv = 0;
        options.save_file = ['demo_' curr_sig_type '_' curr_site '_'];
        PulseAnalyse6(rel_sig, options);
        close all
        
    end
    
end

end

function make_haem_characteristics_table(PATHS)

fprintf('\n - Making table of haemodynamic characteristics')

%% load parameters
load(PATHS.exported_data_mat_pwdb_data)

%% Skip if this doesn't have all the required ages
if ~isequal(unique(data.pw_inds.age(:)'), 25:10:75)
    return
end

%% extract haemodynamic parameters

% select which parameters to extract
cardiac_params = {'HR', 'SV', 'CO', 'LVET', 'dPdt', 'PFT', 'RFV'};
arterial_params = {'SBP_a', 'DBP_a', 'MBP_a', 'PP_a', 'SBP_b', 'DBP_b', 'MBP_b', 'PP_b', 'SBP_f', 'DBP_f', 'MBP_f', 'PP_f', 'PP_amp', 'AP_a', 'AI_a', 'Tr_a', 'PWV_a', 'PWV_cf', 'PWV_br', 'PWV_fa', 'dia_asc_a', 'dia_desc_thor_a', 'dia_abd_a', 'len_prox_a', 'MBP_drop_finger', 'MBP_drop_ankle', 'AGI_mod', 'RI', 'SI', 'IAD'};
vascular_params = {'svr'};
if sum(strcmp(fieldnames(data.haemods), 'pvc'))
    vascular_params = [vascular_params, {'pvc','tau'}];
end
params = { 'age', cardiac_params{:}, arterial_params{:}, vascular_params{:}};

% extract each of the parameters in turn
for param_no = 1 : length(params)
    eval(['table_data.' params{param_no} ' = extractfield(data.haemods, ''' params{param_no} ''');'])
end
if sum(strcmp(fieldnames(data.haemods), 'pvc'))
    % change units
    table_data.pvc = 10e9*table_data.pvc;
end


% rename some of the variables
vars.old = {'AI_a', 'AP_a', 'Tr_a', 'pvc'};
vars.new = {'AIx', 'AP', 'Tr', 'pvc_adj'};
for var_no = 1 : length(vars.old)
    if strcmp(fieldnames(table_data), vars.old{var_no})
        eval(['table_data.' vars.new{var_no} ' = table_data.' vars.old{var_no} ';']);
        eval(['table_data = rmfield(table_data, ''' vars.old{var_no} ''');']);
        if sum(strcmp(arterial_params, vars.old{var_no}))
            arterial_params{strcmp(arterial_params, vars.old{var_no})} = vars.new{var_no};
        else
            vascular_params{strcmp(vascular_params, vars.old{var_no})} = vars.new{var_no};
        end
    end
end
clear vars var_no

%% Find out how many physiologically plausible subjects in each age group
ages = unique(data.config.age);
no_subjs.all = sum(data.plausibility.plausibility_log);
for age_no = 1 : length(ages)
    curr_age = ages(age_no);
    temp_no_subjs = sum(data.plausibility.plausibility_log(data.config.age == curr_age));
    eval(['no_subjs.age' num2str(curr_age) ' = temp_no_subjs;']);
    clear temp_no_subjs curr_age
end
clear age_no

%% Create Table

for s = 1 :2
    
    if s == 1, do_range = false; else, do_range = true; end
    csv_text = ['Haemodynamic Characteristic, All Subjects, 25, 35, 45, 55, 65, 75',newline];
    csv_text = [csv_text, 'n, ', num2str(no_subjs.all), ', ', num2str(no_subjs.age25), ', ', num2str(no_subjs.age35), ', ', num2str(no_subjs.age45), ', ', num2str(no_subjs.age55), ', ', num2str(no_subjs.age65), ', ', num2str(no_subjs.age75), newline];
    table_text = ['\\textbf{Cardiac} & & & & & & & \\\\', newline];
    csv_text = [csv_text, 'Cardiac , , , , , , ,', newline];
    [table_text, csv_text] = add_params_table_text(table_text, csv_text, cardiac_params, table_data, data.plausibility.plausibility_log, do_range);
    table_text = [table_text, '\\hline', newline];
    table_text = [table_text, '\\textbf{Arterial} & & & & & & & \\\\', newline];
    csv_text = [csv_text, 'Arterial , , , , , , ,', newline];
    [table_text, csv_text] = add_params_table_text(table_text, csv_text, arterial_params, table_data, data.plausibility.plausibility_log, do_range);
    table_text = [table_text, '\\hline', newline];
    table_text = [table_text, '\\textbf{Vascular Beds} & & & & & & & \\\\', newline];
    csv_text = [csv_text, 'Vascular Beds , , , , , , ,', newline];
    [table_text, csv_text] = add_params_table_text(table_text, csv_text, vascular_params, table_data, data.plausibility.plausibility_log, do_range);
    
    if ~do_range
        fid = fopen(PATHS.characteristics_table, 'w');
        fid2 = fopen([PATHS.characteristics_table,'.csv'], 'w');
    else
        fid = fopen(PATHS.characteristics_table_range, 'w');
        fid2 = fopen([PATHS.characteristics_table_range, '.csv'], 'w');
    end
    fprintf(fid, table_text);
    fprintf(fid2, csv_text);
    fclose(fid);
    fclose(fid2);
end

end

function [table_text, csv_text] = add_params_table_text(table_text, csv_text, curr_params, table_data, plausibility_log, do_range)

if nargin<4
    do_range = false;
end

for param_no = 1 : length(curr_params)
    
    curr_param = curr_params{param_no};
    
    if strcmp(curr_param(1:2), 'EX')
        continue
    end
    
    % Make label
    label = make_param_label(curr_param);
    table_line = ['- ', label];
    
    % Calculate stats
    eval(['rel_vals.v = table_data.' curr_param ';'])
    rel_vals.age = table_data.age;
    
    % - all subjects
    rel_els = ~isnan(rel_vals.v) & plausibility_log(:)';
    mean_val = mean(rel_vals.v(rel_els));
    std_val = std(rel_vals.v(rel_els));
    range_val = range(rel_vals.v(rel_els));
    if strcmp(curr_param, 'dPdt')
        mean_val = round(mean_val);
        std_val = round(std_val);
    end
    if sum(strcmp(curr_param, {'PP_amp', 'CO', 'tau', 'RI', 'SI', 'AGI_mod'}))
        if ~do_range
            table_line = [table_line, ' & ' num2str(mean_val, '%.2f') ' $\\pm$ '  num2str(std_val, '%.2f') ];
        else
            table_line = [table_line, ' & ' num2str(mean_val, '%.2f') ' $($'  num2str(range_val, '%.2f') '$)$' ];
        end
    else
        if ~do_range
            table_line = [table_line, ' & ' num2str(mean_val, '%.1f') ' $\\pm$ '  num2str(std_val, '%.1f') ];
        else
            table_line = [table_line, ' & ' num2str(mean_val, '%.1f') ' $($'  num2str(range_val, '%.1f') '$)$' ];
        end
    end
    
    
    % - individual ages
    ages = unique(rel_vals.age);
    for age_no = 1 : length(ages)
        curr_age = ages(age_no);
        rel_els = find(rel_vals.age == curr_age & ~isnan(rel_vals.v) & plausibility_log(:)');
        mean_val = mean(rel_vals.v(rel_els));
        std_val = std(rel_vals.v(rel_els));
        if strcmp(curr_param, 'dPdt')
            mean_val = round(mean_val);
            std_val = round(std_val);
        end
        if length(rel_els) ~= 2
            range_val = range(rel_vals.v(rel_els));
        else
            range_val = rel_vals.v(rel_els(2)) - rel_vals.v(rel_els(1));
        end
        if sum(strcmp(curr_param, {'PP_amp', 'CO', 'tau', 'RI', 'SI', 'AGI_mod'}))
            if ~do_range
                table_line = [table_line, ' & ' num2str(mean_val, '%.2f') ' $\\pm$ '  num2str(std_val, '%.2f') ];
            else
                table_line = [table_line, ' & ' num2str(mean_val, '%.2f') ' $($'  num2str(range_val, '%.2f') '$)$' ];
            end
        else
            if ~do_range
                table_line = [table_line, ' & ' num2str(mean_val, '%.1f') ' $\\pm$ '  num2str(std_val, '%.1f') ];
            else
                table_line = [table_line, ' & ' num2str(mean_val, '%.1f') ' $($'  num2str(range_val, '%.1f') '$)$' ];
            end
        end
    end
    table_line = [table_line, '\\\\', newline];
    
    table_text = [table_text, table_line];
    
    % create csv line
    csv_line = strrep(table_line, ' &', ',');
    csv_line = strrep(csv_line, '$\\pm$', '+/-');
    csv_line = strrep(csv_line, ', systolic', '');
    csv_line = strrep(csv_line, ', diastolic', '');
    csv_line = strrep(csv_line, ', mean', '');
    csv_line = strrep(csv_line, ', pulse pressure', '');
    csv_line = strrep(csv_line, '\\\\', '');
    %csv_line = strrep(csv_line, '\\%', '%%');
    csv_line = strrep(csv_line, '\\hspace{4.1cm}', '                                            ');
    csv_line = strrep(csv_line, '\\hspace{4.5cm}', '                                                ');
    csv_line = strrep(csv_line, '\\hspace{4.43cm}', '                                                ');
    csv_line = strrep(csv_line, '\\hspace{2.7cm}', '                                ');
    csv_text = [csv_text, csv_line];
    clear table_line csv_line
    
end

end

function make_simulated_vs_prescribed_figure(PATHS)

fprintf('\n - Making simulated vs prescribed figure')

only_baseline_sims = true;

% Load data
load(PATHS.exported_data_mat_pwdb_data)
ages = data.config.age;

% Prescribed characteristics
cardiac_params = {'HR', 'SV', 'LVET', 'PFT', 'RFV'};
arterial_params = {'MBP', 'PWV_cf', 'PWV_br', 'PWV_fa', 'dia_asc_a', 'dia_desc_thor_a', 'dia_abd_a', 'dia_car', 'len_prox_a'};
params = {cardiac_params{:}, arterial_params{:}};
for param_no = 1 : length(params)
    eval(['prescribed.' params{param_no} ' = nan(length(data.config.age),1);'])
    eval(['simulated.' params{param_no} ' = nan(length(data.config.age),1);'])
end

%% Extract prescribed and simulated values of each characteristic
for sim_no = 1 : length(data.config.age)
    for param_no = 1 : length(params)
        
        curr_param = params{param_no};
        
        % Extract simulated value of this parameter
        if sum(strcmp(fieldnames(data.haemods), params{param_no}))
            % if this parameter has already been extracted in the table_data
            eval(['simulated.' curr_param '(sim_no) = data.haemods(sim_no).' curr_param ';'])
        end
        
        % Extract prescribed value of this parameter
        if sum(strcmp(lower(fieldnames(data.config)), lower(params{param_no})))
            eval(['prescribed.' curr_param '(sim_no) = data.config.' lower(curr_param) '(sim_no);'])
            rel_cols = find(~strcmp(lower(data.config.variations.param_names), lower(params{param_no})));
            eval(['prescribed.' curr_param '_rel(sim_no) = sum(abs(data.config.variations.params(sim_no,rel_cols)))==0;'])
            %eval(['prescribed.' curr_param '_base(sim_no) = sum(abs(data.config.variations.params(sim_no)))==0;'])
        else
            eval(['prescribed.' curr_param '(sim_no) = data.config.desired_chars.' lower(curr_param) '(sim_no);'])
            rel_param_name = strrep(curr_param, 'PWV_cf', 'pwv');
            rel_param_name = strrep(rel_param_name, 'PWV_br', 'pwv');
            rel_param_name = strrep(rel_param_name, 'PWV_fa', 'pwv');
            rel_param_name = strrep(rel_param_name, 'dia_asc_a', 'dia');
            rel_param_name = strrep(rel_param_name, 'dia_desc_thor_a', 'dia');
            rel_param_name = strrep(rel_param_name, 'dia_abd_a', 'dia');
            rel_param_name = strrep(rel_param_name, 'dia_car', 'dia');
            rel_param_name = strrep(rel_param_name, 'len_prox_a', 'len');
            rel_cols = find(~strcmp(data.config.variations.param_names, rel_param_name));
            eval(['prescribed.' curr_param '_rel(sim_no) = sum(abs(data.config.variations.params(sim_no,rel_cols)))==0;'])
            %eval(['prescribed.' curr_param '_base(sim_no) = sum(abs(data.config.variations.params(sim_no,:)))==0;'])
        end
        
    end
end

%% Make plots of simulated characteristics

ftsize = 20; lwidth = 2;
paper_size = [500,350];

for param_no = 1 : length(params)
        
    % extract data for this param
    curr_param = params{param_no};
    [~, units] = make_param_label(curr_param);
    eval(['rel_data.p = prescribed.' curr_param ';']);
    eval(['rel_data.p_rel = prescribed.' curr_param '_rel;']);
    eval(['rel_data.s = simulated.' curr_param ';']);
    
    % find mean and std of parameter at each age
    unique_ages = unique(data.config.age);
    for age_no = 1 : length(unique_ages)
        if only_baseline_sims
            rel_els = ages == unique_ages(age_no) & rel_data.p_rel(:);
        else
            rel_els = ages == unique_ages(age_no);
        end
        rel_data.p_mean(age_no) = mean(rel_data.p(rel_els));
        rel_data.p_std(age_no) = std(rel_data.p(rel_els));
        rel_data.s_mean(age_no) = mean(rel_data.s(rel_els));
        rel_data.s_std(age_no) = std(rel_data.s(rel_els));
    end
    
    % make prescribed change with age plot
    plot(unique_ages, rel_data.p_mean, 'o-k', 'LineWidth', lwidth, 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'), hold on
    plot(unique_ages, rel_data.p_mean - rel_data.p_std, '--k', 'LineWidth', lwidth)
    plot(unique_ages, rel_data.p_mean + rel_data.p_std, '--k', 'LineWidth', lwidth)
    
    temp.min = min(rel_data.p(:));
    temp.max = max(rel_data.p(:));
    temp.range = temp.max-temp.min;
    lims = [floor(temp.min-0.1*temp.range), ceil(temp.max+0.1*temp.range)]; clear temp
    if length(unique(lims)) == 1
        lims = [lims(1)-1 lims(2)+1];
    end
    
    % tidy-up
    xlim([min(unique_ages)-5, max(unique_ages)+5])
    ylim(lims)
    xlabel('Age [yrs]', 'FontSize', ftsize)
    set(gca, 'FontSize', ftsize)
    title([strrep(curr_param,'_', ' '), ' [', units,']'], 'FontSize', ftsize)
    
    % save
    PrintFigs(gcf, paper_size/70, [PATHS.Analysis_figures, 'pres_', curr_param])
        
end

%% Make plots of prescribed vs simulated characteristics

ftsize = 20; lwidth = 2;
paper_size = [1200,350];
only_baseline_sims = false;

for param_no = 1 : length(params)
        
    % extract data for this param
    curr_param = params{param_no};
    [~, units] = make_param_label(curr_param);
    eval(['rel_data.p = prescribed.' curr_param ';']);
    eval(['rel_data.p_rel = prescribed.' curr_param '_rel;']);
    eval(['rel_data.s = simulated.' curr_param ';']);
    
    % make scatter plot
    figure('Position', [20,20,paper_size])
    subplot(1,3,1)
    plot([0,10000],[0,10000], '--k'), hold on % line of identity
    plot(rel_data.p, rel_data.s, 'xk', 'LineWidth', lwidth, 'MarkerSize', 13)
    
    % tidy-up
    temp.min = min([rel_data.p(:); rel_data.s(:)]);
    temp.max = max([rel_data.p(:); rel_data.s(:)]);
    temp.range = temp.max-temp.min;
    lims = [floor(temp.min-0.1*temp.range), ceil(temp.max+0.1*temp.range)]; clear temp
    if length(unique(lims)) == 1
        lims = [lims(1)-1 lims(2)+1];
    end
    xlim(lims)
    ylim(lims)
    xlabel('Prescribed', 'Color', 'b', 'FontSize', ftsize)
    ylabel('Simulated', 'Color', 'r', 'FontSize', ftsize)
    set(gca, 'FontSize', ftsize)
    title([strrep(curr_param,'_', ' '), ' [', units,']'], 'FontSize', ftsize)
    
    % find mean and std of parameter at each age
    unique_ages = unique(ages);
    for age_no = 1 : length(unique_ages)
        if only_baseline_sims
            rel_els = ages == unique_ages(age_no) & rel_data.p_rel(:);
        else
            rel_els = ages == unique_ages(age_no);
        end
        rel_data.p_mean(age_no) = mean(rel_data.p(rel_els));
        rel_data.p_std(age_no) = std(rel_data.p(rel_els));
        rel_data.s_mean(age_no) = mean(rel_data.s(rel_els));
        rel_data.s_std(age_no) = std(rel_data.s(rel_els));
    end
    
    % make prescribed change with age plot
    subplot(1,3,2)
    plot(unique_ages, rel_data.p_mean, 'o-b', 'LineWidth', lwidth, 'MarkerSize', 7, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'), hold on
    plot(unique_ages, rel_data.p_mean - rel_data.p_std, '--b', 'LineWidth', lwidth)
    plot(unique_ages, rel_data.p_mean + rel_data.p_std, '--b', 'LineWidth', lwidth)
    
    % tidy-up
    xlim([min(unique_ages)-5, max(unique_ages)+5])
    ylim(lims)
    xlabel('Age [yrs]', 'FontSize', ftsize)
    ylabel([strrep(curr_param,'_', ' '), ' [', units,']'], 'FontSize', ftsize)
    set(gca, 'FontSize', ftsize)
    title('Prescribed', 'Color', 'b', 'FontSize', ftsize)
    
    % make simulated change with age plot
    subplot(1,3,3)
    plot(unique_ages, rel_data.s_mean, 'o-r', 'LineWidth', lwidth, 'LineWidth', lwidth, 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'), hold on
    plot(unique_ages, rel_data.s_mean - rel_data.s_std, '--r', 'LineWidth', lwidth)
    plot(unique_ages, rel_data.s_mean + rel_data.s_std, '--r', 'LineWidth', lwidth)
    
    % tidy-up
    xlim([min(unique_ages)-5, max(unique_ages)+5])
    ylim(lims)
    xlabel('Age [yrs]', 'FontSize', ftsize)
    ylabel([strrep(curr_param,'_', ' '), ' [', units,']'], 'FontSize', ftsize)
    set(gca, 'FontSize', ftsize)
    title('Simulated', 'Color', 'r', 'FontSize', ftsize)
    
    % save
    PrintFigs(gcf, paper_size/70, [PATHS.Analysis_figures, 'sim_vs_pres_', curr_param])
        
end


end

function make_baseline_waves_figure(PATHS, sim_no)

% load data
load(PATHS.exported_data_mat_pwdb_data)

% Setup figure for Pressure and Flow vel
paper_size = [350, 600];
figure('Position', [20,20,paper_size])
fig_settings.lwidth = 2;
fig_settings.ftsize = 12;
fig_settings.req_sites =  {'Carotid', 'AorticRoot', 'Brachial', 'Radial', 'Digital', 'Femoral', 'AntTibial'};
fig_settings.sig_types = {'P', 'U'};
fig_settings.colors = {'r', 'r'};
fig_settings.ylims = {[60,125], [-0.15, 0.75]};

% Plot
plot_baseline_signals(fig_settings, data)

% Save figure
PrintFigs(gcf, paper_size/70, PATHS.baseline_waves_fig_a)

% Setup figure for PPG and Area
figure('Position', [20,20,paper_size])
fig_settings.req_sites =  {'SupTemporal', 'Carotid', 'Brachial', 'Radial', 'Digital', 'Femoral', 'AntTibial'};
fig_settings.sig_types = {'A', 'PPG'};
fig_settings.colors = {'b', 'b'};
fig_settings.ylims = {'auto', [0, 1.1]};

% Plot
plot_baseline_signals(fig_settings, data)

% Save figure
PrintFigs(gcf, paper_size/70, PATHS.baseline_waves_fig_b)

end

function make_simulated_vs_literature_figure(PATHS)

fprintf('\n - Making simulated vs literature figure')

only_baseline_sims = true;

% Load data
load(PATHS.exported_data_mat_pwdb_data)
ages = data.config.age;

% Simulated characteristics
cardiac_params = {'CO'};
arterial_params = {'DBP_a', 'AP_a', 'SBP_b', 'DBP_b', 'PP_b', 'MBP_b', 'SBP_a', 'PP_a', 'PP_amp', 'AI_a', 'Tr_a'};
params = {cardiac_params{:}, arterial_params{:}};

%% Extract literature values of each characteristic
for param_no = 1 : length(params)
    
    curr_param = params{param_no};
    
    % Extract literature values of this parameter
    switch curr_param
        case 'CO' % Taken from Le2016 (male data, Table 2)
            literature_data.age = 25:10:65;
            literature_data.mean = [6.6, 6.2, 5.8, 5.4, 5.0];
            literature_data.sd = [8.5-4.65, 8.1-4.25, 7.7-3.85, 7.3-3.45, 6.9-3.05]./(2*1.96);
        case 'SBP_b' % taken from McEniery2005 (male data, Table 1)
            literature_data.age = 15:10:85;
            literature_data.mean = [123, 124, 123, 125, 125, 126, 127, 130];
            literature_data.sd = [10,10,9,9,9,9,9,8];
        case 'DBP_a' % taken from McEniery2005 (male data, Table 1)
            %%%%%%%%%%%%%% THIS IS BRACHIAL RATHER THAN AORTIC DATA %%%%%%%
            literature_data.age = 15:10:85;
            literature_data.mean = [73,75,77,79,79,78,76,75];
            literature_data.sd = [8,10,9,9,9,9,9,8];
        case 'DBP_b' % taken from McEniery2005 (male data, Table 1)
            literature_data.age = 15:10:85;
            literature_data.mean = [73,75,77,79,79,78,76,75];
            literature_data.sd = [8,10,9,9,9,9,9,8];
        case 'PP_b' % taken from McEniery2005 (male data, Table 1)
            literature_data.age = 15:10:85;
            literature_data.mean = [50,49,47,46,46,49,51,55];
            literature_data.sd = [9,9,8,7,8,8,8,9];
        case 'MBP_b' % taken from McEniery2005 (male data, Table 1)
            literature_data.age = 15:10:85;
            literature_data.mean = [88,89,92,95,95,94,93,92];
            literature_data.sd = [8,8,8,7,7,7,7,8];
        case 'SBP_a' % taken from McEniery2005 (male data, Table 1)
            literature_data.age = 15:10:85;
            literature_data.mean = [103,105,109,113,115,117,118,120];
            literature_data.sd = [8,8,9,9,9,9,9,8];
        case 'PP_a' % taken from McEniery2005 (male data, Table 1)
            literature_data.age = 15:10:85;
            literature_data.mean = [29,30,31,34,35,39,42,45];
            literature_data.sd = [5,6,6,6,7,7,7,9];
        case 'PP_amp' % taken from McEniery2005 (male data, Table 1)
            literature_data.age = 15:10:85;
            literature_data.mean = [1.72, 1.7 , 1.50, 1.39, 1.33, 1.26, 1.24, 1.25];
            literature_data.sd =   [0.11, 0.14, 0.18, 0.15, 0.16, 0.13, 0.12, 0.15];
        case 'AP_a' % taken from McEniery2005 (male data, Table 1)
            % This is aortic rather than carotid data
            literature_data.age = 15:10:85;
            literature_data.mean = [-1,1,4,7,9,11,13,14];
            literature_data.sd = [3,4,5,4,5,5,5,5];
        case 'AI_a' % taken from McEniery2005 (male data, Table 1)
            % This is aortic rather than carotid data
            literature_data.age = 15:10:85;
            literature_data.mean = [-2,2,12,19,24,28,30,30];
            literature_data.sd = [8,11,13,10,10,9,9,10];
        case 'Tr_a' % taken from McEniery2005 (male data, Table 1)
            literature_data.age = 15:10:85;
            literature_data.mean = [150,154,151,148,143,141,136,133];
            literature_data.sd = [17,21,21,16,15,12,12,16];
    end
    
    % store literature values
    eval(['literature.' curr_param ' = literature_data;']);
    
    % Extract simulated values
    simulated_data.age = unique(data.config.age(:)');
    for age_no = 1 : length(simulated_data.age)
        curr_age = simulated_data.age(age_no);
        rel_els = data.config.age == curr_age & data.plausibility.plausibility_log;
        param_data = extractfield(data.haemods, curr_param);
        rel_vals = param_data(rel_els); clear param_data rel_els curr_age
        %%%%%%%%%%%%% THIS IS IGNORING NANs %%%%%%%%%%%%%%%%
        simulated_data.mean(age_no) = mean(rel_vals(~isnan(rel_vals)));
        simulated_data.sd(age_no) = std(rel_vals(~isnan(rel_vals)));
        clear rel_vals
    end
    clear age_no
    
    % store simualted values
    eval(['simulated.' curr_param ' = simulated_data;']);
    
    clear literature_data simulated_data
end

%% Make plots of literature vs simulated characteristics
up.ylim_offset = 0.1;
req_color = [0,0,0];
age_ticks = 20:10:80;
ftsize = 22; lwidth = 2;

for param_no = 1 : length(params)
    
    curr_param = params{param_no};
    eval(['literature_data = literature.' curr_param ';'])
    eval(['simulated_data = simulated.' curr_param ';'])
        
    paper_size = [900,400];
    figure('Position', [100,100,paper_size])
    subplot('Position', [0.12,0.16,0.38,0.72])
    
    % Find ylims
    temp(1) = min([literature_data.mean-literature_data.sd, simulated_data.mean-simulated_data.sd]);
    temp(2) = max([literature_data.mean+literature_data.sd, simulated_data.mean+simulated_data.sd]);
    rel_ylims(1) = temp(1) - 0.1*(range(temp));
    rel_ylims(2) = temp(2) + 0.1*(range(temp));
    
    % - Literature
    
    % plot SD
    no_sd = 1;
    plot_sds(literature_data.age, literature_data.mean, literature_data.sd, no_sd, req_color, up);
    hold on
    
    % plot mean value
    plot(literature_data.age, literature_data.mean, 'LineWidth', lwidth, 'Color', req_color), hold on,
    
    % tidy up
    [~, unit, abbr, ~, graph_title_no_units] = make_param_label(curr_param);
    unit = strrep(unit, '%%', '%');
    abbr = strrep(abbr, 'PP_{amp}', 'abc');
    temp = strfind(abbr, '_');
    if ~isempty(temp), abbr = abbr(1:temp-1); end
    abbr = strrep(abbr, 'abc', 'PP_{amp}');
    xlabel('Age [years]', 'FontSize', ftsize)
    ylab = ylabel([abbr, ' [', unit, ']'], 'FontSize', ftsize);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.15, 0.5, 0]);
    set(gca, 'FontSize', ftsize, 'XTick', age_ticks)
    ylim(rel_ylims);
    xlim([20,80])
    grid on
    box off
    
    % - Simulated
    subplot('Position', [0.59,0.16,0.38,0.72])
    
    % plot SD
    no_sd = 1;
    plot_sds(simulated_data.age, simulated_data.mean, simulated_data.sd, no_sd, req_color, up);
    hold on
    
    % plot mean value
    plot(simulated_data.age, simulated_data.mean, 'o-', 'LineWidth', lwidth, 'Color', req_color), hold on,
    
    % tidy up
    xlabel('Age [years]', 'FontSize', ftsize)
    set(gca, 'FontSize', ftsize, 'XTick', age_ticks)
    ylim(rel_ylims);
    xlim([20,80])
    grid on
    box off
    
    % - title
    str = graph_title_no_units;
    dim = [0.2,0.7,0.6,0.3];
    annotation('textbox',dim,'String',str,'LineStyle', 'none', 'FontSize', ftsize+8, 'HorizontalAlignment', 'center');
    
    % save
    PrintFigs(gcf, paper_size/70, [PATHS.Analysis_figures, 'lit_vs_sim_', curr_param])
        
end

end

function make_parameters_characteristics_figure(PATHS)

fprintf('\n - Making characteristics vs parameters figures')


% load data
load(PATHS.exported_data_mat_pwdb_data)

% identify baseline simulation for baseline age
baseline_sim_no = find(data.config.baseline_sim_for_all);
baseline_age = data.config.age(baseline_sim_no);

% for each parameter, identify simulations where only that parameter changed without anything else changing.
all_params = data.config.variations.params;
param_names = data.config.variations.param_names;

% add age to list of parameters
param_names = [{'age'}; param_names];
all_params = [zeros(size(all_params,1),1), all_params];

% find relevant simulations for each parameter
[rel_sims, param_variations] = deal(cell(length(param_names),1));
for param_no = 1 : length(param_names)
    columns_for_other_params = setxor(1:length(param_names), param_no);
    temp = all_params(:,columns_for_other_params);
    if ~strcmp(param_names{param_no}, 'age')
        rel_sims{param_no} = find(~any(temp,2) & data.config.age == baseline_age);
    else
        rel_sims{param_no} = find(~any(temp,2));
    end
    param_variations{param_no} = all_params(rel_sims{param_no},param_no);
end
clear temp

%% Make plots of parameters vs characteristics

% Simulated characteristics
cardiac_chars = {'CO'};
arterial_chars = {'PWV_a', 'MBP_drop_finger', 'MBP_drop_ankle', 'SBP_b', 'DBP_b', 'PP_b', 'MBP_b', 'SBP_a', 'PP_a', 'AP_a', 'PP_amp', 'AI_a', 'Tr_a'};
resultant_chars = {cardiac_chars{:}, arterial_chars{:}};

% cycle through input parameters
for param_no = 1 : length(param_names)
    curr_param = param_names{param_no};    
    
    % rel sims
    curr_rel_sims = rel_sims{param_no};
    
    if length(curr_rel_sims) < 2
        continue
    end
    
    % extract data for this parameter
    curr_param_vals = nan(length(curr_rel_sims),1);
    eval(['curr_param_vals = data.config.' curr_param '(curr_rel_sims);']);
    true_val = true;
    
    % make plot of characteristic vs parameter
    
    % cycle through simulated characteristics
    paper_size = [200,200,700,300];
    for res_char_no = 1 : length(resultant_chars)
        curr_res_char = resultant_chars{res_char_no};
        
        % extract data for this characteristic
        eval(['all_char_values = extractfield(data.haemods, ''' curr_res_char ''');']);
        char_values = all_char_values(curr_rel_sims);
        
        % sort according to param values
        [new_curr_param_vals, order] = sort(curr_param_vals);
        new_char_values = char_values(order); clear order char_values
        
        % make plot of characteristic vs parameter
        up.ylim_offset = 0.1;
        req_color = [0,0,0];
        age_ticks = 20:10:80;
        ftsize = 22; lwidth = 2;
        paper_size = [550,400];
        figure('Position', [100,100,paper_size])
        subplot('Position', [0.2,0.2,0.77,0.76])
        
        plot(new_curr_param_vals, new_char_values, 'o-', 'MarkerSize', 7, 'MarkerFaceColor', req_color, 'LineWidth', lwidth, 'Color', req_color)
        if range(all_char_values) == 0
            temp = 0.05*unique(all_char_values);
        else
            temp = 0.1*range(all_char_values);
        end
        ylim([min(all_char_values)-temp, max(all_char_values)+temp])
        
        % tidy up
        [~, unit, abbr, label] = make_param_label(curr_param);
        if true_val
            x_label_text = label;
        else
            temp = strfind(label, '[');
            x_label_text = [label(1:temp-2), '  [no SDs]'];
        end
        xlabel(x_label_text, 'FontSize', ftsize)
        [~, unit, abbr, label] = make_param_label(curr_res_char);
        ylab = ylabel(label, 'FontSize', ftsize);
        set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
        set(gca, 'FontSize', ftsize)
        grid on
        
        % save
        PrintFigs(gcf, paper_size/70, [PATHS.Analysis_figures, 'param_', curr_param, '_vs_char_', curr_res_char])
        
        clear new*
    end
    
    clear curr_rel_sims res_char_no curr_res_param char_values true_val curr_param_vals rel_col rel_sim_no
end

end

function make_ppg_estimation_figure(PATHS)

fprintf('\n - Making PPG estimation figure')

% load data
load(PATHS.exported_data_mat_pwdb_data)

% identify baseline simulation for baseline age
baseline_sim_no = find(data.config.baseline_sim_for_all);
baseline_age = data.config.age(data.config.baseline_sim_for_all);


%% Make plots of waves vs parameters

% Setup plotting
rel_sites = {'Carotid', 'Radial', 'Digital'};
ftsize = 24;
lwidth = 2;
offset = 0;
paper_size = [550,350];

% make plots for each wave measurement site
for site_no = 1 : length(rel_sites)
    curr_site = rel_sites{site_no};
    artery_name = find_artery_name(curr_site); % use PPG names
    
    % Make figure
    figure('Position', [20,20, paper_size])
    subplot('Position', [0.08,0.19,0.90,0.75])
    
    leg_labels = {}; max_t = 0;
    units = 'normalised';
    
    % Plot each of the waves in turn (corresponding to the three ages)
    colors = [0,0,0];
    
    % extract data to plot
    eval(['curr.ppg = data.waves.PPG_' curr_site '{baseline_sim_no};']);
    eval(['curr.p = data.waves.P_' curr_site '{baseline_sim_no};']);
    curr.fs = data.waves.fs;
    curr.t = [0:length(curr.p)-1]/curr.fs;
    
    % plot
    curr.p = (curr.p - min(curr.p))/range(curr.p);
    curr.ppg = (curr.ppg - min(curr.ppg))/range(curr.ppg);
    plot(curr.t, curr.p, '--', 'Color', colors(1,:), 'LineWidth', lwidth), hold on
    plot(curr.t, curr.ppg +offset, 'Color', colors(1,:), 'LineWidth', lwidth), hold on
    ylim([-0.1 1.1+offset]);
    set(gca, 'YTick', [])
    
    % store legend label
    leg_labels = {'Pressure', 'PPG'};
    max_t = max(curr.t);
    
    % tidy up
    xlim([0, max_t])
    xlabel('Time [s]', 'FontSize', ftsize)
    ylab = ylabel('Pulse Waves [norm.]', 'FontSize', ftsize);
    %ylab = ylabel({curr_wave_type, ['[' units ']']}, 'FontSize', ftsize, 'Rotation', 0);
    %set(ylab, 'Units', 'Normalized', 'Position', [-0.17, 0.5, 0]);
    set(gca, 'FontSize', ftsize)
    legend(leg_labels), clear leg_labels
    dim = [0.15,0.7,0.7,0.3];
    annotation('textbox',dim,'String',artery_name,'LineStyle', 'none', 'FontSize', ftsize, 'HorizontalAlignment', 'center');
    box off
    
    % save
    PrintFigs(gcf, paper_size/70, [PATHS.Analysis_figures, 'estimate_PPG_at_' strrep(artery_name, ' ', '_')])
    
end
clear curr_rel_sims res_char_no curr_res_param char_values true_val curr_param_vals rel_col rel_sim_no

end

function make_gamma_plots(PATHS)

fprintf('\n - Making gamma plots')

ftsize = 24;
lwidth = 2;
offset = 0.25;
paper_size = [600,400];

% load data
load(PATHS.exported_data_mat_pwdb_data)

rel_sims = find(data.config.baseline_sim_for_all);

% Make figure
figure('Position', [20,20, paper_size])

for sim_no = 1 : length(rel_sims)
    
    curr_sim_no = rel_sims(sim_no);
    
    % extract input data
    in = data.config.constants;
    r = mean([data.config.network.inlet_radius(curr_sim_no,:); data.config.network.outlet_radius(curr_sim_no,:)]);
    r = sort(r);
    a = pi*r.^2;
    gamma = 0.001*((in.gamma_b1(curr_sim_no)./(100*2*r))+in.gamma_b0(curr_sim_no))./a;
    
    % plot
    plot(2*r/100,gamma, 'LineWidth', lwidth), hold on,
    
    % tidy up
    xlabel('Diameter [cm]', 'FontSize', ftsize)
    ylab = ylabel('Gamma', 'FontSize', ftsize);
    set(gca, 'FontSize', ftsize)
    
end

% save
PrintFigs(gcf, paper_size/70, [PATHS.Analysis_figures, 'gamma'])
    
end

function make_pressure_area_plots(PATHS)

fprintf('\n - Making pressure-area plots')

%% Mynard data

% Import waves from Mynard's article
loadpath = '/Users/petercharlton/Google Drive/Work/Code/nektar/Mynard2015_data/ABME2015data.mat';
load(loadpath);

% Setup plotting
rel_domains = [1,15,21,22,42, 46, 49, 87];
rel_sites = {'AorticRoot', 'Carotid', 'Brachial', 'Radial', 'CommonIliac', 'Femoral', 'AntTibial', 'SupTemporal', 'Digital'};
ftsize = 24;
lwidth = 2;
offset = 0.25;
paper_size = [600,400];
% 
% % make plots for each measurement site
%     
% for domain_no = 1 : length(rel_domains)
%     curr_domain_no = rel_domains(domain_no);
%     artery_name = find_artery_name(curr_domain_no);
%     
%     % Make figure
%     figure('Position', [20,20, paper_size])
%     
%     % extract relevant wave from Mynard's data
%     lit_wave.t = a115lb.t;
%     lit_wave.fs = 1000;
%     switch curr_domain_no
%         case 1
%             mynard_name = 'AoRt';
%         case 15
%             mynard_name = 'Lcar';
%         case 21
%             mynard_name = 'LBrach';
%         case 22
%             mynard_name = 'LRadI';
%         case 42
%             mynard_name = 'LInIl';
%         case 46
%             mynard_name = 'Lfem';
%         case 49
%             mynard_name = 'LATib';
%         case 87
%             mynard_name = 'LSupTemp';
%         case 112
%             mynard_name = 'AoRt';
%     end
%     wave_el = find(strcmp(a115lb.monitor.name, mynard_name));
%     lit_wave.p = a115lb.tnode.p(:,wave_el)/1333.3;
%     lit_wave.a = a115lb.tnode.A(:,wave_el);
%     
%     % Find best fit line
%     f=fit(lit_wave.p,lit_wave.a,'poly1');
%     fit_line = feval(f,lit_wave.p);
%     
%     % plot
%     plot(lit_wave.p,fit_line, '--r', 'LineWidth', lwidth), hold on,
%     plot(lit_wave.p, lit_wave.a, 'b', 'LineWidth', lwidth), hold on
%         
%     % tidy up
%     xlabel('Pressure [mmHg]', 'FontSize', ftsize)
%     ylab = ylabel('Area [cm^2]', 'FontSize', ftsize);
%     set(gca, 'FontSize', ftsize)
%     title(artery_name, 'FontSize', ftsize+4)
%     
%     % legend
%     legend({'Best fit line', 'Experimental'}, 'Location', 'SouthEast')
%     
%     % save
%     PrintFigs(gcf, paper_size/70, [PATHS.ProcessedData, 'Mynard_p_a_', strrep(artery_name, ' ', '_')])
%     
% end
% clear a115lb



%% Our data

% load data
load(PATHS.exported_data_mat_pwdb_data)

% identify baseline simulation for baseline age
baseline_sim_no = find(data.config.baseline_sim_for_all);

% make plots for each measurement site
for site_no = 1 : length(rel_sites)
    curr_site = rel_sites{site_no};
    
    % Make figure
    figure('Position', [20,20, paper_size])
    
    % extract data to plot
    eval(['wave.p = data.waves.P_' curr_site '{baseline_sim_no};']);
    eval(['wave.a = data.waves.A_' curr_site '{baseline_sim_no}*(100^2);']);  % convert from m2 to cm2
    
    % Find best fit line
    f=fit(wave.p,wave.a,'poly1');
    fit_line = feval(f,wave.p);
    
    % plot
    plot(wave.p,fit_line, '--r', 'LineWidth', lwidth), hold on,
    plot(wave.p, wave.a, 'b', 'LineWidth', lwidth), hold on
        
    % tidy up
    xlabel('Pressure [mmHg]', 'FontSize', ftsize)
    ylab = ylabel('Area [cm^2]', 'FontSize', ftsize);
    set(gca, 'FontSize', ftsize)
    title(curr_site, 'FontSize', ftsize+4)
    
    % legend
    legend({'Best fit line', 'Experimental'}, 'Location', 'SouthEast')
    
    % save
    PrintFigs(gcf, paper_size/70, [PATHS.Analysis_figures, 'Simulated_p_a_', curr_site])
    
end

end

function make_parameters_waves_figure(PATHS)

fprintf('\n - Making waves vs parameters figures')

% load data
load(PATHS.exported_data_mat_pwdb_data)

% identify baseline simulation for baseline age
baseline_sim_no = find(data.config.baseline_sim_for_all);
baseline_age = data.config.age(baseline_sim_no);

% for each parameter, identify simulations where only that parameter changed without anything else changing.
all_params = data.config.variations.params;
param_names = data.config.variations.param_names;

% add age to list of parameters
param_names = [{'age'}; param_names];
all_params = [zeros(size(all_params,1),1), all_params];

% find relevant simulations for each parameter
[rel_sims, param_variations] = deal(cell(length(param_names),1));
for param_no = 1 : length(param_names)
    columns_for_other_params = setxor(1:length(param_names), param_no);
    temp = all_params(:,columns_for_other_params);
    if ~strcmp(param_names{param_no}, 'age')
        rel_sims{param_no} = find(~any(temp,2) & data.config.age == baseline_age);
    else
        rel_sims{param_no} = find(~any(temp,2));
    end
    param_variations{param_no} = all_params(rel_sims{param_no},param_no);
end
clear temp

%% Make plots of waves vs parameters

% Setup plotting
wave_types = {'Q', 'A', 'PPG', 'U', 'P'};
rel_sites.P = {'AorticRoot', 'Carotid', 'Brachial', 'Radial', 'Femoral', 'Digital', 'SupMidCerebral'};
rel_sites.PPG = {'Digital', 'Carotid', 'Brachial', 'Radial', 'SupTemporal'};
rel_sites.U = {'AorticRoot', 'Digital', 'Carotid', 'Brachial', 'Radial', 'Femoral'};
rel_sites.Q = {'AorticRoot', 'Carotid', 'SupMidCerebral'};
rel_sites.A = {'AorticRoot', 'Carotid'};
ftsize = 20;
lwidth = 2;
plot_model_data_normalised = false;
offset = 0.25;
paper_size = [700,300];
    
% cycle through input parameters
for param_no = 1 : length(param_names)
    curr_param = param_names{param_no};    
    
    % rel sims
    curr_rel_sims = rel_sims{param_no};
    
    % adjust to only give three waves for age
%     if strcmp(curr_param, 'age')
%         curr_rel_sims = [curr_rel_sims(1); curr_rel_sims(end)];
%     end
    
    if length(curr_rel_sims) < 2 & ~strcmp(curr_param, 'age')
        continue
    end
    
    % extract data for this parameter
    eval(['curr_param_vals = data.config.' curr_param '(curr_rel_sims);']);
    
    % order from high to low
    [curr_param_vals, order] =sort(curr_param_vals, 'descend');
    curr_rel_sims = curr_rel_sims(order);
    
    % make plot of wave vs parameter
    
    % make plots for each wave type and each measurement site
    for wave_type_no = 1 : length(wave_types)
        
        curr_wave_type = wave_types{wave_type_no};
        eval(['rel_wave_sites = rel_sites.' curr_wave_type ';'])
        
        for site_no = 1 : length(rel_wave_sites)
            curr_site = rel_wave_sites{site_no};
            
            % check that this site has been exported
            temp = fieldnames(data.waves);
            if ~sum(~cellfun(@isempty,strfind(temp, curr_site)))
                continue
            end
                
            if strcmp(curr_wave_type, 'PPG')
                artery_name = find_artery_name(curr_site); % use PPG names
            else
                artery_name = curr_site;
            end
            
            %% Figure showing baseline wave
            
            if strcmp(curr_param, 'age')
                
                % Make figure
                if strcmp(curr_wave_type, 'Q') && strcmp(curr_site, 'AorticRoot')
                    paper_size = [400,300];
                    ftsize = 24;
                    axis_ftsize = 20;
                else
                    paper_size = [500,300];
                    ftsize = 32;
                    axis_ftsize = 24;
                end
                figure('Position', [20,20, paper_size])
                
                % Plot each of the waves in turn (corresponding to the three ages)
                rel_sim_no = curr_rel_sims(curr_param_vals == baseline_age);
                
                % extract data to plot
                if strcmp(curr_wave_type, 'Q')
                    eval(['curr_wav.u = data.waves.U_' curr_site '{rel_sim_no};']);
                    eval(['curr_wav.a = data.waves.A_' curr_site '{rel_sim_no};']);
                    curr_wav.v = (1e6)*curr_wav.u.*curr_wav.a;
                    units = 'ml/s'; % converted from 'm3/s';
                elseif strcmp(curr_wave_type, 'A')
                    eval(['curr_wav.v = data.waves.' curr_wave_type '_' curr_site '{rel_sim_no};']);
                    curr_wav.v = (100*100)*curr_wav.v;
                    units = 'cm2'; % converted from 'm3/s';
                else
                    eval(['curr_wav.v = data.waves.' curr_wave_type '_' curr_site '{rel_sim_no};']);
                    eval(['units = data.waves.units.' curr_wave_type ';'])
                end
                units = strrep(units, '2', '^2');
                units = strrep(units, '3', '^3');
                curr_wav.fs = data.waves.fs;
                curr_wav.t = [0:length(curr_wav.v)-1]/curr_wav.fs;
                
                % plot
                plot(curr_wav.t, curr_wav.v, 'k', 'LineWidth', lwidth), hold on
                
                % tidy up
                if ylim == [0 1]
                    ylim([-0.1 1.1])
                end
                xlim([0, max(curr_wav.t)])
                xlabel('Time [s]', 'FontSize', ftsize)
                ylab = ylabel([curr_wave_type, '  [' units ']'], 'FontSize', ftsize);
                set(gca, 'FontSize', axis_ftsize)
                box off
                curr_range = range(curr_wav.v);
                if ~strcmp(curr_wave_type, 'PPG')
                    ylim([min(curr_wav.v)-0.1*curr_range, max(curr_wav.v)+0.1*curr_range])
                else
                    ylim([0 1])
                end
                
                % annotate if aortic flow rate wave
                if strcmp(curr_wave_type, 'Q') && strcmp(curr_site, 'AorticRoot')
                    
                    plot(curr_wav.t(end), curr_wav.v(end), 'ok', 'MarkerSize', 8), 
                    
                    ylim([3*min(curr_wav.v), 1.2*max(curr_wav.v)])
                    
                    % Shade areas
                    rel_els = 1 : (1+find(curr_wav.v(2:end)<0 & curr_wav.v(1:end-1)>=0, 1));
                    patch([curr_wav.t(rel_els), flipud(curr_wav.t(rel_els))],[curr_wav.v(rel_els); zeros(length(rel_els),1)],[1,0.5,0.5], 'LineStyle', 'none')
                    rel_els = rel_els(end) : rel_els(end) + find(curr_wav.v(rel_els(end)+1:end-1)<0 & curr_wav.v(rel_els(end)+2:end)>=0);
                    patch([curr_wav.t(rel_els), flipud(curr_wav.t(rel_els))],[curr_wav.v(rel_els); zeros(length(rel_els),1)],[0.5,0.5,1], 'LineStyle', 'none')
                    
%                     % annotate parameters
%                     dim = [.30 .58 .1 .1];
%                     str = 'SV';
%                     annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'FontSize', ftsize-6);
%                     dim = [.41 .35 .1 .1];
%                     str = 'RFV';
%                     annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'FontSize', ftsize-6);
%                     
%                     annotation('doublearrow',[0.135,0.205],0.85*ones(1,2))
%                     dim = [.135 .85 .1 .1];
%                     str = 'PFT';
%                     annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'FontSize', ftsize-6);
%                     
%                     dim = [.8 .35 .1 .1];
%                     str = '60/HR';
%                     annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'FontSize', ftsize-6);
%                     
%                     annotation('doublearrow',[0.135,0.395],0.25*ones(1,2))
%                     dim = [.22 .24 .1 .1];
%                     str = 'LVET';
%                     annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'FontSize', ftsize-6);
                    
                end
                
                % save
                PrintFigs(gcf, paper_size/70, [PATHS.Analysis_figures, 'waves_', curr_wave_type, '_' artery_name, '_baseline'])
            end
            
            if length(curr_rel_sims) < 2 & strcmp(curr_param, 'age')
                continue
            end
            
            %% Figure showing changes with parameter
            
            % Setup figure
            paper_size = [700,300];
            ftsize = 20;
            figure('Position', [20,20, paper_size])
            leg_labels = {}; max_t = 0;
            
            % Plot each of the waves in turn (corresponding to the three ages)
            colors = linspace(0.3,0.7, length(curr_param_vals))'*ones(1,3);
            colors = [0, 0, 1; 0, 0, 0; 1, 0, 0; 0.5, 0, 0; 0, 0, 0.5; 0.6, 0.6, 0.6];
            temp1 = linspace(0,1, length(curr_param_vals))';
            temp2 = linspace(0,0, length(curr_param_vals))';
            temp3 = linspace(1,0, length(curr_param_vals))';
            colors = [temp1,temp2,temp3];
            min_val = inf; max_val = -inf;
            for sim_no = 1 : length(curr_param_vals)
                rel_sim_no = curr_rel_sims(sim_no);
                
                % extract data to plot
                if ~strcmp(curr_wave_type, 'Q')
                    eval(['curr_wav.v = data.waves.' curr_wave_type '_' curr_site '{rel_sim_no};']);
                    eval(['units = data.waves.units.' curr_wave_type ';'])
                else
                    eval(['curr_wav.u = data.waves.U_' curr_site '{rel_sim_no};']);
                    eval(['curr_wav.a = data.waves.A_' curr_site '{rel_sim_no};']);
                    curr_wav.v = curr_wav.u.*curr_wav.a*1e6;
                    units = 'ml/s';
                end
                
                if plot_model_data_normalised
                    units = 'normalised';
                else
                    units = strrep(units, '2', '^2');
                    units = strrep(units, '3', '^3');
                end
                curr_wav.fs = data.waves.fs;
                curr_wav.t = [0:length(curr_wav.v)-1]/curr_wav.fs;
                                
                % plot
                if ~plot_model_data_normalised
                    plot(curr_wav.t, curr_wav.v, 'Color', colors(sim_no,:), 'LineWidth', lwidth), hold on
                    if strcmp(curr_wave_type, 'U') || strcmp(curr_wave_type, 'Q')
                        plot(curr_wav.t(end), curr_wav.v(end), 'o', 'Color', colors(sim_no,:), 'LineWidth', lwidth, 'MarkerSize', 8,'HandleVisibility','off')
                    end
                    curr_min = min(curr_wav.v);
                    curr_max = max(curr_wav.v);
                    if curr_min < min_val, min_val = curr_min; end
                    if curr_max > max_val, max_val = curr_max; end
                else
                    curr_wav.v = (curr_wav.v - min(curr_wav.v))/range(curr_wav.v);
                    plot(curr_wav.t, curr_wav.v +(offset*(param_val_no-1)), 'Color', colors(sim_no,:), 'LineWidth', lwidth), hold on
                    ylim([-0.1 1.1+(offset*(length(curr_param_vals)-1))]);
                    set(gca, 'YTick', [])
                    min_val = 0; max_val = 1;
                end
                
                % store legend label
                if strcmp(curr_param, 'lvet')
                    temp2 = round(curr_param_vals(sim_no));
                else
                    temp2 = curr_param_vals(sim_no);
                end
                if ~strcmp(curr_param, 'lvet')
                    leg_labels = [leg_labels, num2str(temp2, 2)];
                else
                    leg_labels = [leg_labels, num2str(temp2, 3)];
                end
                max_t = max([max_t, max(curr_wav.t)]);
                
            end
            
            % tidy up
            ylim([min_val - (0.05*(max_val-min_val)), max_val + (0.1*(max_val-min_val))])
            xlim([0, max_t])
            xlabel('Time [s]', 'FontSize', ftsize)
            ylab = ylabel([curr_wave_type, '  [' units ']'], 'FontSize', ftsize);
            set(gca, 'FontSize', ftsize)
            leg_h = legend(leg_labels, 'Location', 'NorthEast');
            [label, units, abbr, graph_title, graph_title_no_units] = make_param_label(curr_param);
            title(leg_h,[strrep(upper(strrep(curr_param, '_', ' ')), 'AGE', 'age'), ' (' units ')'],'FontWeight','Normal');
            clear temp
            temp.upper_els = find(isstrprop(artery_name,'upper'));
            if length(temp.upper_els)>1
                temp.artery_name = [artery_name(1:temp.upper_els(2)-1), ' ', artery_name(temp.upper_els(2):end)];
            else
                temp.artery_name = artery_name;
            end
            title(['Changes in ' temp.artery_name, ' ', curr_wave_type ' with ' strrep(upper(strrep(curr_param, '_', ' ')), 'AGE', 'age')], 'FontSize', ftsize,'FontWeight','Normal')
            clear temp leg_labels
            box off 
            
            % save
            PrintFigs(gcf, paper_size/70, [PATHS.Analysis_figures, 'waves_', curr_wave_type, '_' strrep(artery_name, ' ', '_'), '_vs_', curr_param])
            
        end 
    end
    clear curr_rel_sims res_char_no curr_res_param char_values true_val curr_param_vals rel_col rel_sim_no
end

end

function make_wrist_ppg_parameters_figure(PATHS)

fprintf('\n - Making PPG waves vs parameters figure')

% load data
load(PATHS.exported_data_mat_pwdb_data)

% identify baseline simulation for baseline age
baseline_sim_no = find(data.config.baseline_sim_for_all);
baseline_age = data.config.age(data.config.baseline_sim_for_all);

% for each parameter, identify simulations where only that parameter changed without anything else changing.
all_params = data.config.variations.params;
param_names = data.config.variations.param_names;

% add age to list of parameters
param_names = [{'age'}; param_names];
all_params = [zeros(size(all_params,1),1), all_params];

% find relevant simulations for each parameter
[rel_sims, param_variations] = deal(cell(length(param_names),1));
for param_no = 1 : length(param_names)
    columns_for_other_params = setxor(1:length(param_names), param_no);
    temp = all_params(:,columns_for_other_params);
    if ~strcmp(param_names{param_no}, 'age')
        rel_sims{param_no} = find(~any(temp,2) & data.config.age == baseline_age);
    else
        rel_sims{param_no} = find(~any(temp,2));
    end
    param_variations{param_no} = all_params(rel_sims{param_no},param_no);
end
clear temp

%% Make plots of waves vs parameters

% Setup plotting
rel_sites = {'Radial'};
ftsize = 10;
lwidth = 2;
offset = 0.25;
paper_size = [1000,300];

% Identify parameters I'm interested in
rel_param_els = ones(size(param_names));
for s = 1 : length(param_names)
    curr_param = param_names{s};
    if strcmp(curr_param, 'dbp') || ...
            strcmp(curr_param, 'alpha') || ...
            strcmp(curr_param, 'ht') || ...
            strcmp(curr_param, 'mu') || ...
            strcmp(curr_param, 'p_out') || ...
            strcmp(curr_param, 'p_drop') || ...
            strcmp(curr_param, 'reg_vol') || ...
            strcmp(curr_param, 'rho') || ...
            strcmp(curr_param, 't_pf') || ...
            strcmp(curr_param, 'time_step') || ...
            strcmp(curr_param, 'visco_elastic_log') || ...
            strcmp(curr_param, 'age')
        rel_param_els(s) = 0;
    end
    
    % check there are variations for this parameter
    curr_rel_sims = rel_sims{s};
    if length(curr_rel_sims) ==1
        rel_param_els(s) = 0;
    end
    clear curr_rel_sims
end

% cycle through input parameters
done_a_plot = false;
for site_no = 1 : length(rel_sites)
    counter_no = 0;
    curr_site = rel_sites{site_no};
    
    % make plot of wave vs parameter
    figure('Position', [20,20, paper_size])
    
    for param_no = 1 : length(param_names)
        curr_param = param_names{param_no};
        
        % skip params I'm not interested in
        if ~rel_param_els(param_no)
            continue
        end
        
        % rel sims
        curr_rel_sims = rel_sims{param_no};
        
        if length(curr_rel_sims) < 2
            continue
        end
        
        done_a_plot = true;
        
        % extract data for this parameter
        eval(['curr_param_vals = data.config.' curr_param '(curr_rel_sims);']);
        
        % order from high to low
        [curr_param_vals, order] =sort(curr_param_vals, 'descend');
        curr_rel_sims = curr_rel_sims(order);
        
        counter_no = counter_no+1;
        subplot(2, ceil(sum(rel_param_els)/2), counter_no)
        artery_name = find_artery_name(curr_site); % use PPG names
        leg_labels = {}; max_t = 0;
        units = data.waves.units.PPG;
        
        % Plot each of the waves in turn (corresponding to the three ages)
        ylims = [inf, -inf];
        param_counter_no = 0;
        for param_val_no = 1 : length(curr_param_vals)
            rel_sim_no = curr_rel_sims(param_val_no);
            
            curr_param_val = curr_param_vals(param_val_no);
            if curr_param_val == min(curr_param_vals)
                curr_color = [0 0 1];
            elseif curr_param_val == max(curr_param_vals)
                curr_color = [1,0,0];
            elseif curr_param_val == median(curr_param_vals)
                curr_color = [0,0,0];
            else
                continue
            end
            param_counter_no = param_counter_no+1;
            
            % extract data to plot
            eval(['curr_wav.v = data.waves.PPG_' curr_site '{rel_sim_no};']);
            curr_wav.fs = data.waves.fs;
            curr_wav.t = [0:length(curr_wav.v)-1]/curr_wav.fs;
            
            % plot
            plot(curr_wav.t, curr_wav.v, 'Color', curr_color, 'LineWidth', lwidth), hold on
            
            % store ylim info
            curr_range = range(curr_wav.v);
            ylims = [min([ylims(1), min(curr_wav.v)-0.1*curr_range]), ...
                max([ylims(2), max(curr_wav.v)+0.1*curr_range])];
            
            % store legend label
            if strcmp(curr_param, 'lvet')
                temp2 = round(curr_param_vals(param_val_no));
            else
                temp2 = curr_param_vals(param_val_no);
            end
            leg_labels = [leg_labels, num2str(temp2, 2)];
            max_t = max([max_t, max(curr_wav.t)]);
            
        end
        
        % tidy up
        ylim([-0.1 1.1])
        xlim([0, max_t])
        if counter_no > ceil(sum(rel_param_els)/2)
            xlabel('Time [s]', 'FontSize', ftsize)
            %            elseif ~strcmp(curr_param, 'hr')
            %                set(gca, 'XTick', [])
        end
        %             if counter_no <= length(wave_types)
        %                 title([strrep(strrep(curr_wave_type, '_', ' '), ' WK', ''), '  [' units ']'], 'FontSize', ftsize);
        %             end
        switch curr_param
            case 'dia'
                title_txt = 'Arterial Diameter';
            case 'hr'
                title_txt = 'Heart Rate';
            case 'len'
                title_txt = 'Aortic Length';
            case 'lvet'
                title_txt = 'Duration of Systole';
            case 'mbp'
                title_txt = 'Mean Blood Pressure';
            case 'pvc'
                title_txt = 'Peripheral Compliance';
            case 'pwv'
                title_txt = 'Arterial Stiffness';
            case 'sv'
                title_txt = 'Stroke Volume';
            case 'gamma'
                title_txt = 'Arterial Wall Viscosity';
        end
        title(title_txt, 'FontSize', ftsize);
        
        % y label
%         if counter_no <= ceil(sum(rel_param_els)/2)
%             [label, units, abbr, graph_title] = make_param_label(curr_param);
%             ylab = ylabel({abbr, ['[' units ']']}, 'FontSize', ftsize);
%             set(ylab, 'Units', 'Normalized', 'Position', [-0.17, 0.5, 0]);
%         end
        
        % tidy up
        set(gca, 'FontSize', ftsize, 'YTick', [])
        %legend(leg_labels, 'Location', 'NorthEast'),
        clear leg_labels
        box off
    end
    
    % save
    if done_a_plot
        PrintFigs(gcf, paper_size/70, [PATHS.Analysis_figures, artery_name, '_ppg_waves_vs_params'])
    end
    
    clear curr_rel_sims res_char_no curr_res_param char_values true_val curr_param_vals rel_col rel_sim_no
    
end

end

function make_initial_parameters_waves_figure(PATHS)

fprintf('\n - Making waves vs parameters figure')

% load data
load(PATHS.exported_data_mat_pwdb_data)

% identify baseline simulation for baseline age
baseline_sim_no = find(data.config.baseline_sim_for_all);
baseline_age = data.config.age(data.config.baseline_sim_for_all);

% for each parameter, identify simulations where only that parameter changed without anything else changing.
all_params = data.config.variations.params;
param_names = data.config.variations.param_names;

% add age to list of parameters
param_names = [{'age'}; param_names];
all_params = [zeros(size(all_params,1),1), all_params];

% find relevant simulations for each parameter
[rel_sims, param_variations] = deal(cell(length(param_names),1));
for param_no = 1 : length(param_names)
    columns_for_other_params = setxor(1:length(param_names), param_no);
    temp = all_params(:,columns_for_other_params);
    if ~strcmp(param_names{param_no}, 'age')
        rel_sims{param_no} = find(~any(temp,2) & data.config.age == baseline_age);
    else
        rel_sims{param_no} = find(~any(temp,2));
    end
    param_variations{param_no} = all_params(rel_sims{param_no},param_no);
end
clear temp

%% Make plots of waves vs parameters

% Setup plotting
wave_types = {'U', 'P', 'PPG'};
rel_sites = {'Carotid', 'Radial'};
ftsize = 10;
lwidth = 2;
offset = 0.25;
paper_size = [400,1000];

% cycle through input parameters
no_rel_params = sum(cellfun(@length, param_variations)>=2);
for site_no = 1 : length(rel_sites)
    counter_no = 0;
    curr_site = rel_sites{site_no};
    
    % make plot of wave vs parameter
    figure('Position', [20,20, paper_size])
    done_plot = 0;
    
    for param_no = 1 : length(param_names)
        curr_param = param_names{param_no};
        
        % rel sims
        curr_rel_sims = rel_sims{param_no};
        
        if length(curr_rel_sims) < 2
            continue
        end
        done_plot = 1;
        
        % extract data for this parameter
        eval(['curr_param_vals = data.config.' curr_param '(curr_rel_sims);']);
        
        % order from high to low
        [curr_param_vals, order] =sort(curr_param_vals, 'descend');
        curr_rel_sims = curr_rel_sims(order);
                
        % make plots for each wave type and each measurement site
        for wave_type_no = 1 : length(wave_types)
            
            counter_no = counter_no+1;
            curr_wave_type = wave_types{wave_type_no};
            subplot(no_rel_params, length(wave_types), counter_no)
            
            if strcmp(curr_wave_type, 'PPG')
                artery_name = find_artery_name(curr_site); % use PPG names
            else
                artery_name = curr_site;
            end
            
            leg_labels = {}; max_t = 0;
            
            eval(['units = data.waves.units.' curr_wave_type ';']);
            
            % Plot each of the waves in turn (corresponding to the three ages)
            colors = linspace(0.3,0.7, length(curr_param_vals))'*ones(1,3);
            colors = [0, 0, 1; 0, 0, 0; 1, 0, 0; 0.5, 0, 0; 0, 0, 0.5; 0.4, 0.6, 0.4];
            ylims = [inf, -inf];
            for param_val_no = 1 : length(curr_param_vals)
                rel_sim_no = curr_rel_sims(param_val_no);
               
                eval(['curr_wav.v = data.waves.' curr_wave_type '_' curr_site '{rel_sim_no};']);
                curr_wav.fs = data.waves.fs;
                curr_wav.t = [0:length(curr_wav.v)-1]/curr_wav.fs;
                                
                % plot
                plot(curr_wav.t, curr_wav.v, 'Color', colors(param_val_no,:), 'LineWidth', lwidth), hold on
                if strcmp(curr_wave_type, 'U') && strcmp(curr_param, 'HR')
                    plot(curr_wav.t(end), curr_wav.v(end), 'o', 'Color', colors(param_val_no,:), 'LineWidth', lwidth, 'MarkerSize', 8,'HandleVisibility','off')
                end
                
                % store ylim info
                curr_range = range(curr_wav.v);
                ylims = [min([ylims(1), min(curr_wav.v)-0.1*curr_range]), ...
                    max([ylims(2), max(curr_wav.v)+0.1*curr_range])];
                
                % store legend label
                if strcmp(curr_param, 'lvet')
                    temp2 = round(curr_param_vals(param_val_no));
                else
                    temp2 = curr_param_vals(param_val_no);
                end
                leg_labels = [leg_labels, num2str(temp2, 2)];
                max_t = max([max_t, max(curr_wav.t)]);
                
            end
            
            % tidy up
            if strcmp(curr_wave_type, 'PPG')
                ylim([-0.1 1.1])
            else
                ylim(ylims)
                if strcmp(curr_wave_type, 'U')
                    
                end
            end
            xlim([0, max_t])
            set(gca, 'FontSize', ftsize-4)
            if counter_no >= length(wave_types)*no_rel_params-2 
                xlabel('Time [s]', 'FontSize', ftsize)
            elseif ~strcmp(curr_param, 'hr')
                set(gca, 'XTick', [])
            end
            if counter_no <= length(wave_types)
                 % Title
                 str = [curr_wave_type, '  [' units ']'];
                 dim = [-0.2+counter_no*0.285,0.81,0.3,0.15];
                 annotation('textbox',dim,'String',str,'LineStyle', 'none', 'FontSize', ftsize, 'HorizontalAlignment', 'center');
                
            end
            
            % y label
            if rem(counter_no,length(wave_types)) == 1
                [label, units, abbr, graph_title] = make_param_label(curr_param);
                abbr = strrep(abbr, '_a', '');
                abbr = strrep(abbr, 'Len', 'Length');
                abbr = strrep(abbr, 'Dia', 'Diam.');
                ylab = ylabel(abbr, 'FontSize', ftsize, 'Rotation', 0);
                set(ylab, 'Units', 'Normalized', 'Position', [-0.38, 0.5, 0]);
            end
            
            % tidy up
            %legend(leg_labels, 'Location', 'NorthEast'),
            clear leg_labels
            box off
            
        end
    end
    
    if done_plot
        % save
        PrintFigs(gcf, paper_size/70, [PATHS.Analysis_figures, 'waves_initial_variations_', strrep(artery_name, ' ', '_')])
    end
    
    clear curr_rel_sims res_char_no curr_res_param char_values true_val curr_param_vals rel_col rel_sim_no
    
    
end

end

function make_sensitivity_analysis_figures(PATHS)

fprintf('\n - Making sensitivity analyses figures')

% load data
load(PATHS.exported_data_mat_pwdb_data)

% for each parameter, identify simulations where only that parameter changed without anything else changing.
all_params = data.config.variations.params;

param_names = data.config.variations.param_names;
[rel_sims, param_variations, baseline_logs, baseline_sims, param_variations, param_values] = deal(cell(length(param_names),1));
% cycle through each simulation input parameter
for param_no = 1 : length(param_names)
    
    % identify those columns which correspond to other parameters (i.e. not this one)
    columns_for_other_params = setxor(1:length(param_names), param_no);
    % extract variations corresponding to the other parameters
    temp = all_params(:,columns_for_other_params);
    % identify simulations in which none of these other parameters were varied from baseline values 
    rel_sims{param_no} = find(~any(temp,2));
    % skip this parameter if it was not varied at all from its baseline value
    if length(rel_sims{param_no}) <= sum(data.config.baseline_sim_for_age)
        continue
    end
    
    % use all the simulations
    rel_sims{param_no} = 1:length(data.config.age);
    % setup variables
    temp_baseline_vals = nan(length(rel_sims{param_no}),1);
    temp_baseline_logs = false(length(rel_sims{param_no}),1);
    % cycle through each simulation
    for s = 1 : length(rel_sims{param_no})
        % identify age of this simulation
        curr_sim = rel_sims{param_no}(s);
        curr_age = data.config.age(curr_sim);
        % extract baseline value of this parameter
        temp_rel_baseline_age_sim = find(data.config.age == curr_age & data.config.baseline_sim_for_age);
        eval(['temp_baseline_vals(s) = data.config.' param_names{param_no} '(temp_rel_baseline_age_sim);']);
        % store the baseline simulation for this age
        if temp_rel_baseline_age_sim == curr_sim
            temp_baseline_logs(s) = true;
        end
        temp_baseline_sims(s) = temp_rel_baseline_age_sim;
    end
    % store details of the baseline simulations corresponding to each simulation
    baseline_logs{param_no} = temp_baseline_logs; clear temp_baseline_logs
    baseline_sims{param_no} = temp_baseline_sims; clear temp_baseline_sims
    baseline_vals{param_no} = temp_baseline_vals; clear temp_baseline_vals
    % store the variations for this parameter (in number of SDs from age-specific mean)
    param_variations{param_no} = all_params(rel_sims{param_no},param_no);
    % store the values of this parameter
    eval(['param_values{param_no} = data.config.' param_names{param_no} '(rel_sims{param_no});']);
end

%% sensitivity analyses

for s = 1 : length(data.haemods)
    data.haemods(s).trial = (data.haemods(s).SBP_b - data.haemods(s).DBP_b);
end
for s = 1 : length(data.haemods)
    data.haemods(s).trial2 = (data.haemods(s).SBP_b - data.haemods(s).DBP_b)./(data.haemods(s).P2pk_a-data.haemods(s).DBP_a);
end
for s = 1 : length(data.haemods)
    data.haemods(s).trial3 = (data.haemods(s).SBP_b - data.haemods(s).DBP_b)./(data.haemods(s).P1in_a-data.haemods(s).DBP_a);
end

% Simulated characteristics
resultant_params = {'trial', 'trial2', 'trial3', 'CO', 'LVET', 'AP_a', 'AI_a', 'P1in_a', 'P2in_a', 'P1pk_a', 'P2pk_a', 'P1in_c', 'P2in_c', 'P1pk_c', 'P2pk_c', 'PP_a', 'SMBP_a', 'MBP_a', 'SBP_diff', 'SBP_a', 'PP_amp', 'IAD', 'SI', 'RI', 'AGI_mod', 'MBP_drop_finger', 'MBP_drop_ankle', 'SBP_b', 'DBP_b', 'PP_b', 'MBP_b', 'AP_c', 'AI_c', 'Tr_a'};

for do_all = 0:1
    
    I = nan(length(param_names), length(resultant_params));
    for param_no = 1 : length(param_names)
        % skip if this parameter wasn't varied indpendently
        if length(param_variations{param_no}) <= 1 % (1 accounts for the baseline sim)
            I(param_no,:) = nan;
            continue
        end
        
        % extract data for this parameter
        curr_rel_sims = rel_sims{param_no};
        curr_rel_variations = param_variations{param_no};
        curr_baseline_vals = baseline_vals{param_no};
        curr_baseline_sim_nos = baseline_sims{param_no};
        curr_baseline_logs = baseline_logs{param_no};
        
        % perform sensitivity analysis
        for res_param_no = 1 : length(resultant_params)
            
            % Extract data for this resultant param
            curr_res_param = resultant_params{res_param_no};
            eval(['all_values = extractfield(data.haemods, ''' curr_res_param ''');']);
            all_values = all_values(:);
            
            % identify relevant simulations (excluding those which were physiologically implausible)
            if do_all
                rel_els = ~curr_baseline_logs & curr_rel_variations~=0 & ~isnan(all_values);
            else
                rel_els = ~curr_baseline_logs & curr_rel_variations~=0 & ~isnan(all_values) & data.plausibility.plausibility_log;
            end
            rel_value_sims = curr_rel_sims(rel_els);
            rel_baseline_sims = curr_baseline_sim_nos(rel_els);
            rel_variations = curr_rel_variations(rel_els);
            
            % identify param values
            res_param_values = all_values(rel_value_sims);
            res_param_baseline_values = all_values(rel_baseline_sims);
            
            % calculate index
            v = rel_variations;
            I(param_no,res_param_no) = 100*mean((res_param_values - res_param_baseline_values)./(abs(res_param_baseline_values).*v));
        end
    end
    
    % Make sensitivity plots for each resultant parameter
    paper_size = [200,200,700,400];
    for res_param_no = 1 : length(resultant_params)
        
        curr_i_vals = I(:, res_param_no);
        curr_param_names = param_names;
        rel_els = ~strcmp(curr_param_names, 'ht') & ~strcmp(curr_param_names, 'rho');
        curr_i_vals = curr_i_vals(rel_els);
        curr_param_names = curr_param_names(rel_els);
        %req_order = {'hr','sv','lvet','t_pf','reg_vol','mbp','pwv','p_out','dia','len','pvc'};
        
        ref_param = resultant_params{res_param_no};
        curr_param_names = strrep(curr_param_names, 'reg_vol', 'Rvol');
        curr_param_names = strrep(curr_param_names, 'p_out', 'Pout');
        
        % ignore parameters which we're not interested in
        rel_no = strcmp(curr_param_names, 'dbp');
        curr_i_vals(rel_no) = nan;
        
        fig_done = make_sens_analysis_fig(curr_param_names, curr_i_vals, 'norm', paper_size, ref_param);
        
        if fig_done ==1
            if do_all
                filename = [PATHS.Analysis_figures, 'sens_', ref_param, '_all'];
            else
                filename = [PATHS.Analysis_figures, 'sens_', ref_param, '_plausible'];
            end
            PrintFigs(gcf, paper_size(3:4)/70, filename)
        end
        
    end
    
end

end

function fig_done = make_sens_analysis_fig(model_params, i_vals, type, paper_size, ref_param)

% exclude variables which weren't varied
rel_els = ~isnan(i_vals);
i_vals = i_vals(rel_els);
model_params = model_params(rel_els);

if isempty(model_params)
    fig_done = 0;
    return
else
    fig_done = 1;
end

% re-arrange
rel_order = {'hr','sv','lvet','t_pf','Rvol','dia','len','pwv','gamma','mbp','pvc'};
rel_order = {'hr','sv','lvet','t_pf','Rvol','dia','len','pwv','mbp','pvc'};
rel_order = {'hr', 'sv', 'lvet', 'dia', 'pwv', 'mbp'};
for s = 1 : length(rel_order)
    order(s) = find(strcmp(model_params,rel_order{s}));    
end
model_params = model_params(order);
i_vals = i_vals(order);

% re-name
model_params = strrep(model_params, 'hr', 'Heart Rate');
model_params = strrep(model_params, 'sv', 'Stroke Volume');
model_params = strrep(model_params, 'lvet', 'Duration Systole');
model_params = strrep(model_params, 't_pf', 'Time to Peak Flow');
model_params = strrep(model_params, 'Rvol', 'Regurgitation Volume');
model_params = strrep(model_params, 'dia', 'Large Art. Diameter');
model_params = strrep(model_params, 'len', 'Prox. Aortic Length');
model_params = strrep(model_params, 'pwv', 'Pulse Wave Velocity');
model_params = strrep(model_params, 'mbp', 'Periph. Vasc. Res.');
model_params = strrep(model_params, 'pvc', 'Periph. Vasc. Comp.');
model_params = strrep(model_params, 'Pout', 'Outflow Pressure');

if strcmp(type, 'abs')
    ylims = [-30, 75];
    i_vals = abs(i_vals);
    ylab_txt = 'abs(I) [%]';
else
    ylims = [-100, 100];
    ylab_txt = 'I [%]';
end

model_params = strrep(model_params, '_', ' ');

figure('Position', paper_size)
subplot('Position',[0.13,0.30,0.86,0.63])
ftsize = 20;
bar(i_vals)
set(gca,'XTickLabel',model_params)
set(gca, 'FontSize', ftsize)
ylab = ylabel(ylab_txt, 'Rotation', 0, 'FontSize', ftsize);
set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
%ylim(ylims)
% [label, units, abbr, graph_title, graph_title_no_units] = make_param_label(ref_param);
% title(graph_title_no_units, 'FontSize', ftsize)
xtickangle(30);
xlim([0.5, length(model_params)+0.5])

% annotations
y_val = 0.73;
dim = [.65   y_val .3 .3];
str = {'Arterial Network'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',ftsize);

dim = [.19 y_val .3 .3];
str = {'Cardiac Properties'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',ftsize);

hold on
ylims = [min(i_vals) - 0.1*range(i_vals), max(i_vals) + 0.22*range(i_vals)];
plot(3.5*ones(1,2), [-100 100], '--k')
ylim(ylims)

end

function rel_ylims = plot_sds(ages, means, sds, no_sd, req_color, up)

do_fill = 0;

temp1 = means - (no_sd*sds);
temp2 = means + (no_sd*sds);

if do_fill
    temp = [temp1(:)', fliplr(temp2(:)')];
    fill_vals.y = temp;
    temp = [ages(:)', fliplr(ages(:)')];
    fill_vals.x = temp;
    clear temp
    h = fill(fill_vals.x, fill_vals.y, 'r'); hold on
    h.FaceColor = req_color;
    h.EdgeColor = 'none';
else
    plot(ages, temp1, '--k', 'LineWidth', 2), hold on
    plot(ages, temp2, '--k', 'LineWidth', 2)
end

min_val = min(temp1);
max_val = max(temp2);
range_val = up.ylim_offset*(max_val - min_val);
rel_ylims = [min_val-range_val, max_val+range_val];

end

function make_wave_speed_figure(PATHS)

fprintf('\n - Making wave speed figure')

% load collated data
load(PATHS.exported_data_mat_pwdb_data)

% Find relevant simulations: (i) baseline simulation for each age
%                           (ii) variations in PWV for baseline age
rel_sims.ages = find(data.config.baseline_sim_for_age);
baseline_age = data.config.age(data.config.baseline_sim_for_all);
non_pwv_cols = ~strcmp(data.config.variations.param_names, 'pwv'); non_pwv_cols = non_pwv_cols(:)';
rel_sims.pwvs = find(data.config.baseline_sim_for_all | ...
    (data.config.age == baseline_age & sum(abs(data.config.variations.params(:,non_pwv_cols)),2) == 0) );

% Setup figure for wave speed curves
paper_size = [800, 600];
fig_settings.lwidth = 2;
fig_settings.ftsize = 20;
fig_settings.colors = [1,1,1];

% cycle through each plot type
plot_types = fieldnames(rel_sims);
for plot_no = 1 : length(plot_types)
    
    figure('Position', [20,20,paper_size])
    labels = {};

    % cycle through each relevant simulation
    eval(['rel_sims.curr = rel_sims.' plot_types{plot_no} ';'])
    for sim_no = 1 : length(rel_sims.curr)
        curr_sim_no = rel_sims.curr(sim_no);
        
        % calculate wave speeds
        rel_network_spec.inlet_radius = data.config.network.inlet_radius(curr_sim_no,:);
        rel_network_spec.outlet_radius = data.config.network.outlet_radius(curr_sim_no,:);
        rel_network_spec.length = data.config.network.length(curr_sim_no,:); rel_network_spec.length = rel_network_spec.length(:);
        k = data.config.constants.k(curr_sim_no,:);
        rho = data.config.constants.rho(curr_sim_no);
        [wave_speed, ave_radius] = calculate_theoretical_wave_speeds(rel_network_spec, k, rho);
        clear k rho rel_network_spec
        [ave_radius, order] = sort(ave_radius);
        wave_speed.all = wave_speed.all(order);
        
        % plot signal
        color_intensity = 0.2+0.5*(sim_no/length(rel_sims.curr));
        curr_color = color_intensity*fig_settings.colors;
        plot(2*ave_radius*1000, wave_speed.all, 'o-', 'Color', curr_color, 'LineWidth', fig_settings.lwidth, 'MarkerSize', 8, 'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color)
        hold on
        
        clear wave_speed ave_radius curr_color color_intensity
        
        if strcmp(plot_types{plot_no}, 'ages')
            labels{end+1} = num2str(data.config.age(curr_sim_no));
        else
            labels{end+1} = num2str(data.config.pwv_SD(curr_sim_no));
        end
        
    end
    
    % tidy up
    set(gca, 'FontSize', fig_settings.ftsize)
    ylab = ylabel('Wave Speed [m/s]', 'FontSize', fig_settings.ftsize);
    %set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
    xlab = xlabel('Diameter [mm]', 'FontSize', fig_settings.ftsize);
    box off
    legend(labels, 'FontSize', fig_settings.ftsize-2), clear labels
    
    % Save figure
    PrintFigs(gcf, paper_size/70, [PATHS.wave_speed_fig, '_', plot_types{plot_no}])
    
end


end

function [wave_speeds, ave_radius] = calculate_theoretical_wave_speeds(rel_network_spec, k, rho)

% find current wave speeds
ave_radius = mean([rel_network_spec.inlet_radius(:), rel_network_spec.outlet_radius(:)], 2);
ave_radius_cm = 100*ave_radius;

% - New version
wave_speed = empirical_wave_speed(ave_radius_cm, k, rho);
wave_speeds.all = wave_speed;

% extract pwvs

% - cf (noting that carotid measurements are taken half way along carotid)
lens_aorta_carotid = [rel_network_spec.length([1,2]); rel_network_spec.length(15)/2];
carotid_radius_cm = 100* (rel_network_spec.inlet_radius(15) - (0.5*(rel_network_spec.inlet_radius(15) - rel_network_spec.outlet_radius(15))));
path_speeds_aorta_carotid = [wave_speed([1,2]); empirical_wave_speed(carotid_radius_cm, k, rho)];

lens_aorta_femoral = rel_network_spec.length([1,2,14,18,27,28,35,37,39,41,42,44]);
path_speeds_aorta_femoral = wave_speed([1,2,14,18,27,28,35,37,39,41,42,44]);

path_len = sum(lens_aorta_femoral)-sum(lens_aorta_carotid);
time_taken = sum(lens_aorta_femoral./path_speeds_aorta_femoral) - sum(lens_aorta_carotid./path_speeds_aorta_carotid);
wave_speeds.aorta = path_len/time_taken;

% - arm (noting that brachial measurements are taken half way along brachial)
lens_aorta_radial = rel_network_spec.length([1,2,14,19,21,22]);
path_speeds_aorta_radial = wave_speed([1,2,14,19,21,22]);
lens_aorta_brachial = [rel_network_spec.length([1,2,14,19]); rel_network_spec.length(21)/2];
brachial_radius_cm = 100* (rel_network_spec.inlet_radius(21) - (0.5*(rel_network_spec.inlet_radius(21) - rel_network_spec.outlet_radius(21))));
path_speeds_aorta_brachial = [wave_speed([1,2,14,19]); empirical_wave_speed(brachial_radius_cm, k, rho)];

path_len = sum(lens_aorta_radial)-sum(lens_aorta_brachial);
time_taken = sum(lens_aorta_radial./path_speeds_aorta_radial) - sum(lens_aorta_brachial./path_speeds_aorta_brachial);
wave_speeds.arm = path_len/time_taken;

% - leg
ankle_els = [1,2,14,18,27,28,35,37,39,41,42,44,46,49];
lens_aorta_ankle = rel_network_spec.length(ankle_els);
path_speeds_aorta_ankle = wave_speed(ankle_els);
femoral_els = [1,2,14,18,27,28,35,37,39,41,42,44];
lens_aorta_femoral = rel_network_spec.length(femoral_els);
path_speeds_aorta_femoral = wave_speed(femoral_els);

path_len = sum(lens_aorta_ankle)-sum(lens_aorta_femoral);
time_taken = sum(lens_aorta_ankle./path_speeds_aorta_ankle) - sum(lens_aorta_femoral./path_speeds_aorta_femoral);
wave_speeds.leg = path_len/time_taken;

end

function wave_speed = empirical_wave_speed(ave_radius_cm, k, rho)

Eh_D0 = (k(1)*exp(k(2)*ave_radius_cm))+k(3); % Eh/D0 (from Mynard's 2015 paper, eqn 3)
c0_squared = (2/3)*(Eh_D0/(rho/1000)); % from Mynard's 2015 paper, eqn 3, converts rho from kg/m3 to g/cm3
wave_speed = sqrt(c0_squared)/100;  % converts from cm/s to m/s.

end

function make_pw_along_path_figure(PATHS)

% load path data
if ~exist(PATHS.exported_data_mat_pwdb_data_w_aorta_finger_path, 'file')
    return
else
    fprintf('\n - Plotting pulse waves along particular paths')
end

do_fid_pts = 1;
do_brach_color = 1;  % whether or not to highlight the brachial artery
%plot_f_b = 0;  % whether or not to plot forward and backward waves

load(PATHS.exported_data_mat_pwdb_data)
ages = [min(data.config.age), max(data.config.age)];

for age_no = 1:length(ages)
    
    curr_age = ages(age_no);
    
    % identify baseline simulation
    baseline_sim_no = data.config.baseline_sim_for_age & data.config.age == curr_age;
    
    % Setup figure for Pressure and Flow vel
    paper_size = [600, 400];
    fig_settings.lwidth = 2;
    fig_settings.ftsize = 18;
    paths = {'arm', 'leg'};
    fig_settings.sig_types = {'P','U'};
    fig_settings.colors = {[1,0,0], [1,0,0]};
    fig_settings.ylims = {[0,140], [-0.1, 1.5]};
    
    for path_no = 1 : length(paths)
        curr_path = paths{path_no};
        switch curr_path
            case 'arm'
                curr_path_var_name = 'aorta_finger';
            case 'leg'
                curr_path_var_name = 'aorta_foot';
        end
        
        for sig_type_no = 1 : length(fig_settings.sig_types)
            curr_sig_type = fig_settings.sig_types{sig_type_no};
            
            % Load relevant data
            switch curr_path
                case 'arm'
                    load(PATHS.exported_data_mat_pwdb_data_w_aorta_finger_path)
                case 'leg'
                    switch curr_sig_type
                        case 'P'
                            load(PATHS.exported_data_mat_pwdb_data_w_aorta_foot_path_p)
                        case 'U'
                            load(PATHS.exported_data_mat_pwdb_data_w_aorta_foot_path_u)
                        case 'A'
                            load(PATHS.exported_data_mat_pwdb_data_w_aorta_foot_path_a)
                    end
            end
            
            
            % extract relevant data
            eval(['rel_data = data.path_waves.' curr_path_var_name '(baseline_sim_no);']);
            eval(['rel_data.sig = rel_data.' curr_sig_type ';']);
            rel_data.p1_el = data.pw_inds.AorticRoot_P1in_T(baseline_sim_no)*data.path_waves.fs+1;
            rel_data.p2_el = data.pw_inds.AorticRoot_P2pk_T(baseline_sim_no)*data.path_waves.fs+1;
            
            figure('Position', [20,20, paper_size])
            if length(rel_data.sig) > 50
                inc = 2;
            else
                inc = 1;
            end
            for s = 1:inc:length(rel_data.sig)
                if do_brach_color && strcmp(rel_data.artery{s}, 'brachial')
                    curr_color = 'b';
                else
                    curr_color = 'k';
                end
                plot3([0:length(rel_data.sig{s})-1]/data.path_waves.fs, ones(length(rel_data.sig{s}),1)*rel_data.dist(s), rel_data.sig{s}, curr_color), hold on
            end
            
            if do_brach_color && strcmp(curr_path, 'arm')
                dist_at_start = min(rel_data.dist(strcmp(rel_data.artery, 'brachial')));
                dist_at_end = max(rel_data.dist(strcmp(rel_data.artery, 'brachial')));
                prop_brachial = 0.75;
                desired_dist = dist_at_start + prop_brachial*(dist_at_end-dist_at_start);
                [~,s] = min(abs(rel_data.dist-desired_dist));
                plot3([0:length(rel_data.sig{s})-1]/data.path_waves.fs, ones(length(rel_data.sig{s}),1)*rel_data.dist(s), rel_data.sig{s}, 'b', 'LineWidth', 2), hold on
            end
            
            % Tidy up
            xlabel('Time (s)', 'FontSize', fig_settings.ftsize)
            ylabel('Distance (cm)', 'FontSize', fig_settings.ftsize)
            if strcmp(curr_sig_type, 'P')
                zlabel('Pressure (mmHg)', 'FontSize', fig_settings.ftsize)
            elseif strcmp(curr_sig_type, 'U')
                zlabel('Flow velocity (m/s)', 'FontSize', fig_settings.ftsize)
            elseif strcmp(curr_sig_type, 'Q')
                zlabel('Flow rate (m^3/s)', 'FontSize', fig_settings.ftsize)
            end
            set(gca, 'FontSize', fig_settings.ftsize -2 )
            view(15, 21)
            %view(0, 21)
            grid on
            
            eval(['xlim([0 length(rel_data.' curr_sig_type '{s})/data.path_waves.fs])'])
            ylim([0 max(rel_data.dist)])
            if strcmp(curr_sig_type, 'P')
                zlim([65, 130])
            end
            
            % annotate fiducial points
            if do_fid_pts && strcmp(curr_sig_type, 'P')
                plot3(rel_data.p1_el/data.path_waves.fs, rel_data.dist(1), rel_data.sig{1}(rel_data.p1_el), 'or', 'MarkerFaceColor', 'r')
                plot3(rel_data.p2_el/data.path_waves.fs, rel_data.dist(1), rel_data.sig{1}(rel_data.p2_el), 'or', 'MarkerFaceColor', 'r')
                if age_no == 1
                    dim1 = [.13 .43 .1 .1];
                    dim2 = [.28 .32 .1 .1];
                else
                    dim1 = [.13 .43 .1 .1];
                    dim2 = [.29 .43 .1 .1];
                end
                str = 'P1_a';
                annotation('textbox',dim1,'String',str,'FitBoxToText','on','LineStyle','none','Color','r','FontSize', fig_settings.ftsize);
                str = 'P2_a';
                annotation('textbox',dim2,'String',str,'FitBoxToText','on','LineStyle','none','Color','r','FontSize', fig_settings.ftsize);
            end
            
            % Save figure
            if age_no == 1
                savepath = [PATHS.pw_path_fig, '_', curr_path, '_' curr_sig_type '_young'];
            else
                savepath = [PATHS.pw_path_fig, '_', curr_path, '_' curr_sig_type '_elderly'];
            end
            savefig(savepath)
            PrintFigs(gcf, paper_size/70, savepath)
            
            clear curr_color curr_sig s
        end
        
    end
    close all
    
end

end

function plot_baseline_signals(fig_settings, data)

% identify baseline simulation data
baseline_sim_no = find(data.config.baseline_sim_for_all);

% cycle through different signals
for sig_type_no = 1 : length(fig_settings.sig_types)
    curr_sig_type = fig_settings.sig_types{sig_type_no};
    
    % cycle through different sites
    for req_site_no = 1 : length(fig_settings.req_sites)
        curr_site = fig_settings.req_sites{req_site_no};
        
        % extract relevant signal at this site
        eval(['rel_sig.v = data.waves.' curr_sig_type '_' curr_site '{baseline_sim_no};'])
        rel_sig.fs = data.waves.fs;
        rel_sig.t = [0:length(rel_sig.v)-1]/rel_sig.fs;
        
        % convert to friendly units if needed
        if sum(strcmp(curr_sig_type, {'P', 'Pe'}))
            rel_sig.units = 'mmHg';
        elseif strcmp(curr_sig_type, 'A')
            rel_sig.v = rel_sig.v*1000*1000;   % m^2 to mm^2
            rel_sig.units = 'mm^2';
        elseif strcmp(curr_sig_type, 'PPG')
            rel_sig.units = 'au';
        elseif strcmp(curr_sig_type, 'U')
            rel_sig.units = 'm/s';
        end
        
        % setup subplot
        if sig_type_no == 1
            subplot(length(fig_settings.req_sites),2,(2*req_site_no)-1)
        else
            subplot(length(fig_settings.req_sites),2,(2*req_site_no))
        end
        
        % plot signal
        plot(rel_sig.t, rel_sig.v, 'Color', fig_settings.colors{sig_type_no}, 'LineWidth', fig_settings.lwidth)
        
        % tidy up
        set(gca, 'XTick', 0:0.2:1)
        set(gca, 'FontSize', fig_settings.ftsize,'XAxisLocation', 'origin')
        xlim([0, rel_sig.t(end)])
        if req_site_no == 1
            dim = [.2+.5*(sig_type_no-1) .77 .2 .2];
            str = [strrep(curr_sig_type, '_WK', ''), ' [', rel_sig.units, ']'];
            annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle', 'none', 'FontSize', fig_settings.ftsize);
%             title(str, 'FontSize', fig_settings.ftsize);
        end
        if req_site_no < length(fig_settings.req_sites)
            set(gca, 'XTickLabel', [])
        end
        if sum(strcmp(fieldnames(fig_settings), 'ylims'))
            ylims = fig_settings.ylims{sig_type_no};
            if ~isnumeric(ylims)
                temp_range = range(rel_sig.v);
                ylims = [min(rel_sig.v)-0.1*temp_range, max(rel_sig.v)+0.1*temp_range];
            end
            ylim(ylims);
        end
        if req_site_no == length(fig_settings.req_sites)
            xlab = xlabel('Time [s]', 'FontSize', fig_settings.ftsize);
            ylims = ylim;
            if ylims(1)<0
                set(xlab, 'Units', 'Normalized', 'Position', [0.6, -0.6, 0]);
            end
        end
        
        box off
        
    end
    
end


end

function comparison_w_mynard_waves_figure(PATHS)

% see if the data is available
loadpath = '/Users/petercharlton/Google Drive/Work/Code/nektar/Mynard2015_data/ABME2015data.mat';
if ~exist(loadpath, 'file')
    return
end

fprintf('\n - Making Figure to compare waves with those from Mynard 2015')

% load data
load(PATHS.exported_data_mat_pwdb_data)

% identify baseline simulation
baseline_sim_no = find(data.config.baseline_sim_for_all);

% Import waves from Mynard's article
load(loadpath);

% Setup plotting
wave_types = {'A', 'Q','P'};
rel_sites =  {'AorticRoot', 'Carotid', 'Brachial', 'Radial', 'CommonIliac', 'Femoral', 'AntTibial', 'SupTemporal'};
ftsize = 24;
lwidth = 2;
plot_model_data_normalised = false;
offset = 0.25;

% make plots for each wave type and each measurement site
for wave_type_no = 1 : length(wave_types)
    
    curr_wave_type = wave_types{wave_type_no};
    
    for site_no = 1 : length(rel_sites)
        curr_site = rel_sites{site_no};
        
        if strcmp(curr_wave_type, 'PPG')
            artery_name = find_artery_name(curr_site); % use PPG names
        else
            artery_name = curr_site;
        end
        
        paper_size = [600,400];
        
        % Make figure
        figure('Position', [20,20, paper_size])
        
        max_t = 0;
        
        if strcmp(curr_wave_type, 'Q')
            units = 'm^3/s';
        elseif strcmp(curr_wave_type, 'P')
            units = 'mmHg';
        elseif strcmp(curr_wave_type, 'A')
            units = 'mm';
        end
        
        if plot_model_data_normalised
            units = 'normalised';
        end
        
        % extract data to plot
        if strcmp(curr_wave_type, 'Q')
            eval(['temp1 = data.waves.U_', curr_site, '{baseline_sim_no};']);
            eval(['temp2 = data.waves.A_', curr_site, '{baseline_sim_no};']);
            curr_wav.v = temp1.*temp2;
        else
            eval(['curr_wav.v = data.waves.' curr_wave_type '_', curr_site, '{baseline_sim_no};']);
        end
        curr_wav.fs = data.waves.fs;
        curr_wav.t = [0:length(curr_wav.v)-1]/curr_wav.fs;
        
        % convert units if necessary
        if strcmp(curr_wave_type, 'A')
            curr_wav.v= 1000*sqrt(curr_wav.v/pi);
        end
        
        % plot
        if ~plot_model_data_normalised
            plot(curr_wav.t, curr_wav.v, 'b', 'LineWidth', lwidth), hold on
        else
            curr_wav.v = (curr_wav.v - min(curr_wav.v))/range(curr_wav.v);
            plot(curr_wav.t, curr_wav.v + offset, 'b', 'LineWidth', lwidth), hold on
            ylim([-0.1 1.1+offset]);
            set(gca, 'YTick', [])
        end
        
        % store legend label
        max_t = max([max_t, max(curr_wav.t)]);
        
        % tidy up
        xlim([0, max_t])
        xlabel('Time [s]', 'FontSize', ftsize)
        ylab = ylabel([strrep(curr_wave_type, '_', ' '), '  [' units ']'], 'FontSize', ftsize);
        set(gca, 'FontSize', ftsize)
        
        % Plot waves from literature
        
        % extract relevant wave from Mynard's data
        lit_wave.t = a115lb.t;
        lit_wave.fs = 1000;
        switch curr_site
            case 'AorticRoot'
                mynard_name = 'AoRt';
            case 'Carotid'
                mynard_name = 'Lcar';
            case 'Brachial'
                mynard_name = 'LBrach';
            case 'Radial'
                mynard_name = 'LRadI';
            case 'CommonIliac'
                mynard_name = 'LcmIlc';
            case 'Femoral'
                mynard_name = 'Lfem';
            case 'AntTibial'
                mynard_name = 'LATib';
            case 'SupTemporal'
                mynard_name = 'LSupTemp';
        end
        wave_el = find(strcmp(a115lb.monitor.name, mynard_name));
        eval(['lit_wave.v = a115lb.tnode.' strrep(lower(curr_wave_type), 'a', 'A') '(:,wave_el);']);
        if strcmp(curr_wave_type, 'P') || strcmp(curr_wave_type, 'A')
            [~,temp] = min(lit_wave.v);
            if temp > 799
                lit_wave.v = lit_wave.v(temp-799:temp);
            elseif temp < 80
                lit_wave.v = lit_wave.v(temp:temp+799);
            else
                len = length(lit_wave.v);
                lit_wave.v = [lit_wave.v(temp:end); lit_wave.v(1:799-(len-temp+1))];
            end
            clear temp
            
        end
        
            % convert units if necessary
        if strcmp(curr_wave_type, 'P')
            lit_wave.v = lit_wave.v/1333.3;
        elseif strcmp(curr_wave_type, 'Q')
            lit_wave.v = lit_wave.v/1000000;
            [r,lags] = xcorr(curr_wav.v, lit_wave.v);
            [~, rel_lag] = max(r);
            rel_lag = -1*(rel_lag-length(lit_wave.v));
            lit_wave.v = [lit_wave.v(rel_lag:end); lit_wave.v(1:rel_lag-1)];
        elseif strcmp(curr_wave_type, 'A')
            lit_wave.v = lit_wave.v/10000;   
            lit_wave.v= 1000*sqrt(lit_wave.v/pi);
        end
            
        % plot
        plot([0:length(lit_wave.v)-1]/lit_wave.fs, lit_wave.v, 'r', 'LineWidth', lwidth), hold on
        
        max_t = max([max_t, max(curr_wav.t)]);
        
        % tidy up
        xlim([0, max_t])
        xlabel('Time [normalised]', 'FontSize', ftsize)
        ylab = ylabel([curr_wave_type, ' [' units ']'], 'FontSize', ftsize);
        set(gca, 'FontSize', ftsize)
        title(artery_name, 'FontSize', ftsize+4)
        
        % save
        PrintFigs(gcf, paper_size/70, [PATHS.Analysis_figures, 'Mynard_vs_', curr_wave_type, '_', strrep(artery_name, ' ', '_')])
        
    end


end

end

function make_changes_in_waves_age_figure(PATHS)

fprintf('\n - Making Figure to show changes in wave shapes with age')


% - Import waves from literature
% Path for loading waves from literature
lit_waves_loadpath = '/Users/petercharlton/Google Drive/Work/Projects/PWV Project/Reports/Pulse_wave_database/Figures/literature_pulse_waves/waves_from_literature.mat';
if ~exist(lit_waves_loadpath, 'file')
    % if this file doesn't exist then insert dummy data
    waves(1).author = 'abc'; waves(1).site = 'abc'; waves(1).wave_type = 'abc'; waves(1).v.y = []; waves(1).v.fs = nan; waves(1).age = -1;
else
    % otherwise load waves from the literature
    load(lit_waves_loadpath);
end

% load data
load(PATHS.exported_data_mat_pwdb_data)

% identify baseline simulation
baseline_sim_nos = find(data.config.baseline_sim_for_age);
ages = data.config.age(baseline_sim_nos);

% Extract relevant data for min, max, and mid ages
[~,min_el] = min(ages);
mid_el = ceil(0.5*(1+length(ages)));
[~,max_el] = max(ages);
rel_sims = [min_el, mid_el, max_el]; clear min_el mid_el max_el
rel_ages = ages(rel_sims); clear ages

% Setup plotting
wave_types = {'PPG', 'P', 'U'};
rel_sites.P =    {'Brachial', 'AorticRoot', 'Carotid', 'Radial', 'CommonIliac', 'Femoral', 'AntTibial', 'SupTemporal', 'Digital'};
rel_sites.PPG =    {'Digital', 'AntTibial', 'SupTemporal', 'Radial'};
if sum(strcmp(fieldnames(data.waves), 'U_SupMidCerebral'))
    rel_sites.U = {'SupMidCerebral', 'AorticRoot', 'Carotid'};
else
    rel_sites.U = {'AorticRoot', 'Carotid'};
end

ftsize = 28;
lwidth = 2;
age_colors = linspace(0.4,0.7, length(rel_ages));
plot_model_data_normalised = true;

% make plots for each wave type and each measurement site
for wave_type_no = 1 : length(wave_types)
    
    curr_wave_type = wave_types{wave_type_no};
    eval(['rel_wave_sites = rel_sites.' curr_wave_type ';'])
    
    for site_no = 1 : length(rel_wave_sites)
        curr_site = rel_wave_sites{site_no};
        if strcmp(curr_site, 'AorticRoot') || strcmp(curr_site, 'Carotid')
            offset = 0.4;
        else
            offset = 0.4; %0.25;
        end
        if strcmp(curr_wave_type, 'PPG')
            artery_name = find_artery_name(curr_site); % use PPG names
        else
            artery_name = curr_site;
        end
        
        % Determine whether there is a corresponding set of waves from the literature
        rel_waves = false(length(waves),1);
        for s = 1 : length(waves)
            if strcmp(lower(waves(s).site), lower(artery_name)) && strcmp(waves(s).wave_type, curr_wave_type)
                rel_waves(s) = true;
            end
            if strcmp(waves(s).site, 'aortic') && strcmp('AorticRoot', artery_name) && strcmp(waves(s).wave_type, curr_wave_type)
                rel_waves(s) = true;
            end
            if strcmp(waves(s).site, 'cerebral') && strcmp('SupMidCerebral', artery_name) && strcmp(waves(s).wave_type, curr_wave_type)
                rel_waves(s) = true;
            end
        end
        clear s
        rel_wave = find(rel_waves); clear rel_waves
        if ~isempty(rel_wave)
            lit_wave = waves(rel_wave);
            paper_size = [1000,400];
            do_lit_wave = true;
        else
            paper_size = [1000,400];
            do_lit_wave = false;
        end
        
        % Make figure
        figure('Position', [20,20, paper_size])
        subplot('Position', [0.55,0.20,0.42,0.68])
        
        leg_labels = {}; max_t = 0;
        
        eval(['units = data.waves.units.' curr_wave_type ';'])
        
        if plot_model_data_normalised
            units = 'normalised';
        end
        
        % Plot each of the waves in turn (corresponding to the three ages)
        leg_h = [];
        for age_no = 1 : length(rel_ages)
            
            % extract data to plot
            eval(['curr_wav.v = data.waves.' curr_wave_type '_', curr_site '{rel_sims(age_no)};']);
            curr_wav.fs = data.waves.fs;
            curr_wav.t = [0:length(curr_wav.v)-1]/curr_wav.fs;
                        
            % plot
            if ~plot_model_data_normalised
                leg_h(end+1) = plot(curr_wav.t, curr_wav.v, 'Color', age_colors(age_no)*ones(1,3), 'LineWidth', lwidth); hold on
            else
                curr_wav.v = (curr_wav.v - min(curr_wav.v))/range(curr_wav.v);
                leg_h(end+1) = plot(curr_wav.t, curr_wav.v +(offset*(age_no-1)), 'Color', age_colors(age_no)*ones(1,3), 'LineWidth', lwidth); hold on
                ylim([-0.1 1.1+(offset*(length(rel_ages)-1))]);
                set(gca, 'YTick', [])
            end
            
            % store legend label
            leg_labels = [leg_labels, num2str(rel_ages(age_no))];
            max_t = max([max_t, max(curr_wav.t)]);
            
            % save baseline fig
            if age_no == 1
                
                % tidy up
                xlim([0, max_t])
                xlabel('Time [s]', 'FontSize', ftsize)
                ylab = ylabel([curr_wave_type, '  [' units ']'], 'FontSize', ftsize);
                set(gca, 'FontSize', ftsize)
                box off
                title([strrep(artery_name, 'cR', 'c R'), ', ' curr_wave_type], 'FontSize', ftsize+4)
                
                % Save fig
                PrintFigs(gcf, paper_size/70, [PATHS.Analysis_figures, curr_wave_type, '_waves_', artery_name, '_baseline'], false)
                
                title('')
            end
            
        end
        
        % tidy up
        xlim([0, max_t])
        xlabel('Time [s]', 'FontSize', ftsize)
        ylab = ylabel([curr_wave_type, '  [' units ']'], 'FontSize', ftsize);
        %ylab = ylabel({curr_wave_type, ['[' units ']']}, 'FontSize', ftsize, 'Rotation', 0);
        %set(ylab, 'Units', 'Normalized', 'Position', [-0.17, 0.5, 0]);
        set(gca, 'FontSize', ftsize)
        legend(fliplr(leg_h), fliplr(leg_labels)), clear leg_labels
        str = [strrep(artery_name, 'cR', 'c R'), ', ' curr_wave_type];
        str = strrep(str, 'SupMidCerebral', 'Mid Cerebral');
        dim = [0.25,0.7,0.5,0.3];
        annotation('textbox',dim,'String',str,'LineStyle', 'none', 'FontSize', ftsize+8, 'HorizontalAlignment', 'center');
        box off
        
        % Plot waves from literature
        subplot('Position', [0.05,0.20,0.42,0.68])
        if ~do_lit_wave
           lit_wave.v(1).y = nan(1,100);
           lit_wave.v(2).y = nan(1,100);
           lit_wave.v(3).y = nan(1,100);
           lit_wave.v(1).fs = 100;
           lit_wave.v(2).fs = 100;
           lit_wave.v(3).fs = 100;
           lit_wave.age = [25,55,75];
        end
        
        leg_labels = {}; max_t = 0;
        if length(lit_wave.v) < 3
            temp_rel_ages = rel_ages([1,3]);
        elseif strcmp(curr_wave_type, 'PPG')
            temp_rel_ages = [25,35,55];
        elseif strcmp(artery_name, 'Brachial')
            temp_rel_ages = [25,35,75];
        else
            temp_rel_ages = rel_ages;
        end
        
        % Plot each of the waves in turn (corresponding to the three ages)
        leg_h = [];
        for age_no = 1 : length(temp_rel_ages)
            
            % identify relevant wave to plot
            [~, rel_wave_el] = min(abs(lit_wave.age - temp_rel_ages(age_no)));
            curr_wav.v = lit_wave.v(rel_wave_el).y;
            curr_wav.fs = lit_wave.v(rel_wave_el).fs;
            curr_wav.t = [0:length(curr_wav.v)-1]/curr_wav.fs;
            curr_wav.age = lit_wave.age(rel_wave_el);
            
            % plot
            leg_h(end+1) = plot(curr_wav.t, curr_wav.v + (offset*(age_no-1)), 'Color', age_colors(age_no)*ones(1,3), 'LineWidth', lwidth); hold on
            ylim([-0.1 1.1+(offset*(length(temp_rel_ages)-1))]);
            
            % store legend label
            leg_labels = [leg_labels, num2str(curr_wav.age,2)];
            max_t = max([max_t, max(curr_wav.t)]);
            
        end
        
        % tidy up
        xlim([0, max_t])
        if strcmp(artery_name, 'SupMidCerebral') || strcmp(artery_name, 'Digital')
            xlabel('Time [normalised]', 'FontSize', ftsize)
        else
            xlabel('Time [s]', 'FontSize', ftsize)
        end
        ylab = ylabel([curr_wave_type, '  [normalised]'], 'FontSize', ftsize);
        %ylab = ylabel({curr_wave_type, ['[' units ']']}, 'FontSize', ftsize, 'Rotation', 0);
        %set(ylab, 'Units', 'Normalized', 'Position', [-0.17, 0.5, 0]);
        set(gca, 'FontSize', ftsize, 'YTick', [])
        legend(fliplr(leg_h), fliplr(leg_labels)), clear leg_labels
%         title([strrep(artery_name, 'cR', 'c R'), ' (\itin vivo\rm\bf)'], 'FontSize', ftsize+4)
        box off
        
        % save
        PrintFigs(gcf, paper_size/70, [PATHS.Analysis_figures, curr_wave_type, '_waves_', strrep(artery_name, ' ', '_')])
        
    end


end

end

function artery_name = find_artery_name(site)

switch site
    case 'Brachial'
        artery_name = 'Arm';
    case 'Radial'
        artery_name = 'Wrist';
    case 'Carotid'
        artery_name = 'Neck';
    case 'Digital'
        artery_name = 'Finger';
    case 'Femoral'
        artery_name = 'UpperLeg';
    case {'Anterior Tibial', 'AntTibial'}
        artery_name = 'Ankle';
    case 'Thyroid'
        artery_name = 'Ear';
    case 'SupTemporal'
            artery_name = 'Ear';
end

end

function PrintFigs(h, paper_size, savepath, close_fig)

if nargin<4
    close_fig = true;
end

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
if close_fig
    close all;
end

% save 
fid = fopen([savepath, '.txt'], 'w');
p = mfilename('fullpath');
p = strrep(p, '\', '\\');
fprintf(fid, ['Figures generated by:\n\n ' p '.m \n\n on ' datestr(today)]);
fclose all;

end
