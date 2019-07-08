function asi_case_study(pwdb_no)
% ASI_CASE_STUDY generates the plots reported in the case study on
% assessing aortic stiffness from digital PPG waves
%
%               asi_case_study
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
% v.1.0    Peter H. Charlton, King's College London

fprintf('\n --- Running Aortic Stiffness Indices (ASI) Case Study ---')

% Setup paths with current simulation paths
PATHS = setup_paths_for_post_processing(pwdb_no);

% Create ASI plots
create_plots = 0;
if create_plots, create_asi_plots(PATHS, pwdb_no); end

% Extract ASIs
collate_ASIs(PATHS);

% % Quantify correlations of all ASIs
analyse_ASIs(PATHS, pwdb_no);

% Make correlation graphs
correlation_graphs(PATHS, pwdb_no);

% Make results figures
% results_figures(PATHS, pwdb_no);

end

function create_asi_plots(PATHS, sim_no)

fprintf('\n - Creating Arterial Stiffness Index Plots')
%% Pulse wave indices

% Load collated data
load(PATHS.collated_data)

% Find baseline simulation
[data.baseline_sim_all, data.baseline_sim_age] = deal(false(length(collated_data),1));
for sim_no = 1 : length(collated_data)
    if sum(collated_data(sim_no).input_data.sim_settings.variations.params ~= 0) == 0
        data.baseline_sim_age(sim_no) = true;
        if collated_data(sim_no).input_data.sim_settings.age == 25
            data.baseline_sim_all(sim_no) = true;
        end
    end
end
sim_to_use = find(data.baseline_sim_all);

% Make plots in turn
sites = {'digital', 'carotid', 'carotid', 'radial'};
%site_domain_no = [15, 15, 22, 112];
signals = {'PPG', 'PPG', 'P','P', 'PPG'};
domains = extractfield(collated_data(1).output_data, 'domain_no');
fs = 1000;

for plot_no = 1 : length(sites)
    
    % Identify relevant data
    switch sites{plot_no}
        case 'carotid'
            dom_no = 15;
            dist_prop = 0.5;
        case 'radial'
            dom_no = 22;
            dist_prop = 1;
        case 'digital'
            dom_no = 112;
            dist_prop = 1;
    end
    rel_row = find(domains == dom_no);
    rel_data = collated_data(sim_to_use).output_data(rel_row);
    req_dist = dist_prop*collated_data(sim_to_use).input_data.sim_settings.network_spec.length(dom_no);
    [~, rel_dist_el] = min(abs(rel_data.distances-req_dist));
    rel_data.P = rel_data.P(:,rel_dist_el); 
    rel_data.PPG = rel_data.PPG_WK(:,rel_dist_el);
    rel_data.t = [0:length(rel_data.P)-1]/fs;
    
    sig.fs = fs;
    eval(['sig.v = rel_data.' signals{plot_no} ';']);
    options.save_folder = '/Users/petercharlton/Google Drive/Work/Publications/In Preparation/2018 BIHS Conference/figure/';
    options.save_file = [sites{plot_no}, '_', signals{plot_no}, '_'];
    options.plot_third_deriv = false;
    PulseAnalyse6(sig, options);
    
end


end

function collate_ASIs(PATHS)

fprintf('\n - Collating Arterial Stiffness Indices and Input Parameters')
%% Pulse wave indices

% Load collated data
load(PATHS.exported_data_mat_pwdb_data)

% settings
pw_param_names = {'AGI_mod', 'AI_c', 'RI', 'SI'};
pwv_names = {'a', 'cf', 'br', 'fa'};
    
% Cycle through each simulation
subj_counter = 0;
for subj_no = 1 : length(data.config.age)
    
    % skip if this is an implausible subject
    if ~data.plausibility.plausibility_log(subj_no)
        continue
    end
    subj_counter = subj_counter+1;
    
    % record sim no
    asi_data.sim(subj_counter,1) = subj_no;
    asi_data.baseline_sim_all(subj_counter,1) = data.config.baseline_sim_for_all(subj_no);
    asi_data.baseline_sim_age(subj_counter,1) = data.config.baseline_sim_for_age(subj_no);
    
    % Extract model input parameters for this simulation
    input_vars = {'hr', 'sv', 'lvet', 'mbp', 'dia', 'pwv', 'age'};
    for input_var_no = 1 : length(input_vars)
        eval(['asi_data.i_' input_vars{input_var_no} '(subj_counter,1) = data.config.' input_vars{input_var_no} '(subj_no);']);
        
        if ~strcmp(input_vars{input_var_no}, 'age')
            rel_col = find(strcmp(input_vars{input_var_no}, data.config.variations.param_names));
            eval(['asi_data.i_' input_vars{input_var_no} '_SD(subj_counter,1) = data.config.variations.params(subj_no, rel_col);']);
        end
        
    end
    
    feat_no = 0;
        
    % Cycle through pulse wave features
    for param_no = 1 : length(pw_param_names)
        curr_param_name = pw_param_names{param_no};
        feat_no = feat_no + 1;
        
        % Extract pulse wave features (at each site, for each signal)
        asi_data.feat_names{1,feat_no} = curr_param_name;
        
        eval(['asi_data.val(subj_counter,feat_no) = data.haemods(subj_no).' curr_param_name ';']);
        
        clear curr_param_name
    end
    clear param_no
    
    
    % Add in Pulse wave velocities
    for pwv_name_no = 1 : length(pwv_names)
        curr_param_name = pwv_names{pwv_name_no};
        eval(['asi_data.pwv_' curr_param_name '(subj_counter,1) = data.haemods(subj_no).PWV_' curr_param_name ';']);
    end
    clear pwv_name_no
    
end
clear subj_no subj_counter

save(PATHS.collated_ASIs, 'asi_data')

end

function correlation_graphs(PATHS, sim_no)

fprintf('\n - Making correlation graphs')

% setup
% rel_min_age = 25; rel_max_age = 55;
rel_min_age = 45; rel_max_age = 45;
lwidth = 2;
ftsize = 16;
markersize = 7;
paper_size = [400,300];

% Load ASI data
load(PATHS.collated_ASIs)
load(PATHS.ASI_results)

% Highlighted age group
for feat_no = 1 : length(asi_data.feat_names)
    
    figure('Position', [20,20,paper_size])
    plot(asi_data.pwv_a, asi_data.val(:,feat_no), '.', 'MarkerSize', 8, 'Color', 0.2*[1 1 1]), hold on
    h = lsline;
    h.Color = 'k';
    
    rel_els = asi_data.i_age >= rel_min_age & asi_data.i_age <= rel_max_age;
    plot(asi_data.pwv_a(rel_els), asi_data.val(rel_els,feat_no), '.r', 'MarkerSize', 8), hold on
    
    set(gca, 'FontSize', ftsize)
    xlabel('Aortic PWV (m/s)', 'FontSize', ftsize)
    ylab_txt = strrep(asi_data.feat_names{feat_no}, '_', ' ');
    ylab_txt = strrep(ylab_txt, 'AI c', 'AIx');
    ylab_txt = strrep(ylab_txt, 'AGI mod', 'AGI mod (unitless)');
    ylab_txt = strrep(ylab_txt, 'SI', 'SI (m/s)');
    ylab_txt = strrep(ylab_txt, 'RI', 'RI (unitless)');
    ylabel(ylab_txt, 'FontSize', ftsize)
    feat_no = find(strcmp(rsq.feat_names, asi_data.feat_names{feat_no}));
    rel_rsq = rsq.v(feat_no);
    if rel_min_age == rel_max_age
        rel_rsq_age = rsq_age.v(feat_no,find(rsq_age.age >= rel_min_age & rsq_age.age <= rel_max_age));
    else
        refs = asi_data.pwv_a;
        vals = asi_data.val(:,feat_no);
        rel_els = find(asi_data.i_age >= rel_min_age & asi_data.i_age <= rel_max_age & ~isnan(vals) );
        temp = corrcoef(refs(rel_els), vals(rel_els));
        rel_rsq_age = temp(1,2)^2; clear temp
        
    end
    
    % R2 for all data
    dim1 = [.2 .65 .3 .3];
    dim2 = [.2 .55 .3 .3];
    str = ['R^2 = ' num2str(rel_rsq,'%0.2f')];
    annotation('textbox',dim1,'String',str,'FitBoxToText','on','LineStyle','None', 'FontSize', ftsize, 'Color', 0.2*[1,1,1]);
    str = ['R^2 = ' num2str(rel_rsq_age,'%0.2f')];
    annotation('textbox',dim2,'String',str,'FitBoxToText','on','LineStyle','None', 'FontSize', ftsize, 'Color', 'r');
    
    box off
    curr_range = range(asi_data.val(:,feat_no));
    ylim([min(asi_data.val(:,feat_no))-0.1*curr_range, max(asi_data.val(:,feat_no))-0.1*curr_range])
    if strcmp(ylab_txt, 'AIx')
        ylim([min(asi_data.val(:,feat_no))-0.1*curr_range, 35])
    end
    curr_range = range(asi_data.pwv_a);
    xlim([min(asi_data.pwv_a)-0.1*curr_range, max(asi_data.pwv_a)-0.1*curr_range])
    PrintFigs(gcf, paper_size/70, [PATHS.CaseStudies, 'ASI_correlation_plot_hl_' asi_data.feat_names{feat_no}])
    
end

% All data
for feat_no = 1 : length(asi_data.feat_names)
    
    figure('Position', [20,20,paper_size])
    plot(asi_data.pwv_a, asi_data.val(:,feat_no), 'x'), hold on
    h = lsline;
    h.Color = 'k';
    
    set(gca, 'FontSize', ftsize)
    xlabel('Aortic PWV [m/s]', 'FontSize', ftsize)
    ylabel('Modified Ageing Index [au]', 'FontSize', ftsize)
    ylabel(strrep(asi_data.feat_names{feat_no}, '_', ' '), 'FontSize', ftsize)
    feat_no = find(strcmp(rsq.feat_names, asi_data.feat_names{feat_no}));
    rel_rsq = rsq.v(feat_no);
    title(['R^2 = ' num2str(rel_rsq,'%0.2f')], 'FontSize', ftsize)
    box off
    curr_range = range(asi_data.val(:,feat_no));
    ylim([min(asi_data.val(:,feat_no))-0.1*curr_range, max(asi_data.val(:,feat_no))-0.1*curr_range])
    curr_range = range(asi_data.pwv_a);
    xlim([min(asi_data.pwv_a)-0.1*curr_range, max(asi_data.pwv_a)-0.1*curr_range])
    PrintFigs(gcf, paper_size/70, [PATHS.CaseStudies, 'ASI_correlation_plot_' asi_data.feat_names{feat_no}])
    
end

end

function analyse_ASIs(PATHS, sim_no)

fprintf('\n - Analysing Performances of Arterial Stiffness Indices')

paper_size = [200,200,700,300];

% Load ASI data
load(PATHS.collated_ASIs)

% calculate correlation coefficients
refs = asi_data.pwv_a;
for feat_no = 1 : length(asi_data.feat_names)
    vals = asi_data.val(:,feat_no);
    rel_els = find(~isnan(vals));
    [temp, temp2] = corrcoef(refs(rel_els), vals(rel_els));
    rsq.v(feat_no,1) = temp(1,2)^2;
    p(feat_no,1) = temp2(2); clear temp temp2
    
    % best-fit line
    f = fit(vals(rel_els), refs(rel_els), 'poly2');
    est_refs = feval(f, vals(rel_els));
    
    % BA stats
    errors = est_refs - refs(rel_els);
    fit_asi_ba.bias(feat_no) = mean(errors);
    fit_asi_ba.sd(feat_no) = std(errors);
    fit_asi_ba.mae(feat_no) = mean(abs(errors));
    
    fit_asi_ba.med_ae(feat_no) = median(abs(errors));
    fit_asi_ba.lq_ae(feat_no) = quantile(abs(errors), 0.25);
    fit_asi_ba.uq_ae(feat_no) = quantile(abs(errors), 0.75);
    
    %clear errors f est_refs
    
end
rsq.feat_names = asi_data.feat_names;
fit_asi_ba.names = asi_data.feat_names(:);

% calculate correlation coefficients for each age group
refs = asi_data.pwv_a;
ages = unique(asi_data.i_age);
for feat_no = 1 : length(asi_data.feat_names)
    for age_no = 1 : length(ages)
        vals = asi_data.val(:,feat_no);
        rel_els = find(~isnan(vals) & asi_data.i_age == ages(age_no));
        temp = corrcoef(refs(rel_els), vals(rel_els));
        rsq_age.v(feat_no,age_no) = temp(1,2)^2; clear temp
    end
end
rsq_age.feat_names = asi_data.feat_names;
rsq_age.age = ages;

% calculate pwv accuracy and precision
temp = fieldnames(asi_data);
pwv_types = temp(~cellfun(@isempty, strfind(temp, 'pwv_')));
pwv_types = pwv_types(~strcmp(pwv_types, 'i_pwv_SD'));
pwv_types = pwv_types(~strcmp(pwv_types, 'pwv_aortic'));
for pwv_no = 1 : length(pwv_types)
    
    % Rsq stats
    eval(['vals = asi_data.' pwv_types{pwv_no} ';']);
    rel_els = ~isnan(vals);
    temp = corrcoef(refs(rel_els), vals(rel_els));
    ba.rsq(pwv_no,1) = temp(1,2)^2; clear temp
    
    % BA stats
    errors = vals(rel_els) - refs(rel_els);
    ba.bias(pwv_no) = mean(errors);
    ba.sd(pwv_no) = std(errors);
    
    % best-fit line
    f = fit(vals(rel_els), refs(rel_els), 'poly2');
    est_refs = feval(f, vals(rel_els));
    
    % BA stats
    errors = est_refs - refs(rel_els);
    fit_ba.bias(pwv_no) = mean(errors);
    fit_ba.sd(pwv_no) = std(errors);
    fit_ba.mae(pwv_no) = mean(abs(errors));
    
    fit_ba.med_ae(pwv_no) = median(abs(errors));
    fit_ba.lq_ae(pwv_no) = quantile(abs(errors), 0.25);
    fit_ba.uq_ae(pwv_no) = quantile(abs(errors), 0.75);
    clear errors f est_refs
    
    % Age-specific Rsq values
    for age_no = 1 : length(ages)
        rel_els = find(~isnan(vals) & asi_data.i_age == ages(age_no));
        temp = corrcoef(refs(rel_els), vals(rel_els));
        ba.rsq_age.v(pwv_no,age_no) = temp(1,2)^2; clear temp
    end
    clear age_no rel_els vals errors
end
ba.feat_names = pwv_types;
fit_ba.names = pwv_types(:);

% De-compose names of pulse wave inds
rsq.ind{feat_no} = rsq.feat_names{feat_no};

save(PATHS.ASI_results, 'rsq', 'rsq_age', 'ba', 'fit_ba', 'fit_asi_ba');

end

function results_figures(PATHS, sim_no)

fprintf('\n - Making results figures')

load(PATHS.ASI_results)

% setup
lwidth = 2;
ftsize = 16;
markersize = 7;

%% Correlation results summary
do_plot = 1;
if do_plot
    
    addpath(genpath('/Users/petercharlton/Google Drive/Work/Code/Tools'))
    paper_size = [450,270];
    figure('Position', [20,20,paper_size])
    subplot('Position', [0.03,0.19,0.95,0.81])
    ftsize = 16;
    corr_res.pwv = ba.rsq(:);
    corr_res.asi = rsq.v(:);
    
    data = [corr_res.pwv;corr_res.asi];
    catIdx = [ones(size(corr_res.pwv));zeros(size(corr_res.asi))];
    plotSpread(data,'categoryIdx',catIdx,...
        'categoryMarkers',{'+','o'},'categoryColors',{'r','b'})
    view([90 -90])
    lgd = legend({'Single-PW ASIs','PWVs'},'Location','SouthEast');
    lgd.Orientation = 'Horizontal';
    ylabel('Correlation with aortic PWV, R^2', 'FontSize', ftsize)
    set(gca, 'FontSize', ftsize, 'XTick', [])
    xlim([0.25 1.6])
    ylim([0 1])

    PrintFigs(gcf, paper_size/70, [PATHS.CaseStudies, 'ASI_analysis_corr_summary'])
    clear corr_res
end

%% Correlation all vs age-groups
do_plot = 1;
if do_plot
    
    paper_size = [540,250];
    figure('Position', [20,20,paper_size])
    ftsize = 20;
    temp_els = [find(strcmp(rsq.feat_names, 'PPG_15_c_div_a')), ...
        find(strcmp(rsq.feat_names, 'P_22_b_div_a')), ...
        find(strcmp(rsq.feat_names, 'PPG_112_d_div_a'))];
    corr_res.pwv.all = ba.rsq(:);
    corr_res.asi.all = rsq.v(temp_els);
    corr_res.pwv.age = ba.rsq_age.v(:);
    corr_res.asi.age = rsq_age.v(temp_els,:); corr_res.asi.age = corr_res.asi.age(:);
    
    data = [corr_res.pwv.age; corr_res.pwv.all];
    gps = [ones(size(corr_res.pwv.age)); 2*ones(size(corr_res.pwv.all))];
    boxplot(data, gps, 'Widths', 0.7)
    box off
    ylim([0 1])
    view([90 -90])
    set(gca,'FontSize', ftsize, 'XTickLabel', {'Age-specific','All data'})
    ylabel('Correlation with Aortic PWV, R^2', 'FontSize', ftsize)
    abc = findobj(gca, 'Type', 'Line');
    set(abc, 'LineWidth', 2)
    PrintFigs(gcf, paper_size/70, [PATHS.CaseStudies, 'ASI_analysis_corr_age_all_PWV'])
    
    figure('Position', [20,20,paper_size])
    data = [corr_res.asi.age; corr_res.asi.all];
    gps = [ones(size(corr_res.asi.age)); 2*ones(size(corr_res.asi.all))];
    boxplot(data, gps, 'Widths', 0.7)
    box off
    ylim([0 1])
    view([90 -90])
    set(gca,'FontSize', ftsize, 'XTickLabel', {'Age-specific','All data'})
    ylabel('Correlation with Aortic PWV, R^2', 'FontSize', ftsize)
    abc = findobj(gca, 'Type', 'Line');
    set(abc, 'LineWidth', 2)
    PrintFigs(gcf, paper_size/70, [PATHS.CaseStudies, 'ASI_analysis_corr_age_all_ASI'])
    
end

%% PWV B-A results
do_plot = 1;
if do_plot
    
    paper_size = [400,270];
    figure('Position', [20,20,paper_size])
    subplot('Position', [0.29,0.17,0.69,0.81])
    pwv_names = ba.feat_names;
    
    [~,order] = sort(ba.bias);
    
    counter = 0;
    for pwv_no = order
        
        % Plot the results for this PWV technique
        counter = counter+1;
        y_val = counter;
        line_coords.x_loa = [ba.bias(pwv_no) - 2*ba.sd(pwv_no), ba.bias(pwv_no) + 2*ba.sd(pwv_no)];
        line_coords.x_bias = ba.bias(pwv_no);
        line_coords.y = y_val*ones(1,2);
        % - LOAs
        plot(line_coords.x_loa, line_coords.y, 'o-k', 'LineWidth', lwidth, 'MarkerFaceColor', 'k', 'MarkerSize', markersize), hold on
        % - Bias
        plot(line_coords.x_bias, line_coords.y(1), 'dk', 'LineWidth', lwidth, 'MarkerFaceColor', 'k', 'MarkerSize', markersize+2), hold on
        
    end
    
    % Tidy-up
    ax = gca;
    set(gca, 'FontSize', ftsize, 'YTick', 1:length(pwv_names), 'YTickLabel', strrep(strrep(pwv_names(order), 'pwv_', ''), '_', '-'))
    xlabel('Error [m/s]', 'FontSize', ftsize)
    ylim([0.5 length(pwv_names)+0.5])
    box off
    legend('LOAs', 'Bias','Location','SouthEast')
    PrintFigs(gcf, paper_size/70, [PATHS.CaseStudies, 'ASI_analysis_PWV_BA'])
    
end

%% Plot Median AEs

do_plot = 0;
if do_plot
    [~,order] = sort(fit_ba.med_ae);
    order = order(1:3);
    
    counter = 0;
    for pwv_no = order
        
        % Plot the results for this PWV technique
        counter = counter+1;
        y_val = counter;
        line_coords.x_loa = [fit_ba.lq_ae(pwv_no), fit_ba.uq_ae(pwv_no)];
        line_coords.x_bias = fit_ba.med_ae(pwv_no);
        line_coords.y = y_val*ones(1,2);
        % - LOAs
        plot(line_coords.x_loa, line_coords.y, 'o-k', 'LineWidth', lwidth, 'MarkerFaceColor', 'k', 'MarkerSize', markersize), hold on
        % - Bias
        plot(line_coords.x_bias, line_coords.y(1), 'dk', 'LineWidth', lwidth, 'MarkerFaceColor', 'k', 'MarkerSize', markersize+2), hold on
        
    end
    
    temp_els = [find(strcmp(fit_asi_ba.names, 'AGI_mod')), ...
        find(strcmp(fit_asi_ba.names, 'P_22_b_div_a')), ...
        find(strcmp(fit_asi_ba.names, 'PPG_112_slope_b_d'))];
    [~,temp] = sort(fit_asi_ba.mae(temp_els)); rel_asi_ba_els = temp_els(temp(1:3));
    clear temp
    for asi_no = 1 : length(rel_asi_ba_els)
        
        curr_el = rel_asi_ba_els(asi_no);
        
        % Plot the results for this ASI technique
        counter = counter+1;
        y_val = counter;
        line_coords.x_loa = [fit_asi_ba.lq_ae(curr_el), fit_asi_ba.uq_ae(curr_el)];
        line_coords.x_bias = fit_asi_ba.med_ae(curr_el);
        line_coords.y = y_val*ones(1,2);
        % - LOAs
        plot(line_coords.x_loa, line_coords.y, 'o-b', 'LineWidth', lwidth, 'MarkerFaceColor', 'b', 'MarkerSize', markersize), hold on
        % - Bias
        plot(line_coords.x_bias, line_coords.y(1), 'db', 'LineWidth', lwidth, 'MarkerFaceColor', 'b', 'MarkerSize', markersize+2), hold on
        
        
        
    end
    
    % Tidy-up
    ax = gca;
    all_names = [fit_ba.names(order); fit_asi_ba.names(rel_asi_ba_els)];
    set(gca, 'FontSize', ftsize, 'YTick', 1:length(all_names), 'YTickLabel', strrep(strrep(all_names, 'pwv_', ''), '_', '-'))
    xlabel('Error [m/s]', 'FontSize', ftsize)
    ylim([0.5 length(all_names)+0.5])
    box off
    legend('IQR', 'Median','Location','SouthEast')
    PrintFigs(gcf, paper_size/70, [PATHS.CaseStudies, 'ASI_analysis_PWV_med_AE_fit'])
end

%% Plot fitted B-A Results

do_plot = 0;
if do_plot
    
    [~,order] = sort(fit_ba.sd);
    
    counter = 0;
    for pwv_no = order
        
        % Plot the results for this PWV technique
        counter = counter+1;
        y_val = counter;
        line_coords.x_loa = [fit_ba.bias(pwv_no) - 2*fit_ba.sd(pwv_no), fit_ba.bias(pwv_no) + 2*fit_ba.sd(pwv_no)];
        line_coords.x_bias = fit_ba.bias(pwv_no);
        line_coords.y = y_val*ones(1,2);
        % - LOAs
        plot(line_coords.x_loa, line_coords.y, 'o-k', 'LineWidth', lwidth, 'MarkerFaceColor', 'k', 'MarkerSize', markersize), hold on
        % - Bias
        plot(line_coords.x_bias, line_coords.y(1), 'dk', 'LineWidth', lwidth, 'MarkerFaceColor', 'k', 'MarkerSize', markersize+2), hold on
        
    end
    
    [~,temp] = sort(fit_asi_ba.sd); rel_asi_ba_els = temp(1:5); clear temp
    for asi_no = 1 : length(rel_asi_ba_els)
        
        curr_el = rel_asi_ba_els(asi_no);
        
        % Plot the results for this ASI technique
        counter = counter+1;
        y_val = counter;
        line_coords.x_loa = [fit_asi_ba.bias(curr_el) - 2*fit_asi_ba.sd(curr_el), fit_asi_ba.bias(curr_el) + 2*fit_asi_ba.sd(curr_el)];
        line_coords.x_bias = fit_asi_ba.bias(curr_el);
        line_coords.y = y_val*ones(1,2);
        % - LOAs
        plot(line_coords.x_loa, line_coords.y, 'o-b', 'LineWidth', lwidth, 'MarkerFaceColor', 'b', 'MarkerSize', markersize), hold on
        % - Bias
        plot(line_coords.x_bias, line_coords.y(1), 'db', 'LineWidth', lwidth, 'MarkerFaceColor', 'b', 'MarkerSize', markersize+2), hold on
        
        
    end
    
    % Tidy-up
    ax = gca;
    set(gca, 'FontSize', ftsize, 'YTick', 1:length(pwv_names), 'YTickLabel', strrep(strrep(pwv_names(order), 'pwv_', ''), '_', '-'))
    xlabel('Error [m/s]', 'FontSize', ftsize)
    ylim([0.5 length(pwv_names)+0.5])
    box off
    legend('95% CI', 'Mean','Location','SouthEast')
    PrintFigs(gcf, paper_size/70, [PATHS.CaseStudies, 'ASI_analysis_PWV_BA_fit'])
    
end

%% Correlation plot for all ages and only one age group

load(PATHS.collated_ASIs)
rel_asi_name = 'AGI_mod';
%rel_asi_name = 'best';
if strcmp(rel_asi_name, 'best')
    [~, order] = sort(rsq.v, 'descend');
    rel_el = order(1);
else
    rel_el = find(strcmp(rsq.feat_names, rel_asi_name));
end

% make plot
figure('Position', [20,20,paper_size])
subplot('Position', [0.20,0.17,0.78,0.80])
plot(asi_data.pwv_a, asi_data.val(:,rel_el),'xb','LineWidth',2), hold on,
lsline
set(gca, 'FontSize', ftsize)
xlabel('Reference aortic PWV [m/s]', 'FontSize', ftsize)
if rel_el == 1
    ylabel(strrep(rel_asi_name, '_', ' '), 'FontSize', ftsize)
else
    ylab = ylabel('ASI', 'FontSize', ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.15, 0.5, 0]);
end

box off
%text(10,-0.9,{'All ages', ['R^2 = ' num2str(ranking(1),3)]}, 'Color', 'b','FontSize', ftsize)
text(10,-0.7,{'All ages', ['R^2 = ' num2str(rsq.v(rel_el),2)]}, 'Color', 'b','FontSize', ftsize)
PrintFigs(gcf, paper_size/70, [PATHS.CaseStudies, 'ASI_analysis_correlation_all_ages'], 0)
%text(6,-0.96,{'Aged 40-50', ['R^2 = ' num2str(rsq_age.v(order(1)),3)]}, 'Color', 'r','FontSize', ftsize)
rel_els = asi_data.i_age == 45;
plot(asi_data.pwv_a(rel_els), asi_data.val(rel_els,rel_el),'xr','LineWidth',2)
text(6,0,{'Aged 40-50', ['R^2 = ' num2str(rsq_age.v(rel_el, 3),1)]}, 'Color', 'r','FontSize', ftsize)
PrintFigs(gcf, paper_size/70, [PATHS.CaseStudies, 'ASI_analysis_correlation_single_age'])


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