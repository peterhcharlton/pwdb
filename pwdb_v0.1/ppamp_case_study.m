function ppamp_case_study(pwdb_no)
% PPAMP_CASE_STUDY generates the plots reported in the case study on the
% determinants of pulse pressure amplification.
%
%               ppamp_case_study
%
%   Inputs:     - the 'pwdb_data.mat' file produced by 'export_pwdb.m'.
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
% v.1.0    Peter H Charlton, King's College London

fprintf('\n --- Running Pulse Pressure Amplification (PP_amp) Case Study ---')

% Setup paths with current simulation paths
PATHS = setup_paths_for_post_processing(pwdb_no);

% Collate data
data = load_data(PATHS);

% Investigate
make_plots(data, PATHS);

end

function data = load_data(PATHS)

fprintf('\n - Loading data')

% Load all data
load(PATHS.exported_data_mat_pwdb_data);

end

function make_plots(data, PATHS)

fprintf('\n - Making Plots')

rel_fields = {'PWV_a', 'age', 'dia_asc_a', 'P1in_a', 'P2in_a', 'SV', 'HR', 'LVET', 'CO', 'PP_amp', 'PWV_cf', 'PWV_br', 'DBP_a', 'P1pk_a', 'P2pk_a', 'PP_a', 'AP_a', 'DBP_b', 'P1pk_b', 'P2pk_b', 'P1in_b', 'P2in_b', 'PP_b', 'PP_f', 'SBP_a', 'SBP_b', 'MBP_a', 'MBP_b', 'dia_asc_a', 'SMBP_a', 'SBP_diff', 'c', 'AI_a', 'AI_c', 'P1in_c'};
for field_no = 1 : length(rel_fields)
    
    curr_field = rel_fields{field_no};
    eval(['d.' curr_field ' = extractfield(data.haemods, ''' curr_field ''');']);
    
end

%% Make plot of PPamp vs age (using P1 and P2)

ages = unique(d.age);
for age_no = 1 : length(ages)
    curr_age = ages(age_no);
    rel_els = d.age == curr_age & data.plausibility.plausibility_log(:)';
    sbpb.v(age_no) = mean(d.SBP_b(rel_els));
    sbpb.sd(age_no) = std(d.SBP_b(rel_els));
    p1a.v(age_no) = mean(d.P1in_a(rel_els));
    p1a.sd(age_no) = std(d.P1in_a(rel_els));
    p1b.v(age_no) = mean(d.P1in_b(rel_els));
    p1b.sd(age_no) = std(d.P1in_b(rel_els));
    p2a.v(age_no) = mean(d.P2pk_a(rel_els));
    p2a.sd(age_no) = std(d.P2pk_a(rel_els));
    p2b.v(age_no) = mean(d.P2pk_b(rel_els));
    p2b.sd(age_no) = std(d.P2pk_b(rel_els));
    dbpb.v(age_no) = mean(d.DBP_b(rel_els));
    dbpb.sd(age_no) = std(d.DBP_b(rel_els));
    sbpa.v(age_no) = mean(d.SBP_a(rel_els));
    sbpa.sd(age_no) = std(d.SBP_a(rel_els));
    dbpa.v(age_no) = mean(d.DBP_a(rel_els));
    dbpa.sd(age_no) = std(d.DBP_a(rel_els));
    apa.v(age_no) = mean(d.AP_a(rel_els));
    apa.sd(age_no) = std(d.AP_a(rel_els));
    ppamp.v(age_no) = mean( (d.SBP_b(rel_els)-d.DBP_b(rel_els))./(d.SBP_a(rel_els)-d.DBP_a(rel_els)) );
    ppamp.sd(age_no) = std( (d.SBP_b(rel_els)-d.DBP_b(rel_els))./(d.SBP_a(rel_els)-d.DBP_a(rel_els)) );
    ppamp_p1.v(age_no) = mean( (d.SBP_b(rel_els)-d.DBP_b(rel_els))./(d.P1in_a(rel_els)-d.DBP_a(rel_els)) );
    ppamp_p1.sd(age_no) = std( (d.SBP_b(rel_els)-d.DBP_b(rel_els))./(d.P1in_a(rel_els)-d.DBP_a(rel_els)) );
    ppamp_p2.v(age_no) = mean( (d.SBP_b(rel_els)-d.DBP_b(rel_els))./(d.P2pk_a(rel_els)-d.DBP_a(rel_els)) );
    ppamp_p2.sd(age_no) = std( (d.SBP_b(rel_els)-d.DBP_b(rel_els))./(d.P2pk_a(rel_els)-d.DBP_a(rel_els)) );
end

ftsize = 16; lwidth = 2;
age_ticks = unique(data.config.age);
ppamp_types = {'ppamp_p1', 'ppamp_p2', 'ppamp'};
req_colors = {[1 0 0], [0 0 1], [0 0 0]};
no_sd = 1; up.ylim_offset = 1.1; jitter = 0.6;
for type_no = 1 : length(ppamp_types)

    req_color= req_colors{type_no};
    eval(['rel_data = ' ppamp_types{type_no} ';']);
    rel_data.age = ages;
    
    errorbar(rel_data.age-(2-type_no)*jitter, rel_data.v, rel_data.sd, '-o', 'Color', req_color, 'LineWidth', 2), hold on
    
    % tidy up
    xlabel('Age [years]', 'FontSize', ftsize)
    ylabel('PP_{amp} [ratio]', 'FontSize', ftsize)
    set(gca, 'FontSize', ftsize, 'XTick', age_ticks)
    xlim([20,80])
    grid on
    box off
    
end
leg_labels = {'PP_b / (P1_a - DBP_a)', 'PP_b / (P2_a - DBP_a)', 'PP_b / (SBP_a - DBP_a) = PP_{amp}'};
legend(leg_labels, 'Location', 'NorthWest', 'FontSize', ftsize-4)
paper_size = [500,350];
PrintFigs(gcf, paper_size/70, [PATHS.CaseStudies, 'PPamp_amp_aug'])

%% Make plots of determinants

paper_size = [500,350];


figure('Position', [20 20 paper_size])
plot(d.PWV_a.*d.LVET/1000, d.PP_b./(d.P2pk_a-d.DBP_a), '.') % keep this
set(gca, 'FontSize', ftsize)
xlabel('Aortic PWV x LVET [m]', 'FontSize', ftsize)
ylabel('PP_b / (P2_a - DBP_a)', 'FontSize', ftsize)
box off
grid on
PrintFigs(gcf, paper_size/70, [PATHS.CaseStudies, 'PPamp_p2'])

figure('Position', [20 20 paper_size])
plot(d.dia_asc_a, d.PP_b./(d.P1in_a-d.DBP_a), '.') % keep this
set(gca, 'FontSize', ftsize)
xlabel('Ascending Aortic Diameter [mm]', 'FontSize', ftsize)
ylabel('PP_b / (P1_a - DBP_a)', 'FontSize', ftsize)
box off
grid on
PrintFigs(gcf, paper_size/70, [PATHS.CaseStudies, 'PPamp_p1'])

figure('Position', [20 20 paper_size])
plot(d.PWV_a.*d.LVET/1000, d.PP_amp, '.') % keep this
set(gca, 'FontSize', ftsize)
xlabel('Aortic PWV x LVET [m]', 'FontSize', ftsize)
ylabel('PP_{amp}', 'FontSize', ftsize)
box off
grid on
PrintFigs(gcf, paper_size/70, [PATHS.CaseStudies, 'PPamp'])


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