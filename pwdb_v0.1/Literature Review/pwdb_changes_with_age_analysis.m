function results_table = pwdb_changes_with_age_analysis
% PWDB_CHANGES_WITH_AGE_ANALYSIS runs an analysis of metadata from a literature
% review of articles describing changes in cardiovascular parameters with
% age.
%
%               pwdb_changes_with_age_analysis
%
%	Inputs:
%		pwdb_changes_with_age_review.txt  -  An input file containing the
%               metadata from the literature review, which is provided as
%               supplementary material with the article (detailed below).
%
%	Outputs:
%       results_table  -  A table containing the results of the literature
%               review.
%           
%   Accompanying Article:
%       This code is provided to facilitate reproduction of the literature
%       review analysis performed in:
%           Charlton P.H. et al. Modelling arterial pulse waves in healthy
%           ageing: a database for in silico evaluation of haemodynamics
%           and pulse wave indices, ~~ under review ~~  
%           DOI: ~~ tbc ~~
%       Further information on the literature review is provided in this
%       article and its supplementary material.
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
% v.1.0  Peter H. Charlton, King's College London%

% Output the licence details
provide_licence_details

% Setup universal parameters
up = setup_up;

% Load data from Excel Spreadsheet
data = load_data(up);

% Eliminate repeated data
data = eliminate_repetitions(data, up);

% Make initial table of literature review results
init_results_table = make_init_table_of_results(data, up);

% Make final table of literature review results
results_table = make_final_table_of_results(init_results_table, up);

end

function provide_licence_details

licence_details = ['\n\n        changes_with_age_analysis', ...
    '\n Copyright (C) 2018  King''s College London',...
    '\n This program comes with ABSOLUTELY NO WARRANTY', ... 
    '\n and is available under the GNU public license.\n'];

fprintf(licence_details)

end

function up = setup_up

close all

filename = 'pwdb_changes_with_age_review.txt';
up.paths.root_folder = [fileparts(mfilename('fullpath')), filesep];
up.paths.raw_data = [up.paths.root_folder, filename];
if ~exist(up.paths.raw_data, 'file')
    uiwait(msgbox(['At the next dialog box, please select the "' filename '" file which was provided as Supplementary Material']))
    [~,up.paths.root_folder] = uigetfile(up.paths.raw_data);
    if isnumeric(up.paths.root_folder) && up.paths.root_folder == 0
        error(['Couldn''t find the input data file (' filename '). Please download it from the article''s supplementary material, and try again.'])
    end
    up.paths.raw_data = [up.paths.root_folder, filename];
end
up.paths.savefolder = up.paths.root_folder;

up.curr_params_for_analysis = {'Dia - Aorta Abd', 'Dia - Aorta Asc', 'Dia - Aorta D Thor', 'Dia - Carotid', 'Ejection time', 'Heart rate', 'Cardiac output', 'Stroke volume'};
up.curr_params_for_change_with_age_plots = {'Stroke volume_pmd', 'Stroke volume_amd', 'Stroke volume_lit', 'Dia_desc_thor_pmd', 'Dia_desc_thor_amd', 'Dia_desc_thor_lit', 'Dia_abd_ao_pmd', 'Dia_abd_ao_amd', 'Dia_abd_ao_lit', 'Dia_asc_ao_pmd', 'Dia_asc_ao_amd', 'Dia_asc_ao_lit', 'Heart rate_pmd', 'Heart rate_amd'};
up.eps_figs = 0;
up.ftsize = 20;

up.mod_ages = 25:80;
up.mod_age_ticks = 20:10:up.mod_ages(end);
up.baseline_age = 25;
up.do_plot = 1;
up.color_2sd = 0.5*ones(1,3);
up.color_1sd = 0.75*ones(1,3);
up.mean_color = 'b';
up.ylim_offset = 0.1;

up.params.change_with_age = {'Dia - Aorta Abd', 'Dia - Aorta D Thor', 'Dia - Aorta Asc', 'Dia - Carotid', 'Length - Aorta Asc', 'Stroke volume', 'Heart rate', 'PWV - Aorta', 'PWV - Arm', 'PWV - Leg', 'Stiffness - Radial', 'Stiffness - Iliac', 'Stiffness - Femoral', 'Stiffness - Carotid', 'Stiffness - Brachial', 'Stiffness - Aorta D Thor', 'Stiffness - Aorta Asc', 'Stiffness - Aorta Abd', 'Sys Compliance', 'Sys Vasc Res', 'Length - Carotid', 'Length - Aorta D Thor', 'Length - Aorta Abd', 'Ejection time'};

end

function init_results_table = make_init_table_of_results(data, up)

fprintf('\n - Making initial table of results')

vars = unique(data.var);

%% Extract the number of articles which reported each type of change for each variable
for var_no = 1 : length(vars)
    
    curr_var = vars{var_no};
    rel_els = find(strcmp(data.var, curr_var));
    
    changes.n(var_no) = length(rel_els);
    changes.inc(var_no) = sum(strcmp(data.change(rel_els), 'increase'));
    changes.dec(var_no) = sum(strcmp(data.change(rel_els), 'decrease'));
    changes.none(var_no) = sum(strcmp(data.change(rel_els), 'none'));
    changes.non_linear(var_no) = sum(strcmp(data.change(rel_els), 'non_linear'));
   
end

%% Find percentage of articles which reported each type of change
temp = fieldnames(changes);
for field_no = 1 : length(temp)
    eval(['perc_changes.' temp{field_no} ' = 100*changes.' temp{field_no} './changes.n;']);
end
y = [perc_changes.dec(:), perc_changes.none(:), perc_changes.inc(:), perc_changes.non_linear(:)];


%% Create strings of articles which investigated each parameter
articles_strings = cell(length(vars),1);
for var_no = 1 : length(vars)
    rel_els = find(strcmp(data.var, vars{var_no}));
    articles_strings{var_no} = '';
    for article_no = 1 : length(rel_els)
        articles_strings{var_no} = [articles_strings{var_no}, data.article{rel_els(article_no)}, ','];
    end
    articles_strings{var_no} = articles_strings{var_no}(1:end-1);
end

%% Create table summarising literature review
var = vars;
n = changes.n';
inc_n = changes.inc';
inc_p = 100*inc_n./n;
dec_n = changes.dec';
dec_p = 100*dec_n./n;
none_n = changes.none';
none_p = 100*none_n./n;
non_linear_n = changes.non_linear';
non_linear_p = 100*non_linear_n./n;

init_results_table = table(var,n,none_p,inc_p,dec_p,non_linear_p,articles_strings);

end

function final_results_table = make_final_table_of_results(init_results_table, up)

fprintf('\n - Making final table of results')

req_params = {'Heart rate', 'Stroke volume', 'Ejection time', 'Peak flow time', 'Reverse Flow Volume', 'PWV Aorta', 'PWV - Arm', 'PWV - Leg', 'Dia - Aorta Asc', 'Dia - Aorta D Thor', 'Dia - Aorta Abd', 'Dia - Carotid', 'Dia - Iliac', 'Dia - Femoral', 'Dia - Brachial', 'Dia - Radial', 'Length - prox aorta', 'Length - dist aorta', 'Length - Carotid', 'Length - Iliac', 'Systemic Vascular Resistance', 'Systemic Vascular Compliance'};

for param_no = 1 : length(req_params)
    
    % Identify current parameter
    curr_param = req_params{param_no};
    
    % See if this corresponds to a single parameter in the literature review
    rel_row_el = find(strcmp(init_results_table.var, curr_param));
    if ~isempty(rel_row_el)
        
        temp = init_results_table(rel_row_el,:);
        articles_strings = strsplit(temp.articles_strings{1,1}, ',');
        temp.n_studies = temp.n;
        temp.n_articles = length(unique(articles_strings));
        temp = [temp(:,1),temp(:,[end-1,end]),temp(:,2:end-2)];
        temp.n = [];
        rel_row = temp;
    else
        % if not, then collect the data from the parameters which
        % correspond to this parameter.
        switch curr_param
            case 'PWV Aorta'
                rel_row_els = ~cellfun(@isempty, strfind(init_results_table.var, 'PWV - Aorta')) & ~strcmp(init_results_table.var, 'PWV - Aorta Leg');
            case 'Length - prox aorta'
                rel_row_els = ~cellfun(@isempty, strfind(init_results_table.var, 'Length - Aorta Asc')) | ~cellfun(@isempty, strfind(init_results_table.var, 'Length - Aorta Arch'));
            case 'Length - dist aorta'
                rel_row_els = ~cellfun(@isempty, strfind(init_results_table.var, 'Length - Aorta Abd')) | ~cellfun(@isempty, strfind(init_results_table.var, 'Length - Aorta D Thor')) | ~cellfun(@isempty, strfind(init_results_table.var, 'Length - Aorta Desc')) | ~cellfun(@isempty, strfind(init_results_table.var, 'Length - Aorta Thor'));
        end
        rel_rows = init_results_table(rel_row_els,:);
        
        % collect these rows into a single row
        var = {curr_param};
        n_studies = sum(rel_rows.n);
        none_p = 100*sum(rel_rows.n.*rel_rows.none_p/100)/n_studies;
        inc_p = 100*sum(rel_rows.n.*rel_rows.inc_p/100)/n_studies;
        dec_p = 100*sum(rel_rows.n.*rel_rows.dec_p/100)/n_studies;
        non_linear_p = 100*sum(rel_rows.n.*rel_rows.non_linear_p/100)/n_studies;
        all_articles = [];
        for s = 1 : height(rel_rows)
            all_articles = [all_articles, strsplit(rel_rows.articles_strings{s}, ',')];
        end
        all_articles = unique(all_articles);
        articles_strings = '';
        for s= 1 : length(all_articles)
            articles_strings = [articles_strings, all_articles{s}, ','];
        end
        articles_strings = {articles_strings(1:end-1)};
        n_articles = length(all_articles);
        rel_row = table(var, n_studies, n_articles, none_p, inc_p, dec_p, non_linear_p, articles_strings);
        clear var n_studies n_articles none_p inc_p dec_p non_linear_p articles_strings rel_rows rel_row_els curr_param
    end
    
    final_results_table(param_no,:) = rel_row;
    clear rel_row rel_row_el
end
clear param_no

% Rename some variables
final_results_table.var = strrep(final_results_table.var, 'PWV - Aorta Leg', 'PWV - Lower Limb');
final_results_table.var = strrep(final_results_table.var, 'PWV - Aorta Arm', 'PWV - Upper Limb');
final_results_table.var = strrep(final_results_table.var, 'Dia - Aorta Asc', 'Dia - Asc Aorta');
final_results_table.var = strrep(final_results_table.var, 'Dia - Aorta D Thor', 'Dia - Desc Thor Aorta');
final_results_table.var = strrep(final_results_table.var, 'Dia - Aorta Abd', 'Dia - Abd Aorta');

end

function rows = load_data(up)

fprintf('\n - Loading data from input file')

filepath = up.paths.raw_data;
if strcmp(filepath(end-2:end),'xls') || strcmp(filepath(end-3:end),'xlsx')
    % If the raw data is presented in Excel Format
    [~, ~, data] = xlsread(up.paths.raw_data);
    headers = data(1,:);
    rel_col = find(strcmp(headers, 'Article'));
    rows.article = data(2:end, rel_col);
    rel_col = find(strcmp(headers, 'Dependent Variable'));
    rows.var = data(2:end, rel_col);
    rel_col = find(strcmp(headers, 'Subgroup'));
    rows.subgroup = data(2:end, rel_col);
    rel_col = find(strcmp(headers, 'Significant Change'));
    rows.change = data(2:end, rel_col);
    rel_col = find(strcmp(headers, 'No subjects'));
    rows.n = data(2:end, rel_col);
    rows.n = convert_to_num(rows.n);
else
    % If it is tab-delimited format
    data = tdfread(up.paths.raw_data);
    rows.article = cellstr(data.Article);
    rows.var = cellstr(data.Dependent_Variable);
    rows.subgroup = cellstr(data.Subgroup);
    rows.change = cellstr(data.Significant_Change);
    rows.n = data.No_subjects;
end

end

function data = eliminate_repetitions(data, up)

fprintf('\n - Eliminating repetitions')

% eliminate rows
elim_rows = [];
% eliminate repetition due to difference between compliance and distensibility (for reasons given in article discussion)
curr_el = find(strcmp(data.article, 'VanderHeijden-Spek2000') & strcmp(data.subgroup, 'compliance'));
elim_rows = [elim_rows; curr_el];
% eliminate repetition due to difference between women and men (taking men)
curr_el = find(strcmp(data.article, 'Fleg1995') & strcmp(data.subgroup, 'women'));
elim_rows = [elim_rows; curr_el];
% eliminate repetitions due to differences between intima and media thickness
curr_el = find(strcmp(data.article, 'Virmani1991') & ~cellfun(@isempty, strfind(data.subgroup, 'media')));
elim_rows = [elim_rows; curr_el];
% eliminate repetitions due to different methods for calculating ejection time (corrected or not)
curr_el = find(strcmp(data.article, 'Gardin1987') & ~cellfun(@isempty, strfind(data.subgroup, 'corrected')));
elim_rows = [elim_rows; curr_el];
curr_el = find(strcmp(data.article, 'Shaw1973') & ~cellfun(@isempty, strfind(data.subgroup, 'Q-wave to PCG sound')));
elim_rows = [elim_rows; curr_el];
% eliminate repetitions due to different methods for measuring ejection time (from pressure wave or echo)
curr_el = find(strcmp(data.article, 'Salvi2018') & ~cellfun(@isempty, strfind(data.subgroup, 'carotid')));
elim_rows = [elim_rows; curr_el];



keep_rows = setxor(1:length(data.article), elim_rows);
data_fields = fieldnames(data);
for field_no = 1 : length(data_fields)
    eval(['data.' data_fields{field_no} ' = data.' data_fields{field_no} '(keep_rows);']);    
end

clear elim_rows keep_rows

% identify entries which are repetitions from the same article.
for s = 1 : length(data.article)
    temp{s} = [data.article{s}, data.var{s}];
end
temp = temp(:);

elim_rows = [];
for s = 2 : length(temp)
    repeat_els = find(strcmp(temp{s}, temp(1:s-1)));
    if ~isempty(repeat_els)
        cand_changes = [data.change{s}; data.change(repeat_els)];
        if length(unique(cand_changes)) > 1
            error('Do something about this');
        end
        elim_rows = [elim_rows; s];
    end
end

% eliminate repetitions
keep_rows = setxor(1:length(temp), elim_rows);
data_fields = fieldnames(data);
for field_no = 1 : length(data_fields)
    eval(['data.' data_fields{field_no} ' = data.' data_fields{field_no} '(keep_rows);']);    
end

end

function new = convert_to_num(old)

% From: https://uk.mathworks.com/matlabcentral/answers/30547-convert-cell-to-matrix-with-mixed-data-types

old(cellfun(@ischar,old)) = {NaN};
new = cell2mat(old);

end
