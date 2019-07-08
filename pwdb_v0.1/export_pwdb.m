function export_pwdb(pwdb_no)
% EXPORT_PWDB exports the pulse wave database in several formats for
% further research.
%
%               export_pwdb
%
%   Inputs:     the 'collated_data.mat', 'haemodynamic_params.mat',
%               'pulse_wave_vels.mat', 'pulse_wave_inds.mat' and
%               'system_chars' files produced by
%               'extract_pwdb_simulation_data.m' and 'pwdb_pwa.m'. 
%
%   Outputs:    - A range of files containing the database in various
%                   formats, stored in the 'exported_data' folder (within
%                   the 'processed_data' folder):
%                      - Matlab format (both a single file and multiple, smaller, files)
%                      - CSV format
%                      - WFDB (waveform database) format, as used in PhysioNet databases
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
%     Permission is hereby granted, free of charge, to any person obtaining a copy
%     of this software and associated documentation files (the "Software"), to deal
%     in the Software without restriction, including without limitation the rights
%     to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%     copies of the Software, and to permit persons to whom the Software is
%     furnished to do so, subject to the following conditions:
%     
%     The above copyright notice and this permission notice shall be included in all
%     copies or substantial portions of the Software.
%     
%     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%     IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%     AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%     LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%     SOFTWARE.
%
%  Copyright (C) 2019  King's College London
% 
% v.1.0    Contributed to by: Marie Willemet, Jordi Alastruey, and Peter H. Charlton

fprintf('\n --- Exporting PWDB in different formats ---')

% Setup paths with current simulation paths
PATHS = setup_paths_for_post_processing(pwdb_no);

% create folders in which to store exported data
create_folders(PATHS)

% load collated data
load(PATHS.collated_data)

% Create files containing PWs
create_pw_files(PATHS, collated_data);

% Create files containing model configs
create_model_config_file(PATHS, collated_data);

% Create files containing model variations
create_model_variation_file(PATHS, collated_data);

% Create files containing haemodynamic params
create_haemod_param_file(PATHS);

% Create files containing model geometry
create_model_geometry_files(PATHS, collated_data);

% Create file containing pulse onset times
create_pulse_onsets_file(PATHS, collated_data);

% Create files containing pulse wave indices
pw_inds = create_pw_ind_file(PATHS, collated_data);

% Creat Matlab file containing all data
create_mat_file_w_all_data(PATHS, collated_data, pw_inds);

fprintf('\n --- Finished exporting PWDB ---')

end

function create_folders(PATHS)

folders_to_make = {PATHS.exported_data, PATHS.exported_data_mat, PATHS.exported_data_csv, PATHS.exported_data_wfdb, PATHS.exported_data_geo};
for s = 1 : length(folders_to_make)
    if ~exist(folders_to_make{s}, 'dir')
        mkdir(folders_to_make{s})
    end
end

end

function create_pw_files(PATHS, collated_data)

fprintf('\n - Creating Pulse Wave Files')
%% Setting up

% settings
sites = {'AorticRoot',   'ThorAorta',  'AbdAorta',  'IliacBif', 'Carotid',  'SupTemporal', 'SupMidCerebral', 'Brachial', 'Radial', 'Digital', 'CommonIliac', 'Femoral', 'AntTibial'};
site_domain_no = [1,            18,          39,          41,        15,             87,                 72,         21,       22,       112,             44,        46,           49];
site_dist_prop = [0,             1,           0,           1,       0.5,              1,                  1,       0.75,        1,         1,            0.5,       0.5,            1];
signals = {'P', 'U', 'A', 'PPG'};   % omit 'Q' as it's U.*A

% identify domains which are available
domains = extractfield(collated_data(1).output_data, 'domain_no');
[~,rel_els,~] = intersect(site_domain_no, domains);
site_domain_no = site_domain_no(rel_els);
sites = sites(rel_els);
site_dist_prop = site_dist_prop(rel_els);

%% Creating pulse wave files (WFDB)
% - WFDB: one PW file per subject
fprintf(': WFDB')
up.dataset_name = PATHS.pwdb_filename_prefix;
fs = collated_data(1).output_data(1).fs;

cd(PATHS.exported_data_wfdb)
% cycle through each subject
no_subjs = length(collated_data);
file_names = cell(no_subjs,1);
for subj_no = 1 : no_subjs
    sig_names = cell(0);
    
    % setup command
    file_names{subj_no} = [up.dataset_name, sprintf('%.4d', subj_no)];
        
    % Extract signals data
    data_mat = nan(2000, 1000); max_no_samps = 0; no_cols = 0;
    for site_no = 1 : length(sites)
        
        curr_site = sites{site_no};
        
        % Identify site-specific info
        curr_domain_no = site_domain_no(site_no);
        curr_domain_el = find(domains == curr_domain_no);
        curr_site_dist_prop = site_dist_prop(site_no);
        
        % Identify the relevant data for this simulation at this site
        rel_wave_data = collated_data(subj_no).output_data(curr_domain_el);
        seg_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(curr_domain_no);
        dists = extractfield(rel_wave_data, 'distances');
        [~,rel_dist_el] = min(abs(curr_site_dist_prop*seg_length-dists)); clear dists
        
        % Extract each signal at this site
        for signal_no = 1 : length(signals)
            if strcmp(signals{signal_no}, 'P')
                factor = 133.33;
            else
                factor = 1;
            end
            curr_sig = signals{signal_no};
            eval(['curr_sig_data = rel_wave_data.' curr_sig '(:,rel_dist_el)/factor;'])
            sig_names{end+1} = [curr_site '_' curr_sig];
            data_mat(1:length(curr_sig_data),no_cols+1) = curr_sig_data;
            max_no_samps = max(max_no_samps, length(curr_sig_data));
            no_cols = no_cols+1;
            clear factor curr_sig curr_sig_data
        end
        clear signal_no rel_wave_data seg_length rel_dist_el curr_domain_no curr_domain_el curr_site_dist_prop curr_site
        
    end
    clear site_no
    
    % Get rid of excess NaNs in data_mat
    data_mat = data_mat(:,1:no_cols); clear no_cols
    data_mat = data_mat(1:max_no_samps,:); clear max_no_samps
    
%     % Get rid of last sample if this was a nan in some recordings
%     if sum(isnan(data_mat(end,:))) > 0
%         data_mat = data_mat(1:end-1,:);
%     end
%     if sum(sum(isnan(data_mat)))
%         error('There are still some nans')
%     end
    
    % Create signals file
    descrip = ['<age>: ' num2str(collated_data(subj_no).input_data.sim_settings.age) ' <sex>: male'];
    units = '';
    for sig_no = 1 : length(sig_names)
        signal_names(sig_no,1) = {[sig_names{sig_no} ', ' ]};
        if ~isempty(strfind(sig_names{sig_no}, '_PPG'))
            rel_unit = 'au';
        elseif ~isempty(strfind(sig_names{sig_no}, '_P'))
            rel_unit = 'mmHg';
        elseif ~isempty(strfind(sig_names{sig_no}, '_A'))
            rel_unit = 'm2';
        elseif ~isempty(strfind(sig_names{sig_no}, '_Q'))
            rel_unit = 'm3_per_sec';
        elseif ~isempty(strfind(sig_names{sig_no}, '_U'))
            rel_unit = 'm_per_s';
        end
        units = [units, rel_unit, '/']; clear rel_unit
    end
    clear sig_no
    units = units(1:end-1);
    % convert this subject's data into WFDB format
    mat2wfdb(data_mat, file_names{subj_no}, fs, [], units, descrip, [], signal_names);
    clear data_mat units descrip signal_names sig_names type subtype 
    
end
clear no_subjs subj_no fs

% create a RECORDS file
file_name = 'RECORDS';
file_path = [PATHS.exported_data_wfdb, file_name];
fid = fopen(file_path, 'w');
formatSpec = '%s\r\n';
for line_no = 1 : length(file_names)
    fprintf(fid, formatSpec, file_names{line_no});
end

% create a DBS file
file_name = 'DBS';
file_path = [PATHS.exported_data_wfdb, file_name];
fid = fopen(file_path, 'w');
fprintf(fid, [up.dataset_name, '\t', 'Preliminary Pulse Wave DataBase Dataset']);

fclose all;    


%% Creating pulse wave files (CSV and Mat)
% - Mat and CSV format: one PW file per site
fprintf(', CSV and Mat')

% Cycle through each arterial site
for site_no = 1 : length(sites)
    
    % Identify site-specific info
    curr_domain_no = site_domain_no(site_no);
    curr_domain_el = find(domains == curr_domain_no);
    curr_site_dist_prop = site_dist_prop(site_no);
    
    %% Extract PWs at this site from all the simulations
    
    % Cycle through each simulation
    for sim_no = 1 : length(collated_data)
        
        % Identify the relevant data for this simulation at this site
        rel_wave_data = collated_data(sim_no).output_data(curr_domain_el);
        seg_length = collated_data(sim_no).input_data.sim_settings.network_spec.length(curr_domain_no);
        dists = extractfield(rel_wave_data, 'distances');
        [~,rel_dist_el] = min(abs(curr_site_dist_prop*seg_length-dists)); clear dists
        
        % Extract each signal
        for signal_no = 1 : length(signals)
            if strcmp(signals{signal_no}, 'P')
                factor = 133.33;
            else
                factor = 1;
            end
            curr_sig = signals{signal_no};
            
            % store signal
            eval(['PWs.' curr_sig '{sim_no} = rel_wave_data.' curr_sig '(:,rel_dist_el)/factor;'])
            
            % calculate pulse onset time
            rel_row = find(domains == 1);
            aortic_root_onset = collated_data(sim_no).output_data(rel_row).start_sample(1);
            if strcmp(curr_sig, 'PPG')
                temp_onset_time = (rel_wave_data.PPG_start_sample(rel_dist_el)-aortic_root_onset)/rel_wave_data.fs;
            else
                temp_onset_time = (rel_wave_data.start_sample(rel_dist_el)-aortic_root_onset)/rel_wave_data.fs;
            end
            % if the time is negative then add on the duration of a cardiac cycle
            if temp_onset_time < 0
                cardiac_cycle_duration = length(rel_wave_data.P(:,1))/rel_wave_data.fs;
                temp_onset_time = temp_onset_time + cardiac_cycle_duration;
            end
            eval(['PWs.onset_times.' curr_sig '(sim_no,1) = temp_onset_time;']); clear temp_onset_time
            
            clear factor curr_sig
        end
        clear signal_no rel_wave_data seg_length rel_dist_el
    end
    clear sim_no curr_site_dist_prop
    
    %% Store PWs for this site for all simulations (CSV and Mat)
        
    % Create filename for this site
    filename = ['PWs_' sites{site_no}];
    
    % Store as CSV and MAT
    file_types = {'csv','mat'};
    sig_names = fieldnames(PWs);
    sig_names = sig_names(~strcmp(sig_names, 'onset_times'));
    sig_names = sig_names(~strcmp(sig_names, 'fs'));
    sig_names = sig_names(~strcmp(sig_names, 'units'));
    for type_no = 1 : length(file_types)
        curr_type = file_types{type_no};
        eval(['filepath = [PATHS.exported_data_' curr_type ', filename, ''.'', curr_type];']);
        switch curr_type
            case 'csv'
                
                % Create CSV file for each signal
                for sig_no = 1 : length(sig_names)
                    curr_sig = sig_names{sig_no};
                    curr_filepath = [filepath(1:end-4), '_', curr_sig, filepath(end-3:end)];
                    
                    % Extract data
                    max_no_samps = max(cellfun(@length, PWs.P));
                    data_mat = nan(length(PWs.P),max_no_samps);
                    for subj_no = 1 : length(PWs.P)
                        eval(['data_mat(subj_no,1:length(PWs.P{subj_no})) = PWs.' curr_sig '{subj_no};']);
                    end
                    clear subj_no max_no_samps
                    
                    % Add in subject numbers
                    temp = 1 :length(collated_data);
                    data_mat = [temp(:), data_mat]; clear temp
                    
                    % Generate header line
                    header_line = 'Subject Number, ';
                    for col_no = 1 : (size(data_mat,2)-1)
                        curr_col = ['pt' num2str(col_no)];
                        header_line = [header_line, curr_col, ', '];                        
                    end
                    clear col_no curr_col
                    header_line = header_line(1:end-2);
                    
                    % Write header line to file
                    fid = fopen(curr_filepath,'w');
                    up.csv.new_line = '\n';
                    fprintf(fid,[header_line, up.csv.new_line]); clear header_line
                    fclose(fid); clear fid
                    
                    % Write to CSV
                    dlmwrite(curr_filepath, data_mat, '-append');
                    clear data_mat curr_filepath curr_sig header_line curr_col col_no
                    
                end
                clear sig_no curr_sig
                
            case 'mat'
                
                % Save to Mat file
                PWs.fs = collated_data(1).output_data(1).fs;
                PWs.units.P = 'mmHg';
                PWs.units.U = 'm/s';
                PWs.units.A = 'm3';
                PWs.units.PPG = 'au';
                save(filepath, 'PWs');
        end
        clear curr_type
    end
    clear type_no
    
end
clear site_no curr_domain_el curr_domain_no

end

function sig_header = find_sig_header(curr_sig)

switch curr_sig
    case 'P'
        sig_header = 'Pressure [mmHg]';
    case 'Q'
        sig_header = 'Flow rate [m3/s]';
    case 'U'
        sig_header = 'Flow velocity [m/s]';
    case 'A'
        sig_header = 'Luminal area [m2]';
    case 'PPG'
        sig_header = 'Photoplethysmogram [unitless]';
end
end

function pw_inds = create_pw_ind_file(PATHS, collated_data)

fprintf('\n - Creating pulse wave indices file: ')

% Load data
load(PATHS.pulse_wave_inds)

% setup
sites = {'AorticRoot',   'ThorAorta',  'AbdAorta',  'IliacBif', 'Carotid',  'SupTemporal', 'SupMidCerebral', 'Brachial', 'Radial', 'Digital', 'CommonIliac', 'Femoral', 'AntTibial'};
site_domain_no = [1,            18,          39,          41,        15,             87,                 72,         21,       22,       112,             44,        46,           49];
site_dist_prop = [0,             1,           0,           1,       0.5,              1,                  1,       0.75,        1,         1,            0.5,       0.5,            1];
site_central_log = [1,           1,           1,           1,         0,              0,                  0,          0,        0,         0,              1,         0,            0]; 
domain_nos = extractfield(collated_data(1).output_data, 'domain_no');

% identify domains which are available
domains = extractfield(collated_data(1).output_data, 'domain_no');
[~,rel_els,~] = intersect(site_domain_no, domains);
site_domain_no = site_domain_no(rel_els);
sites = sites(rel_els);
site_dist_prop = site_dist_prop(rel_els);
site_central_log = site_central_log(rel_els);

% Identify parameters, and extract data
params_at_each_site = {'SBP', 'DBP', 'MBP', 'PP', 'Qmax', 'Qmin', 'Qmean', 'Qtotal', 'Umax', 'Umin', 'Umean', 'Amax', 'Amin', 'Amean', 'P1in', 'P2in', 'P1pk', 'P2pk', 'Psys', 'Pms', 'PPGa', 'PPGb', 'PPGc', 'PPGd', 'PPGe', 'PPGsys', 'PPGdia', 'PPGdic', 'PPGms', 'AI', 'AP', 'RI', 'SI', 'AGI_mod', 'PTT','VP1in','VP2in','VP1pk','VP2pk','UP1in','UP2in','UP1pk','UP2pk'};
params = {'subj_no', 'age'};
for site_no = 1 : length(sites)
    fprintf([sites{site_no}, ' ']);
    domain_el = find(domain_nos == site_domain_no(site_no));
    % cycle through each parameter
    for site_param_no = 1 : length(params_at_each_site)
        % skip if these params are only to be calculated for aortic root
        if ~strcmp(sites{site_no}, 'AorticRoot')
            if sum(strcmp(params_at_each_site{site_param_no}, {'VP1in','VP2in','VP1pk','VP2pk','UP1in','UP2in','UP1pk','UP2pk'}))
                continue
            end
        end
        % skip if this is the PPG at a non-peripheral site
        if site_central_log(site_no) == 1
            if length(params_at_each_site{site_param_no})>2 && strcmp(params_at_each_site{site_param_no}(1:3), 'PPG')
                continue
            end
            if strcmp(params_at_each_site{site_param_no}(1:2), 'RI') || ...
                    strcmp(params_at_each_site{site_param_no}(1:2), 'SI') || ...
                    (length(params_at_each_site{site_param_no})>6 && strcmp(params_at_each_site{site_param_no}, 'AGI_mod'))
                continue
            end
        end
        
        % note down this parameter's name
        params{end+1} = [sites{site_no}, '_', params_at_each_site{site_param_no}, '_V'];
        params{end+1} = [sites{site_no}, '_', params_at_each_site{site_param_no}, '_T'];
        % extract data for it
        for subj_no = 1 : length(collated_data)
            param_data(subj_no, 1) = subj_no;
            param_data(subj_no, 2) = collated_data(subj_no).input_data.sim_settings.age;
            switch params_at_each_site{site_param_no}
                case 'SBP'
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).P(:,rel_dist_el);
                    [~,el] = max(wav);
                    v = wav(el)/133.33;
                    t = (el-1)/collated_data(subj_no).output_data(domain_el).fs;
                case 'DBP'
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).P(:,rel_dist_el);
                    v = min(wav)/133.33;
                    t = -1;
                case 'MBP'
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).P(:,rel_dist_el);
                    v = mean(wav)/133.33;
                    t = -1;
                case 'PP'
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).P(:,rel_dist_el);
                    v = (max(wav)-min(wav))/133.33;
                    t = -1;
                case 'Qmax'
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).U(:,rel_dist_el).*collated_data(subj_no).output_data(domain_el).A(:,rel_dist_el);
                    [~,el] = max(wav);
                    v = wav(el);
                    t = (el-1)/collated_data(subj_no).output_data(domain_el).fs;
                case 'Qmin'
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).U(:,rel_dist_el).*collated_data(subj_no).output_data(domain_el).A(:,rel_dist_el);
                    [~,el] = min(wav);
                    v = wav(el);
                    t = (el-1)/collated_data(subj_no).output_data(domain_el).fs;
                case 'Qmean'
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).U(:,rel_dist_el).*collated_data(subj_no).output_data(domain_el).A(:,rel_dist_el);
                    v = mean(wav);
                    t = -1;
                case 'Qtotal'
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).U(:,rel_dist_el).*collated_data(subj_no).output_data(domain_el).A(:,rel_dist_el);
                    v = sum(wav)/collated_data(subj_no).output_data(domain_el).fs;
                    t = -1;
                case 'Umax'
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).U(:,rel_dist_el);
                    [~,el] = max(wav);
                    v = wav(el);
                    t = (el-1)/collated_data(subj_no).output_data(domain_el).fs;
                case 'Umin'
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).U(:,rel_dist_el);
                    [~,el] = min(wav);
                    v = wav(el);
                    t = (el-1)/collated_data(subj_no).output_data(domain_el).fs;
                case 'Umean'
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).U(:,rel_dist_el);
                    v = mean(wav);
                    t = -1;
                case 'Amax'
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).A(:,rel_dist_el);
                    [~,el] = max(wav);
                    v = wav(el);
                    t = (el-1)/collated_data(subj_no).output_data(domain_el).fs;
                case 'Amin'
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).A(:,rel_dist_el);
                    [~,el] = min(wav);
                    v = wav(el);
                    t = -1;
                case 'Amean'
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).A(:,rel_dist_el);
                    v = mean(wav);
                    t = -1;
                case 'VP1in'
                    % find Q waveform
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).A(:,rel_dist_el).*collated_data(subj_no).output_data(domain_el).U(:,rel_dist_el);
                    % find P1in point
                    pwa_pt = 'p1in'; pwa_sig = 'P';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    samp_no = (t*collated_data(subj_no).output_data(domain_el).fs)+1;
                    % find volume ejected up to P1
                    try v = sum(wav(1:samp_no)); catch, v = nan; end
                    t = -1;
                case 'VP2in'
                    % find Q waveform
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).A(:,rel_dist_el).*collated_data(subj_no).output_data(domain_el).U(:,rel_dist_el);
                    % find P2in point
                    pwa_pt = 'p2in'; pwa_sig = 'P';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    samp_no = (t*collated_data(subj_no).output_data(domain_el).fs)+1;
                    % find volume ejected up to P1
                    try v = sum(wav(1:samp_no)); catch, v = nan; end
                    t = -1;
                case 'VP1pk'
                    % find Q waveform
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).A(:,rel_dist_el).*collated_data(subj_no).output_data(domain_el).U(:,rel_dist_el);
                    % find P1pk point
                    pwa_pt = 'p1pk'; pwa_sig = 'P';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    samp_no = (t*collated_data(subj_no).output_data(domain_el).fs)+1;
                    % find volume ejected up to P1
                    try v = sum(wav(1:samp_no)); catch, v = nan; end
                    t = -1;
                case 'VP2pk'
                    % find Q waveform
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).A(:,rel_dist_el).*collated_data(subj_no).output_data(domain_el).U(:,rel_dist_el);
                    % find P2pk point
                    pwa_pt = 'p2pk'; pwa_sig = 'P';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    samp_no = (t*collated_data(subj_no).output_data(domain_el).fs)+1;
                    % find volume ejected up to P1
                    try v = sum(wav(1:samp_no)); catch, v = nan; end
                    t = -1;
                case 'UP1in'
                    % find U waveform
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).U(:,rel_dist_el);
                    % find P1in point
                    pwa_pt = 'p1in'; pwa_sig = 'P';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    samp_no = (t*collated_data(subj_no).output_data(domain_el).fs)+1;
                    % find U at this point
                    try, v = wav(samp_no); catch v = nan; end
                    t = -1;
                case 'UP2in'
                    % find U waveform
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).U(:,rel_dist_el);
                    % find P2in point
                    pwa_pt = 'p2in'; pwa_sig = 'P';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    samp_no = (t*collated_data(subj_no).output_data(domain_el).fs)+1;
                    % find U at this point
                    try, v = wav(samp_no); catch v = nan; end
                    t = -1;                 
                case 'UP1pk'
                    % find U waveform
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).U(:,rel_dist_el);
                    % find P1pk point
                    pwa_pt = 'p1pk'; pwa_sig = 'P';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    samp_no = (t*collated_data(subj_no).output_data(domain_el).fs)+1;
                    % find U at this point
                    try, v = wav(samp_no); catch v = nan; end 
                    t = -1;                
                case 'UP2pk'
                    % find U waveform
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    wav = collated_data(subj_no).output_data(domain_el).U(:,rel_dist_el);
                    % find P2pk point
                    pwa_pt = 'p2pk'; pwa_sig = 'P';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    samp_no = (t*collated_data(subj_no).output_data(domain_el).fs)+1;
                    % find U at this point
                    try, v = wav(samp_no); catch v = nan; end
                    t = -1;                 
                case 'P1in'
                    pwa_pt = 'p1in'; pwa_sig = 'P';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    v = find_pw_point_value(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                case 'P2in'
                    pwa_pt = 'p2in'; pwa_sig = 'P';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    v = find_pw_point_value(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                case 'P1pk'
                    pwa_pt = 'p1pk'; pwa_sig = 'P';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    v = find_pw_point_value(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                case 'P2pk'
                    pwa_pt = 'p2pk'; pwa_sig = 'P';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    v = find_pw_point_value(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                case 'Psys'
                    pwa_pt = 's'; pwa_sig = 'P';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    v = find_pw_point_value(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                case 'Pms'
                    pwa_pt = 'ms'; pwa_sig = 'P';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    rel_cv_ind = 'ms';
                    rel_row = find(strcmp(pulse_wave_inds(1).cv_ind_names, rel_cv_ind));
                    eval(['v = pulse_wave_inds(subj_no).' pwa_sig '_pwa(domain_el).cv_inds(rel_row,rel_dist_el);'])
                    v = v/133.33;
                case 'PPGms'
                    pwa_pt = 'ms'; pwa_sig = 'PPG';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    rel_cv_ind = 'ms';
                    rel_row = find(strcmp(pulse_wave_inds(1).cv_ind_names, rel_cv_ind));
                    eval(['v = pulse_wave_inds(subj_no).' pwa_sig '_pwa(domain_el).cv_inds(rel_row,rel_dist_el);'])
                case 'Pdic'
                    pwa_pt = 'dic'; pwa_sig = 'P';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    v = find_pw_point_value(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                case 'PPGa'
                    pwa_pt = 'a'; pwa_sig = 'PPG';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    rel_cv_ind = 'a_div_amp';
                    rel_row = find(strcmp(pulse_wave_inds(1).cv_ind_names, rel_cv_ind));
                    eval(['v = pulse_wave_inds(subj_no).' pwa_sig '_pwa(domain_el).cv_inds(rel_row,rel_dist_el);'])
                case 'PPGb'
                    pwa_pt = 'b'; pwa_sig = 'PPG';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    rel_cv_ind = 'b_div_amp';
                    rel_row = find(strcmp(pulse_wave_inds(1).cv_ind_names, rel_cv_ind));
                    eval(['v = pulse_wave_inds(subj_no).' pwa_sig '_pwa(domain_el).cv_inds(rel_row,rel_dist_el);'])
                case 'PPGc'
                    pwa_pt = 'c'; pwa_sig = 'PPG';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    rel_cv_ind = 'c_div_amp';
                    rel_row = find(strcmp(pulse_wave_inds(1).cv_ind_names, rel_cv_ind));
                    eval(['v = pulse_wave_inds(subj_no).' pwa_sig '_pwa(domain_el).cv_inds(rel_row,rel_dist_el);'])
                case 'PPGd'
                    pwa_pt = 'd'; pwa_sig = 'PPG';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    rel_cv_ind = 'd_div_amp';
                    rel_row = find(strcmp(pulse_wave_inds(1).cv_ind_names, rel_cv_ind));
                    eval(['v = pulse_wave_inds(subj_no).' pwa_sig '_pwa(domain_el).cv_inds(rel_row,rel_dist_el);'])
                case 'PPGe'
                    pwa_pt = 'e'; pwa_sig = 'PPG';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    rel_cv_ind = 'e_div_amp';
                    rel_row = find(strcmp(pulse_wave_inds(1).cv_ind_names, rel_cv_ind));
                    eval(['v = pulse_wave_inds(subj_no).' pwa_sig '_pwa(domain_el).cv_inds(rel_row,rel_dist_el);'])
                case 'PPGsys'
                    pwa_pt = 's'; pwa_sig = 'PPG';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    v = find_pw_point_value(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                case 'PPGdic'
                    pwa_pt = 'dic'; pwa_sig = 'PPG';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    v = find_pw_point_value(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                case 'PPGdia'
                    pwa_pt = 'dia'; pwa_sig = 'PPG';
                    t = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                    v = find_pw_point_value(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig);
                case 'AI'
                    pwa_sig = 'P';
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    rel_cv_ind = 'AI';
                    rel_row = find(strcmp(pulse_wave_inds(1).cv_ind_names, rel_cv_ind));
                    eval(['v = pulse_wave_inds(subj_no).' pwa_sig '_pwa(domain_el).cv_inds(rel_row,rel_dist_el);'])
                    t = -1;
                case 'AP'
                    pwa_sig = 'P';
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    rel_cv_ind = 'AP';
                    rel_row = find(strcmp(pulse_wave_inds(1).cv_ind_names, rel_cv_ind));
                    eval(['v = pulse_wave_inds(subj_no).' pwa_sig '_pwa(domain_el).cv_inds(rel_row,rel_dist_el);'])
                    v = v/133.33;
                    t = -1;
                case 'RI'
                    pwa_sig = 'PPG';
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    rel_cv_ind = 'RI';
                    rel_row = find(strcmp(pulse_wave_inds(1).cv_ind_names, rel_cv_ind));
                    eval(['v = pulse_wave_inds(subj_no).' pwa_sig '_pwa(domain_el).cv_inds(rel_row,rel_dist_el);'])
                    t = -1;
                case 'SI'
                    pwa_sig = 'PPG';
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    rel_cv_ind = 'SI';
                    rel_row = find(strcmp(pulse_wave_inds(1).cv_ind_names, rel_cv_ind));
                    eval(['v = pulse_wave_inds(subj_no).' pwa_sig '_pwa(domain_el).cv_inds(rel_row,rel_dist_el);'])
                    t = -1;
                case 'AGI_mod'
                    pwa_sig = 'PPG';
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    rel_cv_ind = 'AGI_mod';
                    rel_row = find(strcmp(pulse_wave_inds(1).cv_ind_names, rel_cv_ind));
                    eval(['v = pulse_wave_inds(subj_no).' pwa_sig '_pwa(domain_el).cv_inds(rel_row,rel_dist_el);'])
                    t = -1;
                case 'PTT'
                    desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
                    [~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));
                    arrival_v = collated_data(subj_no).output_data(domain_el).start_sample(rel_dist_el);
                    start_v = collated_data(subj_no).output_data(1).start_sample(1);
                    v = (arrival_v - start_v)/collated_data(subj_no).output_data(1).fs;
                    t = -1;
            end
            
            % store value
            param_data(subj_no,size(params,2)-1) = v;
            if t~=-1
                param_data(subj_no,size(params,2)) = t;
            elseif subj_no == length(collated_data)
                params = params(1:end-1);
                params{end} = params{end}(1:end-2);
            end
            clear v pwa_sig desired_length rel_dist_el arrival_v start_v v t rel_cv_ind rel_row wav el
            
        end
    end
end

% Make header line
header_line = 'Subject Number, Age, ';
for param_no = 3 : length(params)
    curr_param = params{param_no};
    % Make Header
    header_line = [header_line, curr_param, ', '];
end
clear param_no curr_param label units abbr graph_title
header_line = header_line(1:end-2);

% Write header line to file
curr_filename = [PATHS.exported_data, PATHS.pwdb_filename_prefix, '_pw_indices.csv'];
fid = fopen(curr_filename,'w');
up.csv.new_line = '\n';
fprintf(fid,[header_line, up.csv.new_line]); clear header_line
fclose(fid); clear fid

% Write parameter data to file
dlmwrite(curr_filename, param_data, '-append');

% Output results (to save in single Matlab file
pw_inds = array2table(param_data, 'VariableNames',params);

end

function create_model_config_file(PATHS, collated_data)

fprintf('\n - Creating model configurations file')

% Setup
params.vars = {'subject_ID', 'base', 'base_age'};
params.sim_settings = {'age', 'hr', 'sv', 't_pf', 'reg_vol', 'dbp', 'mbp', 'mu', 'alpha', 'p_drop', 'pvc', 'p_out', 'rho', 'lvet', 'pvr', 'gamma_b0', 'gamma_b1'};
params.network_spec = {'k1', 'k2', 'k3'};
params.all = [params.vars, params.sim_settings, params.network_spec];

% Check names of gamma variables (this is required for the ppwdb)
if ~sum(strcmp(fieldnames(collated_data(1).input_data.sim_settings), 'gamma_b0'))
    do_alternative_gamma = true;
else
    do_alternative_gamma = false;
end

% identify baseline simulation for baseline age and for each age
[base, base_age] = deal(false(length(collated_data),1));
baseline_age = 25;
for sim_no = 1 : length(collated_data)
    if sum(abs(collated_data(sim_no).input_data.sim_settings.variations.params)) == 0
        base_age(sim_no) = true;
        if collated_data(sim_no).input_data.sim_settings.age == baseline_age
            base(sim_no) = true;
        end
    end    
end
subject_ID = 1:length(collated_data);

param_data = nan(length(collated_data),length(params.all));

for subj_no = 1 : length(collated_data)
    
    param_no_counter = 0;
    
    % Insert data for each parameter stored as a variable
    for param_no = 1 : length(params.vars)
        param_no_counter = param_no_counter+1;
        
        curr_param_name = params.vars{param_no};
        
        eval(['curr_param_val = ' curr_param_name '(subj_no);'])
        
        param_data(subj_no,param_no_counter) = curr_param_val;
    end
    
    % Extract this subject's input data
    rel_input_data = collated_data(subj_no).input_data;
    
    % Extract data for each parameter in "sim_settings"
    for param_no = 1 : length(params.sim_settings)
        param_no_counter = param_no_counter+1;
        
        curr_param_name = params.sim_settings{param_no};
        
        % load gamma from a different location
        if do_alternative_gamma & length(curr_param_name)>5 & strcmp(curr_param_name(1:5), 'gamma')
            curr_param_name = strrep(curr_param_name, 'g', 'G');
            eval(['curr_param_val = rel_input_data.sim_settings.fixed.wall.' curr_param_name ';'])
        else % for all other parameters
            eval(['curr_param_val = rel_input_data.sim_settings.' curr_param_name ';'])
        end
        
        param_data(subj_no,param_no_counter) = curr_param_val; clear curr_param_val
    end
    
    % Extract data for each parameter in "network_spec"
    for param_no = 1 : length(params.network_spec)
        curr_param_name = params.network_spec{param_no};
        switch curr_param_name
            case 'k1'
                curr_param_val = rel_input_data.sim_settings.network_spec.k(1);
            case 'k2'
                curr_param_val = rel_input_data.sim_settings.network_spec.k(2);
            case 'k3'
                curr_param_val = rel_input_data.sim_settings.network_spec.k(3);
        end
        param_no_counter = param_no_counter+1;
        param_data(subj_no,param_no_counter) = curr_param_val;
    end
    clear param_no_counter param_no curr_param_name
end

% Change units of some variables
% - Pa to mmHg
rel_param_no = find(strcmp(params.all, 'dbp'));
param_data(:,rel_param_no) = param_data(:,rel_param_no)/133.33;
rel_param_no = find(strcmp(params.all, 'p_out'));
param_data(:,rel_param_no) = param_data(:,rel_param_no)/133.33;

% Change names of some variables
rel_param_no = find(strcmp(params.all, 't_pf'));
params.all{rel_param_no} = 'pft';
rel_param_no = find(strcmp(params.all, 'reg_vol'));
params.all{rel_param_no} = 'rfv';
rel_param_no = find(strcmp(params.all, 'gamma_b0'));
params.all{rel_param_no} = 'b0';
rel_param_no = find(strcmp(params.all, 'gamma_b1'));
params.all{rel_param_no} = 'b1';
rel_param_no = find(strcmp(params.all, 'rho'));
params.all{rel_param_no} = 'density';
rel_param_no = find(strcmp(params.all, 'mu'));
params.all{rel_param_no} = 'viscosity';

% Make header line
header_line = 'Subject Number, ';
for param_no = 2 : length(params.all)
    curr_param = params.all{param_no};
    if length(curr_param)>=4 && strcmp(curr_param(1:4), 'base')
        header_line = [header_line, curr_param, ', '];
    else
        [label, units, abbr, graph_title] = make_param_label(curr_param);
        header_line = [header_line, curr_param, ' [' units, '], '];
    end
end
clear param_no curr_param label units abbr graph_title
header_line = header_line(1:end-2);

% Write header line to file
curr_filename = PATHS.exported_data_model_configs;
fid = fopen(curr_filename,'w');
up.csv.new_line = '\n';
fprintf(fid,[header_line, up.csv.new_line]); clear header_line
fclose(fid); clear fid

% Write parameter data to file
dlmwrite(curr_filename, param_data, '-append');
clear param_data curr_filename

end

function create_pulse_onsets_file(PATHS, collated_data)

fprintf('\n - Creating Onsets File')
%% Setting up

% settings
sites = {'AorticRoot',   'ThorAorta',  'AbdAorta',  'IliacBif', 'Carotid',  'SupTemporal', 'SupMidCerebral', 'Brachial', 'Radial', 'Digital', 'CommonIliac', 'Femoral', 'AntTibial'};
site_domain_no = [1,            18,          39,          41,        15,             87,                 72,         21,       22,       112,             44,        46,           49];
site_dist_prop = [0,             1,           0,           1,       0.5,              1,                  1,       0.75,        1,         1,            0.5,       0.5,            1];
signals = {'P', 'U', 'A', 'PPG'};   % omit 'Q' as it's U.*A

% identify domains which are available
domains = extractfield(collated_data(1).output_data, 'domain_no');
[~,rel_els,~] = intersect(site_domain_no, domains);
site_domain_no = site_domain_no(rel_els);
sites = sites(rel_els);
site_dist_prop = site_dist_prop(rel_els);

%% Creating pulse onsets file (CSV)
up.dataset_name = PATHS.pwdb_filename_prefix;
fs = collated_data(1).output_data(1).fs;

% Cycle through each arterial site
for site_no = 1 : length(sites)
    
    % Identify site-specific info
    curr_site = sites{site_no};
    curr_domain_no = site_domain_no(site_no);
    curr_domain_el = find(domains == curr_domain_no);
    curr_site_dist_prop = site_dist_prop(site_no);
    
    %% Extract PWs at this site from all the simulations
    
    % Cycle through each simulation
    for sim_no = 1 : length(collated_data)
        
        % Identify the relevant data for this simulation at this site
        rel_wave_data = collated_data(sim_no).output_data(curr_domain_el);
        seg_length = collated_data(sim_no).input_data.sim_settings.network_spec.length(curr_domain_no);
        dists = extractfield(rel_wave_data, 'distances');
        [~,rel_dist_el] = min(abs(curr_site_dist_prop*seg_length-dists)); clear dists
        
        % Extract each signal
        for signal_no = 1 : length(signals)
            curr_sig = signals{signal_no};
            
            % calculate pulse onset time
            rel_row = find(domains == 1);
            aortic_root_onset = collated_data(sim_no).output_data(rel_row).start_sample(1);
            if strcmp(curr_sig, 'PPG')
                temp_onset_time = (rel_wave_data.PPG_start_sample(rel_dist_el)-aortic_root_onset)/rel_wave_data.fs;
            else
                temp_onset_time = (rel_wave_data.start_sample(rel_dist_el)-aortic_root_onset)/rel_wave_data.fs;
            end
            % if the time is negative then add on the duration of a cardiac cycle
            if temp_onset_time < 0
                cardiac_cycle_duration = length(rel_wave_data.P(:,1))/rel_wave_data.fs;
                temp_onset_time = temp_onset_time + cardiac_cycle_duration;
            end
            eval(['onset_times.' curr_site '_' curr_sig '(sim_no,1) = temp_onset_time;']);
            clear temp_onset_time aortic_root_onset curr_sig
        end
        clear signal_no rel_wave_data seg_length rel_dist_el
    end
    clear sim_no curr_site_dist_pro
    
    
    
end
clear site_no curr_domain_el curr_domain_no

% Make header line and extract data
col_names = fieldnames(onset_times);
header_line = 'Subject Number, ';
data_mat = nan(length(collated_data),length(col_names)+1);
data_mat(:,1) = 1:length(collated_data);
for col_no = 1 : length(col_names)
    curr_col = col_names{col_no};
    header_line = [header_line, curr_col, ', '];
    eval(['data_mat(:,col_no+1) = onset_times.' curr_col ';']);
end
clear col_no curr_col
header_line = header_line(1:end-2);

% Write header line to file
curr_filename = PATHS.exported_data_onset_times;
fid = fopen(curr_filename,'w');
up.csv.new_line = '\n';
fprintf(fid,[header_line, up.csv.new_line]); clear header_line
fclose(fid); clear fid

% Write parameter data to file
dlmwrite(curr_filename, data_mat, '-append');
clear param_data curr_filename

end

function create_model_geometry_files(PATHS, collated_data)

fprintf('\n - Creating model geometry files: ')

% Setup
params.vars = {'subject_ID'};
params.a = 1;
params.sim_settings = {'age', 'hr', 'sv', 't_pf', 'reg_vol', 'dbp', 'mbp', 'mu', 'alpha', 'p_drop', 'pvc', 'p_out', 'rho', 'lvet', 'pvr'};
params.network_spec = {'seg_no', 'inlet_node', 'outlet_node', 'length', 'inlet_radius', 'outlet_radius', 'c', 'r_sum'};
params.all = [params.network_spec];

% Load data
load(PATHS.peripheral_boundarys)

subject_ID = 1:length(collated_data);

% cycle through each virtual subject
for sim_no = 1 : length(collated_data)
    
    if rem(sim_no,10) == 0
        fprintf([num2str(sim_no), ', '])
    end
    
    % Insert data for each parameter for each segment
    param_no_counter = 0;
    param_data = nan(length(collated_data(1).input_data.sim_settings.network_spec.seg_no),length(params.all));
    for param_no = 1 : length(params.network_spec)
        param_no_counter = param_no_counter+1;
        curr_param_name = params.network_spec{param_no};
        if ~strcmp(curr_param_name, 'c') && ~strcmp(curr_param_name, 'r_sum')
            eval(['curr_param_val = collated_data(sim_no).input_data.sim_settings.network_spec.' curr_param_name ';'])
        else
            % For windkessel boundary characteristics
            eval(['curr_param_val = peripheral_chars.' curr_param_name '(sim_no,:);'])
            curr_param_val(isnan(curr_param_val)) = 0;
        end
        param_data(:,param_no_counter) = curr_param_val; clear curr_param_val
    end
    clear param_no_counter param_no curr_param_name
    
    % Make header line
    header_line = '';
    for param_no = 1 : length(params.all)
        curr_param = params.all{param_no};
        curr_param = strrep(curr_param, 'c', 'peripheral_c');
        curr_param = strrep(curr_param, 'r_sum', 'peripheral_r');
        % Make Header
        header_line = [header_line, curr_param, ', '];
    end
    clear param_no curr_param
    header_line = header_line(1:end-2);
    
    % Write header line to file
    curr_filename = [PATHS.exported_data_geo, PATHS.pwdb_filename_prefix, '_geo_' sprintf('%.4d', sim_no) '.csv'];
    fid = fopen(curr_filename,'w');
    up.csv.new_line = '\n';
    fprintf(fid,[header_line, up.csv.new_line]); clear header_line
    fclose(fid); clear fid
    
    % Write parameter data to file
    dlmwrite(curr_filename, param_data, '-append');
    clear param_data curr_filename
    
end

end

function create_model_variation_file(PATHS, collated_data)

fprintf('\n - Creating model variations file')

% Setup
params = collated_data(1).input_data.sim_settings.variations.param_names;
variations = nan(length(collated_data),length(params));
for subj_no = 1 : length(collated_data)
    variations(subj_no, :) = collated_data(subj_no).input_data.sim_settings.variations.params;
    ages(subj_no,1) = collated_data(subj_no).input_data.sim_settings.age;
end
rel_params = sum(abs(variations))>0;
params = params(rel_params);
variations = variations(:,rel_params); clear rel_params subj_no

subjs = 1:length(collated_data); subjs = subjs(:);
data_mat = [subjs, ages, variations]; clear variations ages subjs

% Make header line
header_line = 'Subject Number, Age, ';
for param_no = 1 : length(params)
    curr_param = params{param_no};
    % Make Header
    header_line = [header_line, curr_param, ', '];
end
clear param_no curr_param params
header_line = header_line(1:end-2);

% Rename some variables
header_line = strrep(header_line, 'reg_vol', 'RFV');
header_line = strrep(header_line, 't_pf', 'PFT');
header_line = upper(header_line);

% Write header line to file
curr_filename = [PATHS.exported_data, PATHS.pwdb_filename_prefix, '_model_variations.csv'];
fid = fopen(curr_filename,'w');
up.csv.new_line = '\n';
fprintf(fid,[header_line, up.csv.new_line]); clear header_line
fclose(fid); clear fid

% Write parameter data to file
dlmwrite(curr_filename, data_mat, '-append');
clear param_data curr_filename

end

function v = find_pw_point_time(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig)

% identify distance along segment
desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
[~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));

% identify the sample corresponding to this fiducial point
rel_row = find(strcmp(pulse_wave_inds(1).fid_pt_names, pwa_pt));
eval(['rel_el = pulse_wave_inds(subj_no).' pwa_sig '_pwa(domain_el).fid_pts(rel_row,rel_dist_el);'])

% calculate the time of this fiducial point
fs = collated_data(subj_no).output_data.fs;
v = (rel_el-1)/fs;

end

function v = find_pw_point_value(collated_data, pulse_wave_inds, subj_no, site_domain_no, site_no, site_dist_prop, domain_el, pwa_pt, pwa_sig)

% identify distance along segment
desired_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(site_domain_no(site_no))*site_dist_prop(site_no);
[~, rel_dist_el] = min(abs(collated_data(subj_no).output_data(domain_el).distances-desired_length));

% identify the sample corresponding to this fiducial point
rel_row = find(strcmp(pulse_wave_inds(1).fid_pt_names, pwa_pt));
eval(['rel_el = pulse_wave_inds(subj_no).' pwa_sig '_pwa(domain_el).fid_pts(rel_row,rel_dist_el);'])

% Extract the value of the signal at this point
if ~isnan(rel_el)
    switch pwa_sig
        case 'PPG'
            v = collated_data(subj_no).output_data(domain_el).PPG(rel_el,rel_dist_el);
        case 'P'
            v = collated_data(subj_no).output_data(domain_el).P(rel_el,rel_dist_el)/133.33;
    end
else
    v = nan;
end

end

function create_haemod_param_file(PATHS)

fprintf('\n - Creating haemodynamic parameters file')

% Load data
load(PATHS.haemodynamic_params)

% setup
params = {'subj_no', 'age', 'HR', 'SV', 'CO', 'LVET', 'dPdt', 'PFT', 'RFV', 'SBP_a', 'DBP_a', 'MBP_a', 'PP_a', 'SBP_b', 'DBP_b', 'MBP_b', 'PP_b', 'PP_amp', 'AP', 'AIx', 'Tr', 'PWV_a', 'PWV_cf', 'PWV_br', 'PWV_fa', 'dia_asc_a', 'dia_desc_thor_a', 'dia_abd_a', 'dia_car', 'len_prox_a', 'MBP_drop_finger', 'MBP_drop_ankle', 'svr'};

% Make header line
header_line = 'Subject Number, ';
for param_no = 2 : length(params)
    curr_param = params{param_no};
    [label, units, abbr, graph_title] = make_param_label(curr_param);
    abbr = strrep(abbr, '{', '');
    abbr = strrep(abbr, '}', '');
    header_line = [header_line, abbr, ' [' units, '], '];
end
clear param_no curr_param label units abbr graph_title
header_line = header_line(1:end-2);

% Write header line to file
curr_filename = [PATHS.exported_data, PATHS.pwdb_filename_prefix, '_haemod_params.csv'];
fid = fopen(curr_filename,'w');
up.csv.new_line = '\n';
fprintf(fid,[header_line, up.csv.new_line]); clear header_line
fclose(fid); clear fid

% Extract parameter data
params_mat = 1:length(haemodynamic_params); params_mat = params_mat(:);
for param_no = 2 : length(params)
    curr_param = params{param_no};
    curr_param = strrep(curr_param, 'AP', 'AP_a');
    curr_param = strrep(curr_param, 'AIx', 'AI_a');
    curr_param = strrep(curr_param, 'Tr', 'Tr_a');
    eval(['params_mat(:,end+1) = extractfield(haemodynamic_params, ''' curr_param ''');'])
end

% Write parameter data to file
dlmwrite(curr_filename, params_mat, '-append'); clear params_mat curr_filename

end

function data = assess_phys_plausibility(data)

%% Extract data required to assess physiological plausibility

p_data = proc_data(data);

%% Assess physiological plausibility

% Compare simulated PW characteristics to literature ones
[implaus_sims, all_implaus_sims, res] = compare_simulated_vs_literature(p_data);

%% Store results
data.plausibility.implaus_sims = implaus_sims;
data.plausibility.all_implaus_sims = all_implaus_sims;
data.plausibility.plausibility_log = false(length(p_data.age),1);
data.plausibility.plausibility_log(setxor(1:length(p_data.age), all_implaus_sims.els)) = true;

end

function p_data = proc_data(data)
% Extract pulse wave indices
p_data.SBP_a = data.pw_inds.AorticRoot_SBP_V;
p_data.PP_a = data.pw_inds.AorticRoot_PP;
p_data.SBP_b = data.pw_inds.Brachial_SBP_V;
p_data.PP_b = data.pw_inds.Brachial_PP;
p_data.DBP_b = data.pw_inds.Brachial_DBP;
p_data.MBP_b = data.pw_inds.Brachial_MBP;

% Extract haemodynamic parameters
p_data.PP_amp = extractfield(data.haemods, 'PP_amp');
p_data.Tr_a = extractfield(data.haemods, 'Tr_a');
p_data.AI_a = extractfield(data.haemods, 'AI_a');
p_data.AI_c = extractfield(data.haemods, 'AI_c');
p_data.AP_a = extractfield(data.haemods, 'AP_a');
p_data.AP_c = extractfield(data.haemods, 'AP_c');
p_data.HR = extractfield(data.haemods, 'HR');
p_data.SV = extractfield(data.haemods, 'SV');
p_data.CO = extractfield(data.haemods, 'CO');
p_data.LVET = extractfield(data.haemods, 'LVET');
p_data.PWV = extractfield(data.haemods, 'PWV_a');
p_data.DIA = extractfield(data.haemods, 'dia_asc_a');
if sum(strcmp(fieldnames(data.haemods), 'c'))
    p_data.C = extractfield(data.haemods, 'c');
end

% Extract configuration
p_data.age = data.config.age;
p_data.variations = data.config.variations;
end

function [implaus_sims, all_implaus_sims, res] = compare_simulated_vs_literature(p_data)

fprintf('\n - Compare simulated vs literature characteristics')

% Simulated characteristics
params = {'AP_a', 'SBP_b', 'DBP_b', 'PP_b', 'MBP_b', 'SBP_a', 'PP_a', 'PP_amp', 'AI_a', 'Tr_a'};
params = {'SBP_b', 'DBP_b', 'PP_b', 'MBP_b', 'SBP_a', 'PP_a', 'PP_amp'};

%% Extract literature values of each characteristic
req_ages = unique(p_data.age(:));
all_implaus_sims.els = [];
all_implaus_sims.params = zeros(length(p_data.age),length(params));
for param_no = 1 : length(params)
    
    curr_param = params{param_no};
    
    % Extract literature values of this parameter
    literature_data.age = 15:10:85;
    switch curr_param
        case 'SBP_b' % taken from McEniery2005 (male data, Table 1)
            literature_data.male.mean = [123, 124, 123, 125, 125, 126, 127, 130];
            literature_data.male.sd = [10,10,9,9,9,9,9,8];
            literature_data.male.n = [172,178,183,258,429,430,280,39];
            literature_data.female.mean = [113,115,115,118,122,126,127,128];
            literature_data.female.sd = [10,10,12,11,11,10,10,10];
            literature_data.female.n = [133,101,165,301,495,509,290,38];
        case 'DBP_b' % taken from McEniery2005 (male data, Table 1)
            literature_data.male.mean = [73,75,77,79,79,78,76,75];
            literature_data.male.sd = [8,10,9,9,9,9,9,8];
            literature_data.male.n = [172,178,183,258,429,430,280,39];
            literature_data.female.mean = [72,73,74,75,75,74,72,70];
            literature_data.female.sd = [8,8,9,8,7,7,8,9];
            literature_data.female.n = [133,101,165,301,495,509,290,38];
        case 'PP_b' % taken from McEniery2005 (male data, Table 1)
            literature_data.male.mean = [50,49,47,46,46,49,51,55];
            literature_data.male.sd = [9,9,8,7,8,8,8,9];
            literature_data.male.n = [172,178,183,258,429,430,280,39];
            literature_data.female.mean = [41,43,41,43,46,51,54,57];
            literature_data.female.sd = [8,7,9,9,9,8,9,11];
            literature_data.female.n = [133,101,165,301,495,509,290,38];
        case 'MBP_b' % taken from McEniery2005 (male data, Table 1)
            literature_data.male.mean = [88,89,92,95,95,94,93,92];
            literature_data.male.sd = [8,8,8,7,7,7,7,8];
            literature_data.male.n = [172,178,183,258,429,430,280,39];
            literature_data.female.mean = [86,86,88,90,93,93,92,90];
            literature_data.female.sd = [8,8,9,9,8,8,8,8];
            literature_data.female.n = [133,101,165,301,495,509,290,38];
        case 'SBP_a' % taken from McEniery2005 (male data, Table 1)
            literature_data.male.mean = [103,105,109,113,115,117,118,120];
            literature_data.male.sd = [8,8,9,9,9,9,9,8];
            literature_data.male.n = [172,178,183,258,429,430,280,39];
            literature_data.female.mean = [98,101,105,109,115,118,119,120];
            literature_data.female.sd = [9,9,11,11,11,10,9,11];
            literature_data.female.n = [133,101,165,301,495,509,290,38];
        case 'PP_a' % taken from McEniery2005 (male data, Table 1)
            literature_data.male.mean = [29,30,31,34,35,39,42,45];
            literature_data.male.sd = [5,6,6,6,7,7,7,9];
            literature_data.male.n = [172,178,183,258,429,430,280,39];
            literature_data.female.mean = [25,27,30,33,38,43,56,49];
            literature_data.female.sd = [6,7,8,8,8,8,8,12];
            literature_data.female.n = [133,101,165,301,495,509,290,38];
        case 'PP_amp' % taken from McEniery2005 (male data, Table 1)
            literature_data.male.mean = [1.72, 1.7 , 1.50, 1.39, 1.33, 1.26, 1.24, 1.25];
            literature_data.male.sd =   [0.11, 0.14, 0.18, 0.15, 0.16, 0.13, 0.12, 0.15];
            literature_data.male.n = [172,178,183,258,429,430,280,39];
            literature_data.female.mean = [1.67,1.59,1.41,1.29,1.22,1.21,1.19,1.18];
            literature_data.female.sd = [0.15,0.2,0.18,0.15,0.11,0.10,0.10,0.11];
            literature_data.female.n = [133,101,165,301,495,509,290,38];
        case 'AP_a' % taken from McEniery2005 (male data, Table 1)
            literature_data.male.mean = [-1,1,4,7,9,11,13,14];
            literature_data.male.sd = [3,4,5,4,5,5,5,5];
            literature_data.male.n = [172,178,183,258,429,430,280,39];
            literature_data.female.mean = [1,3,6,10,13,15,16,18];
            literature_data.female.sd = [3,4,5,5,5,5,5,7];
            literature_data.female.n = [133,101,165,301,495,509,290,38];
        case 'AI_a' % taken from McEniery2005 (male data, Table 1)
            literature_data.male.mean = [-2,2,12,19,24,28,30,30];
            literature_data.male.sd = [8,11,13,10,10,9,9,10];
            literature_data.male.n = [172,178,183,258,429,430,280,39];
            literature_data.female.mean = [5,9,20,28,33,34,35,37];
            literature_data.female.sd = [10,14,12,10,9,9,9,10];
            literature_data.female.n = [133,101,165,301,495,509,290,38];
        case 'Tr' % taken from McEniery2005 (male data, Table 1)
            literature_data.male.mean = [150,154,151,148,143,141,136,133];
            literature_data.male.sd = [17,21,21,16,15,12,12,16];
            literature_data.male.n = [172,178,183,258,429,430,280,39];
            literature_data.female.mean = [145,143,140,136,133,131,129,125];
            literature_data.female.sd = [16,13,16,14,15,14,12,12];
            literature_data.female.n = [133,101,165,301,495,509,290,38];
    end
    
    [literature_data.overall_mean, literature_data.overall_sd] = calc_overall_stats(literature_data);
    
    interp_literature_data.mean = interp1(literature_data.age, literature_data.overall_mean, req_ages);
    interp_literature_data.sd = interp1(literature_data.age, literature_data.overall_sd, req_ages);
    interp_literature_data.age = req_ages;
    
    % store literature values
    eval(['literature.' curr_param ' = interp_literature_data;']);
    
    % Extract simulated values
    simulated_data.age = unique(p_data.age(:)');
    param_implaus_els_high = [];
    param_implaus_els_low = [];
    
    for age_no = 1 : length(interp_literature_data.age)
        curr_age = p_data.age(age_no);
        
        % 99% ranges
        rel_lit_el = find(interp_literature_data.age == curr_age);
        max_val = interp_literature_data.mean(rel_lit_el)+ 2.575*interp_literature_data.sd(rel_lit_el);
        min_val = interp_literature_data.mean(rel_lit_el)- 2.575*interp_literature_data.sd(rel_lit_el);
        clear rel_lit_el
        
        % Relevant simulated values
        rel_els = find(p_data.age == curr_age);
        eval(['param_data = p_data.' curr_param ';']);
        rel_vals = param_data(rel_els); clear param_data curr_age
        
        % simulations outside plausible range
        if ~exist('implaus_sims', 'var') || ~sum(strcmp(fieldnames(implaus_sims), curr_param))
            eval(['implaus_sims.' curr_param '.low_els = [];']);
            eval(['implaus_sims.' curr_param '.low_vars = [];']);
            eval(['implaus_sims.' curr_param '.high_els = [];']);
            eval(['implaus_sims.' curr_param '.high_vars = [];']);
        end
        
        % simulations in which the value of this parameter is lower than expected
        curr_implaus_els = rel_els(rel_vals < min_val);
        eval(['implaus_sims.' curr_param '.low_els = [implaus_sims.' curr_param '.low_els; curr_implaus_els];']);
        curr_implaus_vars = p_data.variations.params(curr_implaus_els, [1:3,7:9]);
        eval(['implaus_sims.' curr_param '.low_vars = [implaus_sims.' curr_param '.low_vars; curr_implaus_vars];']);
        param_implaus_els_low = [param_implaus_els_low; curr_implaus_els];
        all_implaus_sims.params(curr_implaus_els,param_no) = -1;
        
        % simulations in which the value of this parameter is higher than expected
        curr_implaus_els = rel_els(rel_vals > max_val);
        eval(['implaus_sims.' curr_param '.high_els = [implaus_sims.' curr_param '.high_els; curr_implaus_els];']);
        curr_implaus_vars = p_data.variations.params(curr_implaus_els, [1:3,7:9]);
        eval(['implaus_sims.' curr_param '.high_vars = [implaus_sims.' curr_param '.high_vars; curr_implaus_vars];']);
        param_implaus_els_high = [param_implaus_els_high; curr_implaus_els];
        all_implaus_sims.params(curr_implaus_els,param_no) = 1;
                
        clear rel_vals max_val min_val  curr_implaus_els curr_age rel_els
    end
    fprintf(['\n      ' curr_param ': ' num2str(length(param_implaus_els_low)) ' implausibly low; ']);
    fprintf([num2str(length(param_implaus_els_high)) ' implausibly high.']);
    all_implaus_sims.els = [all_implaus_sims.els; param_implaus_els_low; param_implaus_els_high];
            
    clear age_no
    
    clear literature_data interp_literature_data
end
implaus_sims.var_names = p_data.variations.param_names([1:3,7:9]);

all_implaus_sims.els = unique(all_implaus_sims.els);
all_implaus_sims.vars = p_data.variations.params(all_implaus_sims.els, [1:3,7:9]);
all_implaus_sims.var_names = p_data.variations.param_names([1:3,7:9]);
all_implaus_sims.age = p_data.age(all_implaus_sims.els);
all_implaus_sims.params = all_implaus_sims.params(all_implaus_sims.els,:);

%% Calculate Results

% Summary of implausible subjects
res.no_implaus_sims = length(all_implaus_sims.els);
res.no_implaus_sims_due_to_PP = sum(sum(abs(all_implaus_sims.params(:,[3,6])),2)~=0);
high_PP_els = sum(all_implaus_sims.params(:,[3,6]),2)>0;
res.no_implaus_sims_due_to_PP_high = sum(high_PP_els);
low_PP_els = sum(all_implaus_sims.params(:,[3,6]),2)<0;
res.no_implaus_sims_due_to_PP_low = sum(low_PP_els);
rel_els = sum(abs(all_implaus_sims.params(:,[3,6])),2)~=0;
other_els = logical(1-rel_els);
high_PPamp_els = sum(all_implaus_sims.params(other_els,[7]),2)>0;
res.no_implaus_sims_due_to_high_PP_amp = sum(high_PPamp_els);
res.prop_implaus_25 = mean(all_implaus_sims.age==25);
res.prop_implaus_75 = mean(all_implaus_sims.age==75);
rel_col = strcmp(implaus_sims.var_names, 'pwv');
res.prop_implaus_high_PWV = mean(all_implaus_sims.vars(:,rel_col)>0);
res.prop_implaus_low_PWV = mean(all_implaus_sims.vars(:,rel_col)<0);

do_check = 0;
if do_check
    % High PP
    rel_data = all_implaus_sims.vars(high_PP_els,:); n = 1; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=2; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=3; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=4; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=5; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=6; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n));
    % Low PP
    rel_data = all_implaus_sims.vars(low_PP_els,:); n = 1; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=2; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=3; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=4; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=5; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=6; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n));
    % Normal PP and High PPamp
    rel_col = strcmp(params, 'PP_amp');
    temp_els = sum(all_implaus_sims.params(:,[3,6]),2)==0 & all_implaus_sims.params(:,rel_col) == 1;
    rel_data = all_implaus_sims.vars(temp_els,:); n = 1; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=2; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=3; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=4; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=5; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=6; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n));
    % Normal PP, High PPamp and low PWV
    rel_ppamp_col = strcmp(params, 'PP_amp');
    rel_pwv_col = strcmp(implaus_sims.var_names, 'pwv');
    temp_els = sum(all_implaus_sims.params(:,[3,6]),2)==0 & all_implaus_sims.params(:,rel_ppamp_col) == 1 & all_implaus_sims.vars(:,rel_pwv_col) == -1;
    rel_data = all_implaus_sims.vars(temp_els,:); n = 1; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=2; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=3; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=4; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=5; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n)); n=6; subplot(3,2,n), histogram(rel_data(:,n)), xlabel(implaus_sims.var_names(n));
end

end

function [overall_mean, overall_sd] = calc_overall_stats(literature_data)

for age_no = 1 : length(literature_data.age)
    
    % From https://www.statstodo.com/CombineMeansSDs_Pgm.php
    
    curr_data = literature_data.female;
    sum_x_f = curr_data.mean(age_no)*curr_data.n(age_no);
    sum_xsq_f = (curr_data.sd(age_no)^2)*(curr_data.n(age_no)-1)+(sum_x_f^2/curr_data.n(age_no));
    curr_data = literature_data.male;
    sum_x_m = curr_data.mean(age_no)*curr_data.n(age_no);
    sum_xsq_m = (curr_data.sd(age_no)^2)*(curr_data.n(age_no)-1)+(sum_x_m^2/curr_data.n(age_no));
    tn = literature_data.male.n(age_no) + literature_data.female.n(age_no);
    tx = sum_x_f + sum_x_m;
    txx = sum_xsq_m + sum_xsq_f;
    overall_mean(age_no) = tx/tn;
    overall_sd(age_no) = sqrt((txx-(tx^2)/tn)/(tn-1));
    
    clear curr_data sum_x_f sum_xsq_f sum_x_m sum_xsq_m tn tx txx
    
end


end

function create_mat_file_w_all_data(PATHS, collated_data, pw_inds)

fprintf('\n - Collating data in single Matlab file')

% Load peripheral Windkessel boundary conditions
load(PATHS.peripheral_boundarys)

% Load haemodynamic parameters
load(PATHS.haemodynamic_params)

% Load system characteristics
if exist(PATHS.system_chars, 'file')
    load(PATHS.system_chars)
end

% settings
sites = {'AorticRoot',   'ThorAorta',  'AbdAorta',  'IliacBif', 'Carotid',  'SupTemporal', 'SupMidCerebral', 'Brachial', 'Radial', 'Digital', 'CommonIliac', 'Femoral', 'AntTibial'};
site_domain_no = [1,            18,          39,          41,        15,             87,                 72,         21,       22,       112,             44,        46,           49];
site_dist_prop = [0,             1,           0,           1,       0.5,              1,                  1,       0.75,        1,         1,            0.5,       0.5,            1];
signals = {'P', 'U', 'A', 'PPG'};

% identify domains which are available
domains = extractfield(collated_data(1).output_data, 'domain_no');
[~,rel_els,~] = intersect(site_domain_no, domains);
site_domain_no = site_domain_no(rel_els);
sites = sites(rel_els);
site_dist_prop = site_dist_prop(rel_els);

% Check names of gamma variables (this is required for the ppwdb)
if ~sum(strcmp(fieldnames(collated_data(1).input_data.sim_settings), 'gamma_b0'))
    do_alternative_gamma = true;
else
    do_alternative_gamma = false;
end

% waves along arterial paths
do_path = 1;
if do_path == 1
    % choose arterial paths
    if ~do_alternative_gamma
        path_names = {'aorta_finger', 'aorta_foot', 'aorta_brain', 'aorta_r_subclavian'};
    else
        path_names = {'aorta_finger', 'aorta_foot'};
    end
    for s = 1 : length(path_names)
        arterial_paths(s).name = path_names{s};
    end
    clear path_names s
    
    % add segment data for each path
    for path_name_no = 1 : length(arterial_paths)
        curr_path_name = arterial_paths(path_name_no).name;
        switch curr_path_name
            case 'aorta_finger'
                arterial_paths(path_name_no).seg_names = {'asc_aorta', 'aortic_arch1', 'aortic_arch2','subclavian','brachial','radial','palmar','digital'};
                arterial_paths(path_name_no).doms = [1,2,14,19,21,22,108,112];
            case 'aorta_foot'
                arterial_paths(path_name_no).seg_names = {'asc_aorta', 'aortic_arch1', 'aortic_arch2','desc_thor_aorta1','desc_thor_aorta2','abd_aorta1','abd_aorta2','abd_aorta3','abd_aorta4','abd_aorta5','common_iliac','external_iliac','femoral','ant_tibial'};
                arterial_paths(path_name_no).doms = [1,2,14,18,27,28,35,37,39,41,42,44,46,49];
            case 'aorta_brain'
                arterial_paths(path_name_no).seg_names = {'asc_aorta', 'aortic_arch1', 'carotid','int_carotid1','int_carotid2','int_carotid_dist1','int_carotid_dist2','mid_cerebral','sup_mid_cerebral'};
                arterial_paths(path_name_no).doms = [1,2,15,16,79,65,96,71,72];
            case 'aorta_r_subclavian'
                arterial_paths(path_name_no).seg_names = {'asc_aorta', 'brachiocephalic', 'subclavian'};
                arterial_paths(path_name_no).doms = [1,3,4];
        end
    end
end

data.config.baseline_sim_for_all = false(length(collated_data),1);
data.config.baseline_sim_for_age = false(length(collated_data),1);
    
% Cycle through each simulation
for sim_no = 1 : length(collated_data)
    
    %% Extract model input parameters for this simulation
    input_vars = {'hr', 'sv', 'lvet', 't_pf', 'reg_vol', 'dbp', 'mbp', 'pvr', 'pvc', 'p_out', 'len', 'dia', 'pwv', 'age'};
    for input_var_no = 1 : length(input_vars)
        
        % rename variable
        curr_var_name = strrep(input_vars{input_var_no}, 't_pf', 'pft');
        curr_var_name = strrep(curr_var_name, 'reg_vol', 'rfv');
        
        try
            eval(['data.config.' curr_var_name '(sim_no,1) = collated_data(sim_no).input_data.sim_settings.' input_vars{input_var_no} ';']);
        catch
            rel_col = find(strcmp(collated_data(sim_no).input_data.sim_settings.variations.param_names,input_vars{input_var_no}));
            eval(['data.config.' curr_var_name '(sim_no,1) = collated_data(sim_no).input_data.sim_settings.variations.params(rel_col);']);
        end
        if ~strcmp(input_vars{input_var_no}, 'age')
            if strcmp(input_vars{input_var_no}, 'pvr')
                rel_col = find(strcmp(collated_data(sim_no).input_data.sim_settings.variations.param_names,'mbp'));
            else
                rel_col = find(strcmp(collated_data(sim_no).input_data.sim_settings.variations.param_names,input_vars{input_var_no}));
            end
            eval(['data.config.' curr_var_name '_SD(sim_no,1) = collated_data(sim_no).input_data.sim_settings.variations.params(rel_col);']);
        end
    end
    clear input_var_no input_vars rel_col
    
    % Note down whether this is a baseline simulation
    if sum(collated_data(sim_no).input_data.sim_settings.variations.params ~= 0) == 0
        data.config.baseline_sim_for_age(sim_no) = true;
        if collated_data(sim_no).input_data.sim_settings.age == 25
            data.config.baseline_sim_for_all(sim_no) = true;
        end
    end
    
    %% Extract Simulated Waves
    
    % Add in Waves along paths
    if do_path
        % cycle through each path
        for path_name_no = 1 : length(arterial_paths)
            curr_path_name = arterial_paths(path_name_no).name;
            
            % Extract data for this path
            running_dist = 0; counter_no = 0;
            domains = extractfield(collated_data(sim_no).output_data, 'domain_no');
            for site_no = 1 : length(arterial_paths(path_name_no).seg_names)
                curr_domain_no = arterial_paths(path_name_no).doms(site_no);
                curr_domain_el = find(domains == curr_domain_no);
                rel_wave_data = collated_data(sim_no).output_data(curr_domain_el);
                seg_length = collated_data(sim_no).input_data.sim_settings.network_spec.length(curr_domain_no);
                
                for dist_el = 1:length(rel_wave_data.distances)
                    curr_dist = running_dist + rel_wave_data.distances(dist_el);
%                     if site_no > 1 && rel_wave_data.distances(dist_el) == 0
%                         continue
%                     end
                    counter_no = counter_no+1;
                    eval(['curr_seg_name = ''path' num2str(counter_no) '_' arterial_paths(path_name_no).seg_names{site_no} ''';']);
                    signals2 = {'P', 'U', 'A'};
                    for signal_no = 1 : length(signals2)
                        if strcmp(signals2{signal_no}, 'P')
                            factor = 133.33;
                        else
                            factor = 1;
                        end
                        curr_sig = signals2{signal_no};
                        
                        eval(['data.path_waves.' curr_path_name '(sim_no).' curr_sig '{counter_no} = rel_wave_data.' curr_sig '(:,dist_el)/factor;'])
                        clear factor
                    end
                    eval(['data.path_waves.' curr_path_name '(sim_no).dist(counter_no,1) = curr_dist;']);
                    eval(['data.path_waves.' curr_path_name '(sim_no).artery{counter_no,1} = arterial_paths(path_name_no).seg_names{site_no};']);
                    eval(['data.path_waves.' curr_path_name '(sim_no).segment_no{counter_no,1} = arterial_paths(path_name_no).doms(site_no);']);
                    eval(['data.path_waves.' curr_path_name '(sim_no).artery_dist(counter_no,1) = rel_wave_data.distances(dist_el);']);
                    
                    % calculate pulse onset time
                    rel_row = find(domains == 1);
                    aortic_root_onset = collated_data(sim_no).output_data(rel_row).start_sample(1); clear rel_row
                    if strcmp(curr_sig, 'PPG')
                        temp_onset_time = (rel_wave_data.PPG_start_sample(dist_el)-aortic_root_onset)/rel_wave_data.fs;
                    else
                        temp_onset_time = (rel_wave_data.start_sample(dist_el)-aortic_root_onset)/rel_wave_data.fs;
                    end
                    % if the time is negative then add on the duration of a cardiac cycle
                    if temp_onset_time < 0
                        cardiac_cycle_duration = length(rel_wave_data.P(:,1))/rel_wave_data.fs;
                        temp_onset_time = temp_onset_time + cardiac_cycle_duration;
                    end
                    % store pulse onset time
                    eval(['data.path_waves.' curr_path_name '(sim_no).onset_time(counter_no,1)  = temp_onset_time;']);
                    clear aortic_root_onset
                    
                end
                running_dist = running_dist + seg_length;
                clear signal_no curr_sig rel_wave_data seg_length curr_domain_el curr_domain_no
            end
            clear site_no
        end
    end
    
    % Add in Waves at each selected site
    domains = extractfield(collated_data(sim_no).output_data, 'domain_no');
    for site_no = 1 : length(sites)
        curr_domain_no = site_domain_no(site_no);
        curr_domain_el = find(domains == curr_domain_no);
        rel_wave_data = collated_data(sim_no).output_data(curr_domain_el);
        seg_length = collated_data(sim_no).input_data.sim_settings.network_spec.length(curr_domain_no);
        [~, rel_dist_el] = min(abs(rel_wave_data.distances - (site_dist_prop(site_no)*seg_length)));
        for signal_no = 1 : length(signals)
            if strcmp(signals{signal_no}, 'P')
                factor = 133.33;
            else
                factor = 1;
            end
            curr_sig = signals{signal_no};
            
            eval(['data.waves.' curr_sig '_' sites{site_no} '{sim_no} = rel_wave_data.' curr_sig '(:,rel_dist_el)/factor;'])
            clear factor
            
            % calculate pulse onset time
            rel_row = find(domains == 1);
            aortic_root_onset = collated_data(sim_no).output_data(rel_row).start_sample(1); clear rel_row
            if strcmp(curr_sig, 'PPG')
                temp_onset_time = (rel_wave_data.PPG_start_sample(rel_dist_el)-aortic_root_onset)/rel_wave_data.fs;
            else
                temp_onset_time = (rel_wave_data.start_sample(rel_dist_el)-aortic_root_onset)/rel_wave_data.fs;
            end
            % if the time is negative then add on the duration of a cardiac cycle
            if temp_onset_time < 0
                cardiac_cycle_duration = length(rel_wave_data.P(:,1))/rel_wave_data.fs;
                temp_onset_time = temp_onset_time + cardiac_cycle_duration;
            end
            % store pulse onset time
            eval(['data.waves.onset_times.' curr_sig '_' sites{site_no} '(sim_no,1) = temp_onset_time;'])
            clear aortic_root_onset
        end
        clear signal_no curr_sig rel_wave_data seg_length curr_domain_el curr_domain_no
    end
    clear site_no
    
    %% Extract variation for this simulation
    a = collated_data(sim_no).input_data.sim_settings.variations.param_names;
    req_order = {'hr', 'sv', 'lvet', 't_pf', 'reg_vol', 'len', 'dia', 'pwv', 'mbp', 'pvc'};
    curr_params = a;
    req_els = nan(length(curr_params),1);
    counter = length(intersect(a,req_order))+1;
    for param_no = 1 : length(curr_params)
        curr_param = curr_params{param_no};
        if sum(strcmp(req_order, curr_param))
            position = find(strcmp(req_order, curr_param));
            req_els(position) = param_no;
        else
            position = counter;
            req_els(position) = param_no;
            counter = counter+1;
        end
    end
    clear req_order var_no
    
    data.config.variations.params(sim_no,:) = collated_data(sim_no).input_data.sim_settings.variations.params(req_els);
    current_names = collated_data(sim_no).input_data.sim_settings.variations.param_names(req_els);
    if sum(strcmp(fieldnames(data.config.variations), 'param_names')) & ~isequal(current_names, data.config.variations.param_names)
        error('check this')
    elseif sim_no == 1
        data.config.variations.param_names = current_names;
    end
    
    %% Extract desired PWVs
    vars = {'pwv_aorta', 'pwv_leg', 'pwv_arm'};
    mod_names = {'pwv_cf', 'pwv_fa', 'pwv_br'};
    for var_no = 1 : length(vars)
        eval(['data.config.desired_chars.' mod_names{var_no} '(sim_no,1) = collated_data(sim_no).input_data.sim_settings.desired_' vars{var_no} ';']);
    end
    clear var_no vars
    
    %% Extract desired aortic geometry
    % Diameters using values measured from the simulated waves
    domain_nos = extractfield(collated_data(1).output_data, 'domain_no');
    % - asc aorta
    rel_domain_nos = find(collated_data(sim_no).input_data.sim_settings.network_spec.asc_aorta);
    [~,domain_els,~] = intersect(domain_nos, rel_domain_nos);
    for s = 1 : length(domain_els)
        curr_domain_el = domain_els(s);
        ave_areas(s,1) = mean([max(collated_data(sim_no).output_data(curr_domain_el).A(:,1)), max(collated_data(sim_no).output_data(curr_domain_el).A(:,end))]);
    end
    lengths = collated_data(sim_no).input_data.sim_settings.network_spec.length(rel_domain_nos);
    data.config.desired_chars.dia_asc_a(sim_no,1) = 1000*2*(sum(lengths.*sqrt(ave_areas./pi)))/sum(lengths);
    clear ave_areas lengths
    % - desc thor aorta
    rel_domain_nos = find(collated_data(sim_no).input_data.sim_settings.network_spec.desc_thor_aorta);
    [~,domain_els,~] = intersect(domain_nos, rel_domain_nos);
    for s = 1 : length(domain_els)
        curr_domain_el = domain_els(s);
        ave_areas(s,1) = mean([max(collated_data(sim_no).output_data(curr_domain_el).A(:,1)), max(collated_data(sim_no).output_data(curr_domain_el).A(:,end))]);
    end
    lengths = collated_data(sim_no).input_data.sim_settings.network_spec.length(rel_domain_nos);
    data.config.desired_chars.dia_desc_thor_a(sim_no,1) = 1000*2*(sum(lengths.*sqrt(ave_areas./pi)))/sum(lengths);
    clear ave_areas lengths
    % - abd aorta
    rel_domain_nos = find(collated_data(sim_no).input_data.sim_settings.network_spec.abd_aorta);
    [~,domain_els,~] = intersect(domain_nos, rel_domain_nos);
    for s = 1 : length(domain_els)
        curr_domain_el = domain_els(s);
        ave_areas(s,1) = mean([max(collated_data(sim_no).output_data(curr_domain_el).A(:,1)), max(collated_data(sim_no).output_data(curr_domain_el).A(:,end))]);
    end
    lengths = collated_data(sim_no).input_data.sim_settings.network_spec.length(rel_domain_nos);
    data.config.desired_chars.dia_abd_a(sim_no,1) = 1000*2*(sum(lengths.*sqrt(ave_areas./pi)))/sum(lengths);
    clear ave_areas lengths
    % - carotid
    rel_domain_nos = find(collated_data(sim_no).input_data.sim_settings.network_spec.both_carotid);
    [~,domain_els,~] = intersect(domain_nos, rel_domain_nos);
    for s = 1 : length(domain_els)
        curr_domain_el = domain_els(s);
        ave_areas(s,1) = mean([mean(collated_data(sim_no).output_data(curr_domain_el).A(:,1)), mean(collated_data(sim_no).output_data(curr_domain_el).A(:,end))]);
    end
    lengths = collated_data(sim_no).input_data.sim_settings.network_spec.length(rel_domain_nos);
    data.config.desired_chars.dia_car(sim_no,1) = 1000*2*(sum(lengths.*sqrt(ave_areas./pi)))/sum(lengths);
    clear ave_areas lengths    
    
    rel_els = collated_data(sim_no).input_data.sim_settings.network_spec.proximal_aorta;
    data.config.desired_chars.len_prox_a(sim_no,1) = 1000*sum( collated_data(sim_no).input_data.sim_settings.network_spec.length(rel_els) );
    
    %% Store arterial network properties
    data.config.network.segment_name = collated_data(sim_no).input_data.sim_settings.network_spec.segment_name;
    data.config.network.inlet_node(sim_no,:) = collated_data(sim_no).input_data.sim_settings.network_spec.inlet_node;
    data.config.network.outlet_node(sim_no,:) = collated_data(sim_no).input_data.sim_settings.network_spec.outlet_node;
    data.config.network.length(sim_no,:) = collated_data(sim_no).input_data.sim_settings.network_spec.length;
    data.config.network.inlet_radius(sim_no,:) = collated_data(sim_no).input_data.sim_settings.network_spec.inlet_radius;
    data.config.network.outlet_radius(sim_no,:) = collated_data(sim_no).input_data.sim_settings.network_spec.outlet_radius;
    data.config.constants.k(sim_no,:) = collated_data(sim_no).input_data.sim_settings.network_spec.k;
    
    %% Store other simulation properties
    data.config.constants.rho(sim_no,:) = collated_data(sim_no).input_data.sim_settings.rho;
    if ~do_alternative_gamma
        data.config.constants.gamma_b0(sim_no,1) = collated_data(sim_no).input_data.sim_settings.gamma_b0;
        data.config.constants.gamma_b1(sim_no,1) = collated_data(sim_no).input_data.sim_settings.gamma_b1;
    else
        data.config.constants.gamma_b0(sim_no,1) = collated_data(sim_no).input_data.sim_settings.fixed.wall.Gamma_b0;
        data.config.constants.gamma_b1(sim_no,1) = collated_data(sim_no).input_data.sim_settings.fixed.wall.Gamma_b1;
    end
    data.config.constants.mu(sim_no,1) = collated_data(sim_no).input_data.sim_settings.mu;
    data.config.constants.alpha(sim_no,1) = collated_data(sim_no).input_data.sim_settings.alpha;
    data.config.constants.p_drop(sim_no,1) = collated_data(sim_no).input_data.sim_settings.p_drop;
    
    %% Extract system chars
    if exist('system_chars', 'var')
        rel_vars = system_chars.Properties.VariableNames;
        for var_no = 1 : length(rel_vars)
            eval(['data.config.system_chars.' rel_vars{var_no} ' = system_chars.' rel_vars{var_no} '(sim_no);'])
        end
    end
end
clear sim_no

% rename variables
data.config.variations.param_names = strrep(data.config.variations.param_names, 'reg_vol', 'rfv');
data.config.variations.param_names = strrep(data.config.variations.param_names, 't_pf', 'pft');
    
% Add in haemodynamic params
data.haemods = haemodynamic_params;
data.waves.fs = collated_data(1).output_data.fs;
data.path_waves.fs = collated_data(1).output_data.fs;
data.waves.units.P = 'mmHg';
data.waves.units.U = 'm/s';
data.waves.units.A = 'm3';
data.path_waves.units = data.waves.units;
%data.waves.units.Q = 'm3/s';
data.waves.units.PPG = 'au';
data.path_waves.units.dist = 'm';

% Add in peripheral boundary conditions
data.config.network.wk_c = peripheral_chars.c;
data.config.network.wk_r = peripheral_chars.r_sum;

% system chars units
if exist('system_chars', 'var')
    data.config.system_chars.units.pvr = 'Pa s /m3';
    data.config.system_chars.units.pvc = 'm3 /Pa';
    data.config.system_chars.units.pvc_iw = 'm3 /Pa';
    data.config.system_chars.units.ac = 'm3 /Pa';
    data.config.system_chars.units.c = 'm3 /Pa';
    data.config.system_chars.units.tau = 's';
end

% Add in pulse wave indices
data.pw_inds = pw_inds;

%% Assess physiological plausibility of each virtual subject

data = assess_phys_plausibility(data);

%% Save files

orig_data = data; clear data

% - create variable for "aorta_finger" path waves
data = rmfield(orig_data, 'waves');
if sum(strcmp(fieldnames(data.path_waves), 'aorta_finger'))
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_r_subclavian')), rmfield(data.path_waves, 'aorta_r_subclavian'); end
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_foot')), rmfield(data.path_waves, 'aorta_foot'); end
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_brain')), rmfield(data.path_waves, 'aorta_brain'); end
    % save
    a = whos('data');
    if a.bytes > 1.8e9, save(PATHS.exported_data_mat_pwdb_data_w_aorta_finger_path, 'data', '-v7.3');
    else save(PATHS.exported_data_mat_pwdb_data_w_aorta_finger_path, 'data'), end
    clear data
end

% - create variable for "aorta_foot_p" path waves
data = rmfield(orig_data, 'waves');
if sum(strcmp(fieldnames(data.path_waves), 'aorta_foot'))
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_r_subclavian')), rmfield(data.path_waves, 'aorta_r_subclavian'); end
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_finger')), rmfield(data.path_waves, 'aorta_finger'); end
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_brain')), rmfield(data.path_waves, 'aorta_brain'); end
    data.path_waves.aorta_foot = rmfield(data.path_waves.aorta_foot, {'U', 'A'});
    % save
    a = whos('data');
    if a.bytes > 1.8e9, save(PATHS.exported_data_mat_pwdb_data_w_aorta_foot_path_p, 'data', '-v7.3');
    else save(PATHS.exported_data_mat_pwdb_data_w_aorta_foot_path_p, 'data'), end
    clear data
end

% - create variable for "aorta_foot_u" path waves
data = rmfield(orig_data, 'waves');
if sum(strcmp(fieldnames(data.path_waves), 'aorta_foot'))
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_r_subclavian')), rmfield(data.path_waves, 'aorta_r_subclavian'); end
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_finger')), rmfield(data.path_waves, 'aorta_finger'); end
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_brain')), rmfield(data.path_waves, 'aorta_brain'); end
    data.path_waves.aorta_foot = rmfield(data.path_waves.aorta_foot, {'P', 'A'});
    % save
    a = whos('data');
    if a.bytes > 1.8e9, save(PATHS.exported_data_mat_pwdb_data_w_aorta_foot_path_u, 'data', '-v7.3');
    else save(PATHS.exported_data_mat_pwdb_data_w_aorta_foot_path_u, 'data'), end
    clear data
end

% - create variable for "aorta_foot_a" path waves
data = rmfield(orig_data, 'waves');
if sum(strcmp(fieldnames(data.path_waves), 'aorta_foot'))
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_r_subclavian')), rmfield(data.path_waves, 'aorta_r_subclavian'); end
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_finger')), rmfield(data.path_waves, 'aorta_finger'); end
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_brain')), rmfield(data.path_waves, 'aorta_brain'); end
    data.path_waves.aorta_foot = rmfield(data.path_waves.aorta_foot, {'U', 'P'});
    % save
    a = whos('data');
    if a.bytes > 1.8e9, save(PATHS.exported_data_mat_pwdb_data_w_aorta_foot_path_a, 'data', '-v7.3');
    else save(PATHS.exported_data_mat_pwdb_data_w_aorta_foot_path_a, 'data'), end
    clear data
end

% - create variable for "aorta_brain" path waves
data = rmfield(orig_data, 'waves');
if sum(strcmp(fieldnames(data.path_waves), 'aorta_brain'))
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_r_subclavian')), rmfield(data.path_waves, 'aorta_r_subclavian'); end
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_foot')), rmfield(data.path_waves, 'aorta_foot'); end
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_finger')), rmfield(data.path_waves, 'aorta_finger'); end
    % save
    a = whos('data');
    if a.bytes > 1.8e9, save(PATHS.exported_data_mat_pwdb_data_w_aorta_brain_path, 'data', '-v7.3');
    else save(PATHS.exported_data_mat_pwdb_data_w_aorta_brain_path, 'data'), end
    clear data
end

% - create variable for "aorta_rsubclavian" path waves
data = rmfield(orig_data, 'waves');
if sum(strcmp(fieldnames(data.path_waves), 'aorta_r_subclavian'))
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_finger')), rmfield(data.path_waves, 'aorta_finger'); end
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_foot')), rmfield(data.path_waves, 'aorta_foot'); end
    if sum(strcmp(fieldnames(data.path_waves), 'aorta_brain')), rmfield(data.path_waves, 'aorta_brain'); end
    % save
    a = whos('data');
    if a.bytes > 1.8e9, save(PATHS.exported_data_mat_pwdb_data_w_aorta_rsubclavian_path, 'data', '-v7.3');
    else save(PATHS.exported_data_mat_pwdb_data_w_aorta_rsubclavian_path, 'data'), end
    clear data
end

% data without the arterial path waves
data = rmfield(orig_data, 'path_waves');
save(PATHS.exported_data_mat_pwdb_data, 'data')

end

function PrintFigs(h, paper_size, savepath)
set(h,'PaperUnits','inches');
set(h,'PaperSize', [paper_size(1), paper_size(2)]);
set(h,'PaperPosition',[0 0 paper_size(1) paper_size(2)]);
set(gcf,'color','w');
print(h,'-dpdf',savepath)
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
close all;

% save 
fid = fopen([savepath, '.txt'], 'w');
p = mfilename('fullpath');
p = strrep(p, '\', '\\');
fprintf(fid, ['Figures generated by:\n\n ' p '.m \n\n on ' datestr(today)]);
fclose all;

end

