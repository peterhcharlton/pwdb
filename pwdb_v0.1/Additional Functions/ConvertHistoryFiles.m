function savepath = ConvertHistoryFiles(up)
% Convert .his files created by Nektar into Matlab format.
%    Inputs:    up.dir      - the directory containing .his (and optionally
%                               .lum) files. If this is not specified then
%                               the user is prompted to select the
%                               directory manually. 
%               up.filename - a cell containing the filename(s) of .his
%                               files to be imported. If this is not
%                               specified then all .his files within the
%                               chosen directory are imported.
%               up.all_beats - a logical indicating whether data should be
%                               extracted from all beats (= true), or just
%                               the final beat of the simulation (= false).
%               all inputs are optional
%
%    Outputs:   history_files_data.m - a single Matlab file containing the
%                               data from all of the imported .his files,
%                               saved in the chosen directory.
%
%   Accompanying Article:
%       This code is provided to facilitate reproduction of the Pulse Wave
%       Database described in: 
%           Charlton P.H. et al. Modelling arterial pulse waves in healthy
%           ageing: in silico evaluation of haemodynamics and
%           cardiovascular indices, ~~ under review ~~  
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
%  Copyright (C) 2018  King's College London
% 
% Contributed to by: Marie Willemet, Jordi Alastruey, and Peter H. Charlton
% v.1.0


% setup universal parameters
if nargin == 0, up = struct; end
up = setup_up(up);

% order files according to simulation number
up.filename = order_according_to_sim(up.filename);

% convert files into matlab data
[data_from_files, filename_root] = convert_files_into_matlab_data(up);

% Eliminate non-relevant files
[data_from_files, filename_root, up] = eliminate_non_relevant_files(up, data_from_files, filename_root);

% separate data according to filename root
data = separate_data_according_to_filename_root(data_from_files, filename_root, up);

% Eliminate unwanted data
data = eliminate_unwanted_data(data, up);

% Ensure that single pulse waves are continuous
data = make_pulse_waves_continuous(data, up);

% save data in matlab file
savepath = [up.save_dir, 'history_files_data'];
save(savepath, 'data', '-v7.3')   % use v.7.3 for large files (> 2 GB)

end

function up = setup_up(up)

% Update dir field if necessary
if sum(strcmp(fieldnames(up), 'dir')) && ~strcmp(up.dir(end), filesep)
    up.dir = [up.dir, filesep];
end

% Identify missing fields
required_fields = {'dir', 'filename', 'all_beats', 'all_data', 'required_domains', 'required_distance_els', 'required_signals', 'save_dir', 'ds_factor', 'continuous_waves', 'find_pw_els'};
current_fields = fieldnames(up);
missing_fields = setxor(required_fields, current_fields);
% make sure required_domains comes before required_distance_els:
rel_el = find(strcmp(missing_fields, 'required_domains'));
missing_fields = [missing_fields(rel_el); missing_fields([1:rel_el-1, rel_el+1:end])];

% Fill in missing fields
for field_no = 1 : length(missing_fields)
    curr_missing_field = missing_fields{field_no};
    switch curr_missing_field
        case 'dir'
            up.dir = uigetdir('', 'Please select the directory containing files for analysis');
            if ~up.dir
                error('No directory selected')
            end
            if ~strcmp(up.dir(end), filesep)
                up.dir = [up.dir, filesep];
            end
        case 'filename'
            temp = dir([up.dir, '*.his']); 
            if ~isempty(temp)
                up.filename = extractfield(temp, 'name');
            else
                error('There aren''t any .his files in this directory')
            end
        case 'all_beats'
            up.all_beats = true;
        case 'all_data'
            up.all_data = true;
        case 'required_domains'
            %up.required_domains =      [1, 15, 21, 22, 42, 46, 49, 84, 87, 112];
            % this extracts points approx every 2cm along aortic-root - digital path
            %up.required_domains =      [1, 2, 14, 15, 19, 21, 22, 42, 46, 49, 84, 87, 108, 112];
            % this is a further extended set
            up.required_domains =      [1,             2, 7, 14, 16, 18, 19, 27, 28, 35, 37,  108,          39,          41,        42, 15,             84, 87,         21,       22,       112,             44,        46,           49, 72,79,65,96,71,3,4];
        case 'required_distance_els'
            % for the original set
            %up.required_distance_els = {1, [], [], 2, [],  2,  3,  2,  1,  3,  3,  2,  [],   3};
            % for the extended set
            up.required_distance_els = {1,             2, 3,  2,  3,  2,  2,  2,  2,  2,  2,    3,           1,           3,         2,  2,              3,  3,          2,        3,         3,              2,         2,            3,  3, 2, 2, 2, 2,2,2};
            % this extracts points approx every 2cm along aortic-root - digital path
            rel_path_domains = [1,2,14,19,21,22,108,112];
            % this extracts points approx every 2cm along aortic-root - digital and aortic-root - foot and aortic-root - brain paths
            rel_path_domains = [1,2,14,19,21,22,108,112, 18,27,28,35,37,39,41,42,44,46,49,15,16,79,65,96,71];
            % Extracts right subclavian
            rel_path_domains = [rel_path_domains,3,4];
            for dom_no = 1 : length(rel_path_domains)
                rel_el = find(up.required_domains == rel_path_domains(dom_no));
                up.required_distance_els{rel_el} = 'all';
            end
            clear dom_no rel_el rel_path_domains
            
        case 'required_signals'
            up.required_signals = {'P', 'U', 'A', 'Q1D', 'Q_out'};  %{'P', 'Pe', 'U', 'Q', 'A', 'Q1D', 'Q_out'};  %{'P', 'Pe', 'U', 'Q', 'A', 'Q1D', 'Q_out'}; % {'P', 'Pe', 'U', 'Q', 'A'};
        case 'save_dir'
            up.save_dir = up.dir;
        case 'ds_factor'
            up.ds_factor = 1;
        case 'continuous_waves'
            up.continuous_waves = 1;
        case 'find_pw_els'
            up.find_pw_els = 0;
    end
end

end

function ordered_filenames = order_according_to_sim(filenames)

counter = 0;
for s = 1 : length(filenames)
    if strcmp(filenames{s}(1:4), 'sim_')
        counter = counter+1;
    end
end
dont_use = true;
sim_nos = nan(length(filenames),1);
if counter == length(filenames)
    dont_use = false;
    for s = 1 : length(filenames)
        temp = strfind(filenames{s}, '_');
        if length(temp) ~=2
            dont_use = true;
        else
            sim_nos(s) = str2double(filenames{s}(temp(1)+1:temp(2)-1));
        end
        clear temp
    end
end
clear counter s

if ~dont_use
    [~, order] = sort(sim_nos);
    ordered_filenames = filenames(order);
end
clear dont_use order sim_nos

end

function [data_from_files, filename_root] = convert_files_into_matlab_data(up)

counter_no = 0;
for file_no = 1 : length(up.filename)
    curr_file = up.filename{file_no};
    
    % see if this file should be skipped
    if ~up.all_data
        temp1 = strfind(curr_file, '_'); temp1 = temp1(end);
        temp2 = strfind(curr_file, '.'); temp2 = temp2(end);
        seg_no = str2double(curr_file(temp1+1:temp2));
        if ~sum(seg_no == up.required_domains)
            temp = strfind(curr_file, '_');
            %filename_root{file_no,1} = curr_file(1:temp(end)-1); clear temp curr_file temp1 temp2 seg_no
            continue
        end
        clear temp temp1 temp2 seg_no
    end
    
    % Import data from this history file
    fprintf(['\n - ' curr_file])
    curr_file_path =[up.dir, curr_file];
    temp = ReadHistoryFiles(curr_file_path, up);
    temp = downsample_data(temp, up);
    data_from_files{file_no} = temp; clear temp field_no temp_fields
    
    % see if there is any lumped parameter data for this file
    temp1 = strfind(curr_file, '_'); temp1 = temp1(end);
    temp2 = strfind(curr_file, '.'); temp2 = temp2(end);
    trial_lum_name = [up.dir, curr_file(1:temp1-1), '_out_' curr_file(temp1+1:temp2), 'lum'];
    if exist(trial_lum_name, 'file')
        lum_data_from_files = ReadLumpedFiles(trial_lum_name, up);  clear trial_lum_name
        lum_data_from_files = downsample_data(lum_data_from_files, up);
        fields = fieldnames(lum_data_from_files); fields = fields(~strcmp(fields,'fs'));
        for field_no = 1 : length(fields)
            eval(['data_from_files{file_no}.' fields{field_no} ' = lum_data_from_files.' fields{field_no} ';']);
        end
        clear fields lum_data_from_files
    end
    
    % Extract data for a single beat if required
    if ~up.all_beats
        data_from_files{file_no} = Extract_data_for_single_beat(data_from_files{file_no}, curr_file, up);
    end
    
    % store filename root
    temp = strfind(curr_file, '_');
    filename_root{file_no,1} = curr_file(1:temp(end)-1); clear temp curr_file_path
    clear curr_file
end
clear file_no

end

function [data_from_files, filename_root, up] = eliminate_non_relevant_files(up, data_from_files, filename_root)

if ~up.all_data
    rel_files = ~cellfun(@isempty, filename_root);
    filename_root = filename_root(rel_files);
    data_from_files = data_from_files(rel_files);
    up.filename = up.filename(rel_files);
    clear rel_files
end

end

function data = eliminate_unwanted_data(data, up)
        
sims = fieldnames(data);
if ~up.all_data
    
    for sim_no = 1 : length(sims)
        curr_sim = sims{sim_no};
        eval(['sim_data = data.' curr_sim ';'])
        
        % eliminate unwanted signals
        curr_signals = fieldnames(sim_data);
        curr_signals = setxor(curr_signals, {'domain_no', 'distances', 'fs', 'units', 'start_sample'});
        unwanted_signals = setxor(curr_signals, up.required_signals); clear curr_signals
        for sig_no = 1 : length(unwanted_signals)
            if sum(strcmp(fieldnames(sim_data), unwanted_signals{sig_no}))
                sim_data = rmfield(sim_data, unwanted_signals{sig_no});
            end
        end
        clear unwanted_signals
        actual_signals = intersect(fieldnames(sim_data), up.required_signals);
        
        % eliminated unwanted domains
        curr_domains = extractfield(sim_data, 'domain_no');
        [~,req_rows,~] = intersect(curr_domains, up.required_domains); clear curr_domains
        sim_data = sim_data(req_rows); clear req_rows
        
        % eliminate unwanted measurement points
        for domain_el = 1 : length(sim_data)
            curr_domain = sim_data(domain_el).domain_no;
            req_measurement_points = up.required_distance_els{up.required_domains == curr_domain}; clear curr_domain
            % skip if all measurement points are to be extracted from this domain
            if ischar(req_measurement_points) && strcmp(req_measurement_points, 'all')
                continue
            end
            % Remove unwanted signal values
            for signal_no = 1 : length(actual_signals)
                curr_signal = actual_signals{signal_no};
                if ~sum(strcmp(curr_signal, {'Q1D', 'Q_out', 'P1D', 'PC'}))
                    eval(['sim_data(domain_el).' curr_signal ' = sim_data(domain_el).' curr_signal '(:,req_measurement_points);'])
                else
                    eval(['sim_data(domain_el).' curr_signal ' = sim_data(domain_el).' curr_signal ';'])
                end
                clear curr_signal
            end
            clear signal_no
            % Remove unwanted distances
            sim_data(domain_el).distances = sim_data(domain_el).distances(req_measurement_points);
            % Remove unwanted start samples
            sim_data(domain_el).start_sample = sim_data(domain_el).start_sample(req_measurement_points);
            clear req_measurement_points
        end
        clear domain_el
        
        % store reduced data
        eval(['data.' curr_sim ' = sim_data;'])
        
    end
    
end

end
    
function temp = downsample_data(temp, up)

% downsample data
temp_fields = fieldnames(temp);
temp_fields = temp_fields(~strcmp(temp_fields, 'fs') & ~strcmp(temp_fields, 'distances'));
for field_no = 1 : length(temp_fields)
    eval(['curr_field_data = temp.' temp_fields{field_no} ';']);
    if up.ds_factor ~=1
        new_field_data = downsample(curr_field_data, up.ds_factor);
        eval(['temp.' temp_fields{field_no} ' = new_field_data;']);
    end
    clear new_field_data curr_field_data
end

% re-calculate sampling frequency
temp.fs = temp.fs/up.ds_factor;

end

function data = separate_data_according_to_filename_root(data_from_files, filename_root, up)

% Identify filename roots
filename_root = strrep(filename_root,'-','_');
[filename_roots, els, ~] = unique(filename_root);
[~, order] = sort(els);
filename_roots = filename_roots(order);

char_filename_roots = filename_roots;
for s = 1 : length(filename_roots)
    if ~isstrprop(filename_roots{s}(1),'alpha')
        char_filename_roots{s} = ['sim_' filename_roots{s}];
    end
end
clear s

for filename_root_no = 1 : length(filename_roots)
    rel_files = find(strcmp(filename_root, filename_roots{filename_root_no}));
    
    for rel_file_no = 1 : length(rel_files)
        
        % store domain no
        domain_no = up.filename{rel_files(rel_file_no)};
        temp = strfind(domain_no, '_');
        temp2 = strfind(domain_no, '.his');
        domain_no = str2double(domain_no(temp(end)+1:temp2-1));
        eval(['data.' char_filename_roots{filename_root_no} '(rel_file_no).domain_no = domain_no;'])
        clear domain_no temp temp2
        
        % store each variable's data in turn (including units)
        temp = data_from_files{rel_files(rel_file_no)};
        vars = fieldnames(temp);
        units = struct;
        for var_no = 1 : length(vars)
            eval(['data.' char_filename_roots{filename_root_no} '(rel_file_no).' vars{var_no} ' = temp.' vars{var_no} ';'])
            eval(['units.' vars{var_no} ' = find_units(''' vars{var_no} ''');']);
        end
        eval(['data.' char_filename_roots{filename_root_no} '(rel_file_no).units = units;']);
        clear units vars var_no temp
        
    end
    clear rel_file_no rel_files
    
    % sort data according to domain no
    eval(['rel_data = data.' char_filename_roots{filename_root_no} ';'])
    domain_nos = extractfield(rel_data, 'domain_no');
    [~, rel_order] = sort(domain_nos);
    rel_data = rel_data(rel_order);
    eval(['data.' char_filename_roots{filename_root_no} ' = rel_data;'])
    clear rel_data rel_order domain_nos
    
end
clear filename_root_no filename_roots filename_root data_from_files char_filename_roots

end

function history_file_data = ReadHistoryFiles(curr_file_path, up)

% Identify header lines and measurement points:
fid = fopen(curr_file_path);
header_line_log = true; line_no = 0; history_file_data.distances = [];
while header_line_log
    curr_line = fgetl(fid);
    line_no = line_no + 1;
    if ~strcmp(curr_line(1), '#')
        header_line_log = false;
    end
    if strcmp(curr_line(3:8), 'Point ')
        history_file_data.distances(end+1) = str2double(curr_line(16:23));
    end
    if strcmp(curr_line(1:3), '# t')
        header_text = strrep(curr_line, '(x,t)', '');
        header_text = strrep(header_text, ' ', '');
        header_text = strrep(header_text, '#', '');
        headers = textscan(header_text, '%s', 'Delimiter', ',');
        headers = headers{1}; clear header_text
    end
end
fclose all; 
no_header_lines = line_no - 1;
clear curr_line header_line_log fid line_no


% Import Nektar data:
raw = importdata(curr_file_path, ' ', no_header_lines);

raw = raw.data;
for col_no = 1 : length(headers)
    eval(['ext_data.' headers{col_no} ' = raw(:,col_no);']);    
end
history_file_data.fs = 1/median(diff(unique(ext_data.t)));
history_file_data.fs = round(1e6*history_file_data.fs)/1e6;
t_col = strcmp(headers, 't');
point_col = strcmp(headers, 'point');
points = raw(:, point_col);
raw = raw(:, ~t_col & ~point_col);
headers = headers(~t_col & ~point_col);

% Separate data according to each measurement location:
for header_no = 1 : length(headers)
    curr_header = headers{header_no};
    temp = raw(:,header_no);
    header_data = nan(ceil(length(temp)/length(history_file_data.distances)), length(history_file_data.distances));
    for point_no = 1 : length(history_file_data.distances)
        point_temp = temp(points == point_no);
        header_data(1:length(point_temp),point_no) = point_temp;
        clear point_temp
    end
    clear temp
    eval(['history_file_data.' curr_header ' = header_data;']); clear header_data
end

end

function lumped_file_data = ReadLumpedFiles(curr_file_path, up)

% Identify header lines and measurement points:
fid = fopen(curr_file_path);
header_line_log = true; line_no = 0;
while header_line_log
    curr_line = fgetl(fid);
    line_no = line_no + 1;
    if ~strcmp(curr_line(1), '#')
        header_line_log = false;
    end
    if strcmp(curr_line(1:3), '# t')
        header_text = curr_line;
        header_text = strrep(header_text, '# ', '');
        old_len = 0;
        while length(header_text) ~= old_len
            old_len = length(header_text);
            header_text = strrep(header_text, '  ', ' ');
        end
        header_text = strrep(header_text, ' ', ',');
        headers = textscan(header_text, '%s', 'Delimiter', ',');
        headers = headers{1}; clear header_text
    end
end
fclose all; 
no_header_lines = line_no - 1;
clear curr_line header_line_log fid line_no


% Import Nektar data:
raw = importdata(curr_file_path, ' ', no_header_lines);
raw = raw.data;
for col_no = 1 : length(headers)
    eval(['ext_data.' headers{col_no} ' = raw(:,col_no);']);    
end
lumped_file_data.fs = 1/median(diff(unique(ext_data.t)));
lumped_file_data.fs = round(1e6*lumped_file_data.fs)/1e6;
t_col = strcmp(headers, 't');
raw = raw(:, ~t_col);
headers = headers(~t_col);

% Extract each variable
for header_no = 1 : length(headers)
    curr_header = headers{header_no};
    eval(['lumped_file_data.' curr_header ' = raw(:,header_no);']); clear header_data
end

end

function units = find_units(var_name)

switch var_name
    case 'P'
        units = 'Pa';
    case 'PC'
        units = 'Pa';
    case 'P1D'
        units = 'Pa';
    case 'Pe'
        units = 'Pa';
    case 'Pext'
        units = 'Pa';
    case 'U'
        units = 'm/s';
    case 'Q'
        units = 'm3/s';
    case 'Q_out'
        units = 'm3/s';
    case 'Q1D'
        units = 'm3/s';
    case 'A'
        units = 'm2';
    case 'distances'
        units = 'm';
    case 'fs'
        units = 'Hz';
    case 'start_sample'
        units = 'no samples';
end

end

function new_data_from_file = Extract_data_for_single_beat(data_from_file, curr_file, up)

pulse_wave_duration_in_samples = nan;

% cycle through each measurement point
for measurement_pt_no = 1 : length(data_from_file.distances)
    clear beat_onsets minima a keep_minima rel_els heights heights_thresh keep_minima minima_no
    
    % identify start and end of final complete beat
    
    % Find minima in pressure (P)
    a = data_from_file.P(:,measurement_pt_no);
    fs = data_from_file.fs;
    
    
    % Find relevant indices
    use_original_method = 0;
    if use_original_method
        [relevant_inds, pulse_wave_duration_in_samples] = find_rel_inds_using_orig_method(a, fs, measurement_pt_no, pulse_wave_duration_in_samples);
    else
        [relevant_inds, pulse_wave_duration_in_samples] = find_rel_inds_using_new_method(a,fs, measurement_pt_no, pulse_wave_duration_in_samples);
    end
    
    % extract and store relevant part of signal (i.e. last compelete beat)
    vars = fieldnames(data_from_file);
    for var_no = 1 : length(vars)
        curr_var = vars{var_no};
        if strcmp(curr_var, 'fs') || strcmp(curr_var, 'distances')
            eval(['new_data_from_file.' curr_var ' = data_from_file.' curr_var ';']);
            continue
        end
        if ~sum(strcmp(curr_var, {'P1D', 'Q1D', 'Q_out', 'PC'}))
            eval(['new_data_from_file.' curr_var '(:,measurement_pt_no) = data_from_file.' curr_var '(relevant_inds,measurement_pt_no);']);
        else
            eval(['new_data_from_file.' curr_var ' = data_from_file.' curr_var '(relevant_inds,1);']);
        end
    end
    
    
    % Note the time of this beat onset
    new_data_from_file.start_sample(measurement_pt_no) = relevant_inds(1);
        
end

if ~exist('new_data_from_file', 'var')
    a = 1;
end

end

function [relevant_inds, pulse_wave_duration_in_samples] = find_rel_inds_using_orig_method(a, fs, measurement_pt_no, pulse_wave_duration_in_samples)
minima = 1 + find( ...
    a(3:end-1)>a(2:end-2) & a(1:end-3)>a(2:end-2) | ...
    (a(4:end)>a(3:end-1) & a(1:end-3)>a(2:end-2) & a(3:end-1) == a(2:end-2)) ...
    );

% Eliminate repeated minima (or those which are very close to each other)
% - repeated
keep_minima = true(size(minima));
for minima_no = 2:length(minima)
    if minima(minima_no) == minima(minima_no-1)+1
        keep_minima(minima_no) = false;
    end
end
minima = minima(keep_minima);
% - close to each other
keep_minima = true(size(minima));
for minima_no = 1:length(minima)-1
    diffs = minima - minima(minima_no);
    tol_samps = fs * 0.03;  % within 0.030 s
    if sum(diffs > 0 & diffs < tol_samps)
        keep_minima(minima_no) = false;
    end
end
minima = minima(keep_minima);

% determine how high each minimum is in relation to the previous 3 s of
% data
time_int = 1.5; % in secs
heights = nan(length(minima),1);
for minimum_no = 1 : length(minima)
    prev_time_int.deb = minima(minimum_no) - round(time_int*fs);
    prev_time_int.fin = minima(minimum_no);
    if prev_time_int.deb <= 0
        heights(minimum_no) = nan;
    else
        heights(minimum_no) = (a(minima(minimum_no))-min(a(prev_time_int.deb:prev_time_int.fin)))/range(a(prev_time_int.deb:prev_time_int.fin));
    end
end

heights_inc = 0.06;
heights_thresh = 0.175-heights_inc;
successful_extraction = 0;
while ~successful_extraction && heights_thresh < 0.5
    heights_thresh = heights_thresh+heights_inc;
    
    % identify reliable beat onsets according to amplitude
    rel_els = heights<heights_thresh;
    if sum(rel_els) < 3
        continue
    end
    
    try
        beat_onsets = minima(rel_els);
        new_heights = heights(rel_els); clear rel_els
        % identify reliable beat onsets according to timing
        thresh_secs = 0.2;  % threshold in secs
        thresh_samps = thresh_secs*fs;
        beat_onsets_to_eliminate = [];
        for beat_onset_no = 1 : length(beat_onsets)-1
            curr_beat_onsets = beat_onsets(beat_onset_no:beat_onset_no+1);
            if range(curr_beat_onsets) < thresh_samps
                % then identify whether to get rid of the first or second of the two:
                no_samps_to_consider = fs*2;
                max_variation = range(a(curr_beat_onsets(1)-no_samps_to_consider:curr_beat_onsets(1)));
                if range(a(curr_beat_onsets)) < 0.1*max_variation
                    % eliminate the first one if they are about the same height
                    temp = 1;
                else
                    [~, temp] = max(a(curr_beat_onsets));
                end
                clear max_variation no_samps_to_consider
                rel_beat_onset_el = temp+beat_onset_no-1;
                beat_onsets_to_eliminate = [beat_onsets_to_eliminate, rel_beat_onset_el];
            end
        end
        rel_els = setxor(1:length(beat_onsets), beat_onsets_to_eliminate);
        beat_onsets = beat_onsets(rel_els);
        new_heights = new_heights(rel_els); clear rel_els
        
        % determine how high the max value is between each minimum and the
        % next one in relation to the prior signal
        time_int = 2; % in secs
        max_heights = nan(length(beat_onsets),1);
        for beat_onset_no = 1 : length(beat_onsets)-1
            current_beat = a(beat_onsets(beat_onset_no):beat_onsets(beat_onset_no+1));
            prev_time_int.deb = beat_onsets(beat_onset_no) - round(time_int*fs);
            prev_time_int.deb = max([1, prev_time_int.deb]);
            prev_time_int.fin = beat_onsets(beat_onset_no);
            max_heights(beat_onset_no) = (max(current_beat)-min(current_beat))/(max(a(prev_time_int.deb:prev_time_int.fin)) - min(a(prev_time_int.deb:prev_time_int.fin)));
        end
        
        % identify reliable beat onsets according to max height between
        rel_els = max_heights > 0.5 | isnan(max_heights);
        beat_onsets = beat_onsets(rel_els);
        
        % check that the last two beat onsets are a similar time apart
        if range(beat_onsets(end-1:end)) < 0.95*range(beat_onsets(end-2:end-1))
            beat_onsets = beat_onsets(1:end-1);
        end
        
        % check that there is some upslope after the final beat onset
        if a(beat_onsets(end)) == max(a(beat_onsets(end):end))
            beat_onsets = beat_onsets(1:end-1);
        end
        
        % check that there is sufficient signal in the relevant data for this
        % last complete beat
        if measurement_pt_no ~= 1
            curr_relevant_inds = beat_onsets(end-1):beat_onsets(end-1)+pulse_wave_duration_in_samples-1;
            if curr_relevant_inds(end) > length(a)
                beat_onsets = beat_onsets(1:end-1);
            end
            clear curr_relevant_inds
        end
        
        % identify two last beat onsets:
        if ~isempty(beat_onsets)
            last_beat_onsets = beat_onsets(end-1:end);
        else
            last_beat_onsets = [1, length(a)];
            fprintf('\n ---- Couldn''t find a reliable complete beat, so outputting all data')
        end
        
        %  Identify relevant indices for the last complete beat
        if measurement_pt_no == 1
            relevant_inds = last_beat_onsets(1):last_beat_onsets(2)-1;
            pulse_wave_duration_in_samples = length(relevant_inds);
        else
            relevant_inds = last_beat_onsets(1):last_beat_onsets(1)+pulse_wave_duration_in_samples-1;
        end
        
        % See whether this one needs manually annotating
        curr_p_wav = a(relevant_inds);
        curr_p_wav = (curr_p_wav - min(curr_p_wav))/range(curr_p_wav);
        t = find(curr_p_wav > 0.1,1)/fs;
        % If the upslope isn't straightaway (i.e. after the first 60 ms)
        if t > 0.06
            % Look up values
            if up.find_pw_els
                relevant_inds = find_PW_els(curr_file, measurement_pt_no);
            end
            % or manually annotate if they're not there
            if isnan(relevant_inds(1))
                initial_el = (last_beat_onsets(1)-(data_from_file.fs*1.2));
                plot(a( initial_el : end) )
                title = 'Couldn''t find onset of this beat';
                dims = [1 35];
                opts.WindowStyle = 'normal';
                beep,pause(0.5), beep
                if measurement_pt_no == 1
                    prompt = {'Enter sample of second last beat onset:', 'Enter sample of last beat onset:'};
                    definput = {'0','0'};
                    answer = inputdlg(prompt,title,dims,definput,opts);
                    last_beat_onsets = [str2double(answer(1)), str2double(answer(2))];
                else
                    prompt = {'Enter sample of second last beat onset:'};
                    definput = {'0'};
                    answer = inputdlg(prompt,title,dims,definput,opts);
                    last_beat_onsets = str2double(answer(1));
                end
                close all
                last_beat_onsets = last_beat_onsets+initial_el-1;
                % re-find relevant indices
                if measurement_pt_no == 1
                    relevant_inds = last_beat_onsets(1):last_beat_onsets(2)-1;
                    pulse_wave_duration_in_samples = length(relevant_inds);
                else
                    relevant_inds = last_beat_onsets(1):last_beat_onsets(1)+pulse_wave_duration_in_samples-1;
                end
                fprintf(['\n - meas pt ' num2str(measurement_pt_no) ', start ' num2str(relevant_inds(1)) ', end ' num2str(relevant_inds(end))])
                clear answer definput prompt opts dims title initial_el
            end
        end
        
        successful_extraction = 1;
        
    catch
        if heights_thresh > 0.6
            error('Not working')
        end
        continue
    end
    
end

end

function data = make_pulse_waves_continuous(data, up)

% Check that this should be done
if ~(~up.all_beats && up.continuous_waves)
    return
end

% Cycle through each simulation
sim_names = fieldnames(data);
for sim_no = 1 : length(sim_names)
    curr_sim_name = sim_names{sim_no};
    eval(['curr_sim_data = data.' curr_sim_name ';'])
        
    % cycle through each domain
    for row_no = 1 : length(curr_sim_data)
        
        % Cycle through each signal
        signals = fieldnames(curr_sim_data);
        signals = signals(~strcmp(signals, 'domain_no') & ~strcmp(signals, 'distances') & ~strcmp(signals, 'fs') & ~strcmp(signals, 'start_sample') & ~strcmp(signals, 'units'));
        for sig_no = 1 : length(signals)
            
            % extract data for this signal
            curr_sig_name = signals{sig_no};
            eval(['curr_sig = data.' curr_sim_name '(row_no).' curr_sig_name ';'])
            
            % for each measurement point
            for meas_pt_no = 1 : size(curr_sig,2)
                curr_pw = curr_sig(:,meas_pt_no);
                cont_pw = calc_static_wave(curr_pw);
                eval(['data.' curr_sim_name '(row_no).' curr_sig_name '(:,meas_pt_no) = cont_pw;'])
            end
            
        end
        
    end
    
end

end

function [relevant_inds, pulse_wave_duration_in_samples] = find_rel_inds_using_new_method(sig, fs, measurement_pt_no, pulse_wave_duration_in_samples)

%% Identify beats

% - find derivative
mov_avg_duration = 0.01;
no_samps = mov_avg_duration*fs;
sig_mov_avg = movmean(sig,5);
deriv = [0;diff(sig_mov_avg)*fs];

% - find maxima in derivative
deriv_maxima = find_max(deriv);

% - identify relevant maxima (those above a threshold value)
thresh = movmax(deriv, 2*fs);
thresh = 0.95*[thresh(2*fs:end); ones(2*fs-1,1)*thresh(end)];
rel_max = deriv_maxima(deriv(deriv_maxima)>thresh(deriv_maxima));

% - refine candidates
finished_refining = 0;
while ~finished_refining
    thresh_secs = 0.5;
    thresh_samps = thresh_secs*fs;
    rel_max_to_eliminate = [];
    for rel_max_no = 1 : length(rel_max)-1
        curr_rel_max = rel_max(rel_max_no:rel_max_no+1);
        if range(curr_rel_max) < thresh_samps
            % then eliminate the one with the lower derivative
            [~, temp] = min(deriv(curr_rel_max));
            clear max_variation no_samps_to_consider
            rel_max_el = temp+rel_max_no-1;
            rel_max_to_eliminate = [rel_max_to_eliminate, rel_max_el];
        end
    end
    rel_els = setxor(1:length(rel_max), rel_max_to_eliminate);
    rel_max = rel_max(rel_els); clear rel_els
    if sum(diff(rel_max)<thresh_samps) == 0
        finished_refining = 1;
    end
end

%% Identify beat peaks
beat_peaks = nan(length(rel_max),1);
for beat_no = 1 : length(rel_max)
    
    if beat_no < length(rel_max)
        rel_range = rel_max(beat_no):rel_max(beat_no+1);
    else
        rel_range = rel_max(beat_no):length(sig);
    end
    [~, beat_peaks(beat_no)] = max(sig(rel_range));
    beat_peaks(beat_no) = beat_peaks(beat_no)+rel_max(beat_no)-1;
    
end

%% Identify beat onsets
minima = find_min(sig);
deriv_minima = find_min(deriv);
pos_zero_cross = find_pos_zero_cross(deriv);
beat_onsets = nan(length(rel_max),1);
for beat_no = 1 : length(rel_max)
    
    % identify minimum before this beat element
    curr_min = minima(find(minima<rel_max(beat_no),1, 'last'));
    if beat_no == 1
        curr_min = 1;
    end
    
    % identify zero first derivative before this beat element
    curr_zero = pos_zero_cross(find(pos_zero_cross<rel_max(beat_no),1, 'last'));
    
    % identify any minima in first derivative which could be beat onsets
    curr_max = beat_peaks(beat_no);
    deriv_minima_heights = (sig(deriv_minima)-sig(curr_min))/(sig(curr_max)-sig(curr_min));
    curr_deriv_minimum = deriv_minima(find(deriv_minima < rel_max(beat_no) & deriv_minima_heights < 0.2,1,'last'));
    
    % choose current beat onset
    curr_beat_onset = max([curr_min, curr_zero, curr_deriv_minimum]);
    
    % store this beat onset
    beat_onsets(beat_no) = curr_beat_onset;
end

%% Identify relevant inds
if ~isnan(pulse_wave_duration_in_samples)
    relevant_inds = beat_onsets(end-1) : beat_onsets(end-1) + pulse_wave_duration_in_samples-1;
else
    relevant_inds = beat_onsets(end-1) : beat_onsets(end)-1;
end

% check that these don't go past the end of the signal
if relevant_inds(end) > length(sig)
    if ~isnan(pulse_wave_duration_in_samples)
        relevant_inds = beat_onsets(end-2) : beat_onsets(end-2) + pulse_wave_duration_in_samples-1;
    else
        relevant_inds = beat_onsets(end-2) : beat_onsets(end-1)-1;
    end
end

% store duration of PW
pulse_wave_duration_in_samples = length(relevant_inds);

end

function maxima = find_max(sig)

maxima1 = find(sig(2:end-1) > sig(1:end-2) & sig(3:end) < sig(2:end-1)) + 1;
maxima2 = find(sig(2:end-2) > sig(1:end-3) & sig(2:end-2) == sig(3:end-1) & sig(4:end) < sig(3:end-1)) + 2;

maxima = unique([maxima1; maxima2]);

end

function minima = find_min(sig)

minima1 = find(sig(2:end-1) < sig(1:end-2) & sig(3:end) > sig(2:end-1)) + 1;
minima2 = find(sig(2:end-2) < sig(1:end-3) & sig(2:end-2) == sig(3:end-1) & sig(4:end) > sig(3:end-1)) + 2;

minima = unique([minima1; minima2]);

end

function pos_zero_cross = find_pos_zero_cross(deriv)

pos_zero_cross = find(deriv(2:end)>0 & deriv(1:end-1) <= 0)+1;

end

function static_wave = calc_static_wave(orig_wave)

orig_wave = orig_wave(:);

% Calculate expected position of next point
next_point = orig_wave(end) + diff(orig_wave(end-1:end));

% If the wave was static, then this next point would be equal to the first point
% So, we can make the wave static by making this next point equal to the first
temp = linspace(0,orig_wave(1)-next_point, length(orig_wave)); temp = temp(:);
static_wave = orig_wave + temp;

end