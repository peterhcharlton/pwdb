function check_pwdb_processing(pwdb_no, collated_data)

fprintf('\n --- Checking PWDB processing ---')

% Setup paths with current simulation paths
PATHS = setup_paths_for_post_processing(pwdb_no);

% setup
up = setup_up;

% load data
if nargin == 2, up.use_exported_data = 0; end
if up.use_exported_data
    load(PATHS.exported_data_mat_pwdb_data)
else
    if ~exist('collated_data', 'var'), load(PATHS.collated_data), end
    data = convert_collated_data(collated_data);
end

% Check PWA
check_pwa(data,up);

% check PW extraction
check_pw_extraction(data, up)

end

function up = setup_up

up.pw_duration_range = [0.5, 1.0];

up.use_exported_data = 1;

end

function check_pwa(data,up)

% check whether haemod params have been extracted
no_subjs = length(data.haemods);
no_nans = nan(no_subjs,1);
for subj_no = 1 : length(data.haemods)
    a = struct2cell(data.haemods); a = cell2mat(a(:,:,1));
    no_nans(subj_no) = sum(isnan(a));
end


end

function check_pw_extraction(data, up)

% Identify names of waves to be checked
wave_names = fieldnames(data.waves);
wave_keep = true(length(wave_names),1);
for s = 1 : length(wave_names)
    curr_wave = wave_names{s};
    if ~sum(strcmp(curr_wave(1), {'P', 'U', 'A'}))
        wave_keep(s) = false;
    end
end
wave_names = wave_names(wave_keep);

% check each wave in turn
error_subjs_aortic = false(length(data.config.age),1);
error_subjs_others = false(length(data.config.age),1);
for wave_no = 1 : length(wave_names)
    lens = nan(length(data.config.age),1);
    for subj_no = 1 : length(data.config.age)
        eval(['curr_wave_data = data.waves.' wave_names{wave_no} '{subj_no};'])
        lens(subj_no) = length(curr_wave_data);        
    end
    tol = up.pw_duration_range*data.waves.fs;
    rel_els = lens > tol(2) | lens < tol(1);
    if ~isempty(strfind(wave_names{wave_no}, 'AorticRoot'))
        error_subjs_aortic(rel_els) = true;
    else
        error_subjs_others(rel_els) = true;
    end
end

% output answer
if sum(error_subjs_aortic)
    fprintf('\n Problems extracting aortic root waves')
end
if sum(error_subjs_others)
    fprintf('\n Problems extracting non-aortic root waves')
end
if ~sum(error_subjs_aortic) & ~sum(error_subjs_others)
    fprintf('\n No problems extracting waves')
end

end

function data = convert_collated_data(collated_data)

fprintf('\nConverting collated data: ')

% - sampling frequency
data.waves.fs = collated_data(1).output_data(1).fs;

% - age
data.config.age = nan(length(collated_data),1);
for subj_no = 1 : length(collated_data)
    data.config.age(subj_no) = collated_data(subj_no).input_data.sim_settings.age;
end

% - waves

% setup
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

for subj_no = 1 : length(collated_data)
    fprintf([num2str(subj_no) ', '])
    for site_no = 1 : length(sites)
        for sig_no = 1 : length(signals)
            
            % identify name of wave at this site
            wave_site_name = [signals{sig_no}, '_', sites{site_no}];
            
            % make variable if required
            if ~sum(strcmp(fieldnames(data.waves), wave_site_name))
                eval(['data.waves.' wave_site_name ' = cell(length(collated_data),1);'])
            end
            
            % - extract wave
            
            % Identify site-specific info
            curr_domain_no = site_domain_no(site_no);
            curr_domain_el = domains == curr_domain_no;
            curr_site_dist_prop = site_dist_prop(site_no);
            
            % Identify the relevant data for this simulation at this site
            seg_length = collated_data(subj_no).input_data.sim_settings.network_spec.length(curr_domain_no);
            dists = extractfield(collated_data(subj_no).output_data(curr_domain_el), 'distances');
            [~,rel_dist_el] = min(abs(curr_site_dist_prop*seg_length-dists)); clear dists
            
            % Extract wave
            if strcmp(signals{sig_no}, 'P')
                factor = 133.33;
            else
                factor = 1;
            end
            eval(['temp = collated_data(subj_no).output_data(curr_domain_el).' signals{sig_no} '(:,rel_dist_el)/factor;'])
            
            % store this wave
            eval(['data.waves.' wave_site_name '{subj_no} = temp;']);
            
            clear temp factor curr_sig curr_sig_data rel_wave_data seg_length rel_dist_el curr_domain_no curr_domain_el curr_site_dist_prop curr_site
        end
    end
end

end