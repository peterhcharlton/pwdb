function pwdb_pwa(pwdb_no)
% PWDB_PWA derives parameters from pulse waves, and saves them in
% Matlab format.
%
%               pwdb_pwa
%
%   Inputs:     the 'collated_data.mat' file produced by 'extract_pwdb_simulation_data.m'.
%
%   Outputs:    - Matlab files containing the pulse wave parameters:
%                       pulse_wave_inds.mat
%                       pulse_wave_vels.mat
%                       haemodynamic_params.mat
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
% v.1.0  Contributed to by: Marie Willemet, Jordi Alastruey, and Peter H. Charlton

fprintf('\n --- Performing Pulse Wave Analysis ---')

%% Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTINGS TO CHANGE: This function specifies where to save the outputs   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PATHS = setup_paths_for_post_processing(pwdb_no);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

props.carotid = 0.5;

% Extract pulse wave velocities
if ~exist(PATHS.pulse_wave_vels, 'file')
    if ~exist('collated_data', 'var'), load(PATHS.collated_data), end
    extract_pulse_wave_velocities(PATHS, props, collated_data);
end

% Extract pulse wave parameters
if ~exist(PATHS.pulse_wave_inds, 'file')
    if ~exist('collated_data', 'var'), load(PATHS.collated_data), end
    extract_pulse_wave_inds(PATHS, props, collated_data);
end

% Extract haemodynamic parameters
if ~exist(PATHS.haemodynamic_params, 'file')
    if ~exist('collated_data', 'var'), load(PATHS.collated_data), end
    extract_haemodynamic_parameters(PATHS, props, collated_data);
end

fprintf('\n --- Finished analysing PWDB pulse waves ---\n')

end

function extract_pulse_wave_velocities(PATHS, props, collated_data)

fprintf('\n Extracting pulse wave velocities:')

% setup
sim_data.domain_nos = extractfield(collated_data(1).output_data, 'domain_no');

%% Extract PWVs
fprintf('\n - PWVs: ')
pwv_types = {'carotid_femoral', 'carotid_radial', 'carotid_ankle', 'carotid_brachial', ...
    'brachial_femoral', 'brachial_radial', 'brachial_ankle', ...
    'radial_femoral', 'radial_ankle', ...
    'femoral_ankle', ...
    'aorta_iliacbif'};

TT_waveform_ind = 1; % for pressure waves
for sim_no = 1 : length(collated_data)
    fprintf([num2str(sim_no) ', ']);
    sim_data.output_data = collated_data(sim_no).output_data;
    sim_data.input_data = collated_data(sim_no).input_data;
    
    for pwv_type_no = 1 : length(pwv_types)
        curr_pwv_type = pwv_types{pwv_type_no};
        
        % identify two sites
        temp = strfind(curr_pwv_type, '_');
        site1.name = curr_pwv_type(1:temp-1);
        site2.name = curr_pwv_type(temp+1:end);
        
        % Extract information for each site
        for site_no = 1 : 2
            % identify this site
            eval(['curr_site = site' num2str(site_no) ';'])
            % extract geometric data for this site
            switch curr_site.name
                case 'carotid'
                    curr_site.preceeding_domains = [1,2];
                    curr_site.domain_no = 15;
                    curr_site.distance_prop = props.carotid;
                case 'femoral'
                    curr_site.preceeding_domains = [1,2,14,18,27,28,35,37,39,41,42,44];
                    curr_site.domain_no = 46;
                    curr_site.distance_prop = 0.5;
                case 'radial'
                    curr_site.preceeding_domains = [1,2,14,19,21];
                    curr_site.domain_no = 22;
                    curr_site.distance_prop = 1;
                case 'ankle'
                    curr_site.preceeding_domains = [1,2,14,18,27,28,35,37,39,41,42,44,46];
                    curr_site.domain_no = 49;
                    curr_site.distance_prop = 1;
                case 'brachial'
                    curr_site.preceeding_domains = [1,2,14,19];
                    curr_site.domain_no = 21;
                    curr_site.distance_prop = 0.75;
                case 'aorta'
                    curr_site.preceeding_domains = [];
                    curr_site.domain_no = 1;
                    curr_site.distance_prop = 0;
                case 'iliacbif'
                    curr_site.preceeding_domains = [1,2,14,18,27,28,35,37,39];
                    curr_site.domain_no = 41;
                    curr_site.distance_prop = 1;
            end
            % calculate path length
            curr_site.len_aorta_site = sum([sim_data.input_data.sim_settings.network_spec.length(curr_site.preceeding_domains); sim_data.input_data.sim_settings.network_spec.length(curr_site.domain_no)*curr_site.distance_prop]);
            
            % extract waveform data for this site
            curr_site.domain_row = find(sim_data.domain_nos == curr_site.domain_no);
            curr_site.seg_len = sim_data.input_data.sim_settings.network_spec.length(curr_site.domain_no);
            [~, curr_site.distance_el] = min(abs(sim_data.output_data(curr_site.domain_row).distances-(curr_site.distance_prop*curr_site.seg_len)));
            curr_site.v = sim_data.output_data(curr_site.domain_row).P(:,curr_site.distance_el);
            curr_site.start_sample = sim_data.output_data(curr_site.domain_row).start_sample(curr_site.distance_el);
            curr_site.no_samps_in_beat = length(sim_data.output_data(curr_site.domain_row).P(:,curr_site.distance_el));
        
            % store data for this site
            eval(['site' num2str(site_no) ' = curr_site;'])
        end
        clear site_no
        path_len = site2.len_aorta_site - site1.len_aorta_site;
        
        % check to see if the length of each wave is similar
        lens = [length(site1.v), length(site2.v)];
        if abs(diff(lens))/min(lens) > 0.1
            fprintf(['Different length PWs for sim ' num2str(sim_no) ', ' curr_pwv_type])
            Foot_TT = nan;
        else
            
            % time-align waves
            if site2.start_sample < site1.start_sample
                % then the peripheral wave has been extracted from the previous simulated wave
                site2.start_sample = site2.start_sample + site2.no_samps_in_beat;
            end
            delay = site2.start_sample - site1.start_sample;
            site2.v = site2.v([end-delay:end,1:(end-delay-1)]);
            
            % check waves are same duration
            if length(site2.v) == length(site1.v)+1
                site2.v = site2.v(1:end-1);
            elseif length(site2.v) == length(site1.v)+2
                site2.v = site2.v(1:end-2);
            elseif length(site1.v) == length(site2.v)+1
                site1.v = site1.v(1:end-1);
            elseif length(site1.v) == length(site2.v)+2
                site1.v = site1.v(1:end-2);
            end
            if length(site2.v) ~= length(site1.v)
                fprintf(['Different length PWs for sim ' num2str(sim_no) ', ' curr_pwv_type])
                Foot_TT = nan;
            else
                % repeat waves
                temp = linspace(site1.v(1),site1.v(end),length(site1.v));
                site1.v = site1.v+site1.v(1)-temp(:);
                site1.v = repmat(site1.v, [5,1]);
                temp = linspace(site2.v(1),site2.v(end),length(site2.v));
                site2.v = site2.v+site2.v(1)-temp(:);
                site2.v = repmat(site2.v, [5,1]);
                
                %plot(site1.v), hold on, plot(site2.v), xlim([0 1000])
                %close all
                
                % calculate PWV
                algo_ind = 1; %   1- foot-to-foot algorithm;   4- least squares algorithm
                Foot_TT = TTAlgorithm([site1.v site2.v],sim_data.output_data(1).fs,algo_ind,TT_waveform_ind,1,0);
            end
            
        end
        curr_pwv = path_len./Foot_TT(1);
        
%         % setup for this pwv type
%         switch curr_pwv_type
%             case 'carotid_femoral'
% %                 % Previously worked:
% %                 wave1_domain_no = 15;
% %                 wave1_distance_prop = 1;
% %                 wave2_domain_no = 46;
% %                 wave2_distance_prop = 0;
% %                 len_aorta_carotid = sum(sim_data.input_data.sim_settings.network_spec.length([1,2,15]));
% %                 len_aorta_femoral = sum(sim_data.input_data.sim_settings.network_spec.length([1,2,14,18,27,28,35,37,39,41,42,44]));
%                 % In keeping with calculation of input parameters:
%                 wave1_domain_no = 15;
%                 wave1_distance_prop = props.carotid;
%                 wave2_domain_no = 46;
%                 wave2_distance_prop = 0.5;
%                 len_aorta_carotid = sum([sim_data.input_data.sim_settings.network_spec.length([1,2]); sim_data.input_data.sim_settings.network_spec.length(15)*props.carotid]);
%                 len_aorta_femoral = sum([sim_data.input_data.sim_settings.network_spec.length([1,2,14,18,27,28,35,37,39,41,42,44]); sim_data.input_data.sim_settings.network_spec.length(46)/2]);
%                 % Calculate path length:
%                 path_len = len_aorta_femoral-len_aorta_carotid;
%                 
%             case 'carotid_radial'
%                 wave1_domain_no = 15;
%                 wave1_distance_prop = props.carotid;
%                 wave2_domain_no = 22;
%                 wave2_distance_prop = 1;
%                 len_aorta_carotid = sum([sim_data.input_data.sim_settings.network_spec.length([1,2]); sim_data.input_data.sim_settings.network_spec.length(15)*props.carotid]);
%                 len_aorta_radial = sum(sim_data.input_data.sim_settings.network_spec.length([1,2,14,19,21,22]));
%                 path_len = len_aorta_radial-len_aorta_carotid;
%                 
%             case 'brachial_radial'
%                 wave1_domain_no = 21;
%                 wave1_distance_prop = 0.75;
%                 wave2_domain_no = 22;
%                 wave2_distance_prop = 1;
%                 len_aorta_radial = sum(sim_data.input_data.sim_settings.network_spec.length([1,2,14,19,21,22]));
%                 len_aorta_brachial = sum(sim_data.input_data.sim_settings.network_spec.length([1,2,14,19])) + ...
%                     (sim_data.input_data.sim_settings.network_spec.length(21)*wave1_distance_prop);
%                 path_len = len_aorta_radial-len_aorta_brachial;
%                 
%             case 'brachial_ankle'
%                 wave1_domain_no = 21;
%                 wave1_distance_prop = 0.75;
%                 wave2_domain_no = 49;
%                 wave2_distance_prop = 1;
%                 len_aorta_brachial = sum(sim_data.input_data.sim_settings.network_spec.length([1,2,14,19])) + ...
%                     (sim_data.input_data.sim_settings.network_spec.length(21)*wave1_distance_prop);
%                 len_aorta_ankle = sum(sim_data.input_data.sim_settings.network_spec.length( [1,2,14,18,27,28,35,37,39,41,42,44,46,49] ));
%                 path_len = len_aorta_ankle-len_aorta_brachial;
%                 
%             case 'femoral_ankle'
%                 wave1_domain_no = 46;
%                 wave1_distance_prop = 0.5;
%                 wave2_domain_no = 49;
%                 wave2_distance_prop = 1;
%                 len_aorta_femoral = sum([sim_data.input_data.sim_settings.network_spec.length([1,2,14,18,27,28,35,37,39,41,42,44]); sim_data.input_data.sim_settings.network_spec.length(46)/2]);
%                 len_aorta_ankle = sum(sim_data.input_data.sim_settings.network_spec.length( [1,2,14,18,27,28,35,37,39,41,42,44,46,49] ));
%                 path_len = len_aorta_ankle-len_aorta_femoral;
%                 
%             case 'aortic'
%                 wave1_domain_no = 1;
%                 wave1_distance_prop = 0;
%                 wave2_domain_no = 41;
%                 wave2_distance_prop = 1;
%                 path_len = sum(sim_data.input_data.sim_settings.network_spec.length([1,2,14,18,27,28,35,37,39,41]));
%         end
%         
%         % extract required data
%         domain_row = find(sim_data.domain_nos == wave1_domain_no);
%         len = sim_data.input_data.sim_settings.network_spec.length(wave1_domain_no);
%         [~, wave1_distance_el] = min(abs(sim_data.output_data(domain_row).distances-(wave1_distance_prop*len)));
%         wave1.v = sim_data.output_data(domain_row).P(:,wave1_distance_el);
%         %wave1.no_samps_in_beat = length(sim_data.output_data(domain_row).P(:,wave1_distance_el));
%         %wave1.start_sample = rem(sim_data.output_data(domain_row).start_sample(wave1_distance_el), wave1.no_samps_in_beat);
%         wave1.start_sample = sim_data.output_data(domain_row).start_sample(wave1_distance_el);
%         %if wave1.start_sample > 0.5*wave1.no_samps_in_beat
%         %    wave1.start_sample = wave1.no_samps_in_beat - wave1.start_sample;
%         %end
%         domain_row = find(sim_data.domain_nos == wave2_domain_no);
%         len = sim_data.input_data.sim_settings.network_spec.length(wave2_domain_no);
%         [~, wave2_distance_el] = min(abs(sim_data.output_data(domain_row).distances-(wave2_distance_prop*len)));
%         wave2.v = sim_data.output_data(domain_row).P(:,wave2_distance_el);
%         wave2.no_samps_in_beat = length(sim_data.output_data(domain_row).P(:,wave2_distance_el));
%         %wave2.start_sample = rem(sim_data.output_data(domain_row).start_sample(wave2_distance_el), wave2.no_samps_in_beat);
%         wave2.start_sample = sim_data.output_data(domain_row).start_sample(wave2_distance_el);
%         %if wave2.start_sample > 0.5*wave2.no_samps_in_beat
%         %    wave2.start_sample = wave2.no_samps_in_beat - wave2.start_sample;
%         %end
%         clear domain_row
%         
%         % check to see if the length of each wave is similar
%         lens = [length(wave1.v), length(wave2.v)];
%         if abs(diff(lens))/min(lens) > 0.1
%             fprintf(['Different length PWs for sim ' num2str(sim_no) ', ' curr_pwv_type])
%             Foot_TT = nan;
%         else
%             
%             % time-align waves
%             if wave2.start_sample < wave1.start_sample
%                 % then the peripheral wave has been extracted from the previous simulated wave
%                 wave2.start_sample = wave2.start_sample + wave2.no_samps_in_beat;
%             end
%             delay = wave2.start_sample - wave1.start_sample;
%             wave2.v = wave2.v([end-delay:end,1:(end-delay-1)]);
%             
%             % check waves are same duration
%             if length(wave2.v) == length(wave1.v)+1
%                 wave2.v = wave2.v(1:end-1);
%             elseif length(wave2.v) == length(wave1.v)+2
%                 wave2.v = wave2.v(1:end-2);
%             elseif length(wave1.v) == length(wave2.v)+1
%                 wave1.v = wave1.v(1:end-1);
%             elseif length(wave1.v) == length(wave2.v)+2
%                 wave1.v = wave1.v(1:end-2);
%             end
%             if length(wave2.v) ~= length(wave1.v)
%                 fprintf(['Different length PWs for sim ' num2str(sim_no) ', ' curr_pwv_type])
%                 Foot_TT = nan;
%             else
%                 % repeat waves
%                 temp = linspace(wave1.v(1),wave1.v(end),length(wave1.v));
%                 wave1.v = wave1.v+wave1.v(1)-temp(:);
%                 wave1.v = repmat(wave1.v, [5,1]);
%                 temp = linspace(wave2.v(1),wave2.v(end),length(wave2.v));
%                 wave2.v = wave2.v+wave2.v(1)-temp(:);
%                 wave2.v = repmat(wave2.v, [5,1]);
%                 
%                 %plot(wave1.v), hold on, plot(wave2.v), xlim([0 1000])
%                 %close all
%                 
%                 % calculate PWV
%                 algo_ind = 1; %   1- foot-to-foot algorithm;   4- least squares algorithm
%                 Foot_TT = TTAlgorithm([wave1.v wave2.v],sim_data.output_data(1).fs,algo_ind,TT_waveform_ind,1,0);
%             end
%             
%         end
%         
% %         if strcmp(curr_pwv_type, 'brachial_radial')
% %             algo_ind = 1; % foot-to-foot algorithm
% %             Foot_TT_old = TTAlgorithm([wave1.v wave2.v],sim_data.output_data(1).fs,algo_ind,TT_waveform_ind,1,0);
% %             100*(Foot_TT_old(1)-Foot_TT(1))/Foot_TT(1)
% %             100*((delay/1000)-Foot_TT(1))/Foot_TT(1)
% %         end
% 
%         curr_pwv = path_len./Foot_TT(1);
        
        % store this PWV
        eval(['pulse_wave_vels(sim_no).pwv.' curr_pwv_type ' = curr_pwv;']);
        eval(['pulse_wave_vels(sim_no).ptt.' curr_pwv_type ' = Foot_TT(1);']);
        
        % calculate theoretical PWV
        [theor_pwv, theor_ptt] = calculate_theoretical_pwv(site1, site2, sim_data);
        eval(['pulse_wave_vels(sim_no).pwv_theor.' curr_pwv_type ' = theor_pwv;']);
        eval(['pulse_wave_vels(sim_no).ptt_theor.' curr_pwv_type ' = theor_ptt;']);
        
        clear wave1* wave2* theor_pwv theor_ptt curr_pwv Foot_TT site1 site2
    end
    clear curr_pwv_type pwv_type_no
end

% save
save(PATHS.pulse_wave_vels, 'pulse_wave_vels')

end

function [theor_pwv, theor_ptt] = calculate_theoretical_pwv(site1, site2, sim_data)

n_mini_segs = 10;
k = sim_data.input_data.sim_settings.network_spec.k;
rho = sim_data.input_data.sim_settings.rho;
domain_nos = extractfield(sim_data.output_data, 'domain_no');

for site_no = 1 :2
    
    % Extract data for this site
    eval(['curr_site = site' num2str(site_no) ';'])
    
    % Extract information on segments leading up to this site
    % - these are the prescribed values
    segs.in_rad = sim_data.input_data.sim_settings.network_spec.inlet_radius([curr_site.preceeding_domains, curr_site.domain_no]);
    segs.out_rad = sim_data.input_data.sim_settings.network_spec.outlet_radius([curr_site.preceeding_domains, curr_site.domain_no]);
    temp = segs.in_rad(end) - (segs.in_rad(end) - segs.out_rad(end))*curr_site.distance_prop;
    segs.out_rad(end) = temp; clear temp
    clear segs

    % - these are the simulated values (assuming linear tapering again)
    no_segs = length(curr_site.preceeding_domains)+1;
    for seg_no = 1 : no_segs
        if seg_no < no_segs
            curr_seg = curr_site.preceeding_domains(seg_no);
            temp_dist_prop = 1;
        else
            curr_seg = curr_site.domain_no;
            temp_dist_prop = curr_site.distance_prop;
        end
        rel_row = domain_nos == curr_seg;
        dists = sim_data.output_data(rel_row).distances;
        rel_dists = dists(dists<= sim_data.input_data.sim_settings.network_spec.length(curr_site.domain_no)*temp_dist_prop);
        rel_cols = 1:length(rel_dists);
        segs.in_rad(seg_no) = sqrt(mean(sim_data.output_data(rel_row).A(:,rel_cols(1)))/pi);
        segs.out_rad(seg_no) = sqrt(mean(sim_data.output_data(rel_row).A(:,rel_cols(end)))/pi);
    end
    
    % - lengths
    segs.len = [sim_data.input_data.sim_settings.network_spec.length(curr_site.preceeding_domains); sim_data.input_data.sim_settings.network_spec.length(curr_site.domain_no)*curr_site.distance_prop];
    
    % calculate theoretical ptts for each segment, and ptt for this site
    for seg_no = 1 : length(segs.in_rad)
        temp = linspace(segs.in_rad(seg_no), segs.out_rad(seg_no), n_mini_segs+1);
        mini_segs.in_rad = temp(1:end-1);
        mini_segs.out_rad = temp(2:end);
        mini_segs.len = ones(n_mini_segs,1)*segs.len(seg_no)/n_mini_segs; 
        mini_segs.mean_rad = mean([mini_segs.in_rad(:), mini_segs.out_rad(:)],2);
        mini_segs.wave_speed = empirical_wave_speed(mini_segs.mean_rad, k, rho);
        mini_segs.ptt = mini_segs.len./mini_segs.wave_speed;
        segs.ptt(seg_no) = sum(mini_segs.ptt);
        clear mini_segs temp
    end    
    clear seg_no
    
    curr_site.theor_ptt = sum(segs.ptt);
    clear segs
    
    % Store data for this site
    eval(['site' num2str(site_no) ' = curr_site;'])
    
end
clear site_no

% calculate theoretical pwv for this pair of sites
theor_ptt = site2.theor_ptt - site1.theor_ptt;
overall_len = site2.len_aorta_site - site1.len_aorta_site;
theor_pwv = overall_len/theor_ptt;

end

function wave_speed = empirical_wave_speed(ave_radius, k, rho)

ave_radius_cm = ave_radius*100;

Eh_D0 = (k(1)*exp(k(2)*ave_radius_cm))+k(3); % Eh/D0 (from Mynard's 2015 paper, eqn 3)
c0_squared = (2/3)*(Eh_D0/(rho/1000)); % from Mynard's 2015 paper, eqn 3, converts rho from kg/m3 to g/cm3
wave_speed = sqrt(c0_squared)/100;  % converts from cm/s to m/s.

end

function extract_pulse_wave_inds(PATHS, props, collated_data)

%% extract pulse indices
fprintf('\n - pulse indices: ')

% setup
options.do_plot = 0;
do_req_domains = 1;
if do_req_domains
    sites = {'AorticRoot',   'ThorAorta',  'AbdAorta',  'IliacBif', 'Carotid',  'SupTemporal', 'SupMidCerebral', 'Brachial', 'Radial', 'Digital', 'CommonIliac', 'Femoral', 'AntTibial'};
    site_domain_no = [1,            18,          39,          41,        15,             87,                 72,         21,       22,       112,             44,        46,           49];
    site_dist_prop = [0,             1,           0,        1,props.carotid,              1,                  1,       0.75,        1,         1,            0.5,       0.5,            1];
end

for sim_no = 1 : length(collated_data)
    fprintf([num2str(sim_no) ', ']);
    sim_data.domain_nos = extractfield(collated_data(sim_no).output_data, 'domain_no');
    sim_data.output_data = collated_data(sim_no).output_data;
    sim_data.input_data = collated_data(sim_no).input_data;
    
    % for each domain
    for domain_no_el = 1 : length(sim_data.domain_nos)
        
        curr_domain_no = sim_data.domain_nos(domain_no_el);
        if do_req_domains && ~sum(curr_domain_no == site_domain_no)
            continue
        end
        
        curr_distance_els = sim_data.output_data(domain_no_el).distances;
        no_distance_els = length(curr_distance_els);
        
        if do_req_domains
            dist_props = sim_data.output_data(domain_no_el).distances./sim_data.input_data.sim_settings.network_spec.length(curr_domain_no);
            curr_site_el = find(site_domain_no == curr_domain_no);
            [~,req_distance_el] = min(abs(dist_props-site_dist_prop(curr_site_el)));
        end
        
        for distance_el = 1 : no_distance_els
            
            if do_req_domains & distance_el ~= req_distance_el
                continue
            end
            
            % For Pressure
            sig.v = sim_data.output_data(domain_no_el).P(:,distance_el);
            sig.fs = sim_data.output_data(domain_no_el).fs;
            options.do_plot = 0; sig.ht = 1.75;
            if domain_no_el == 7, options.do_plot = 0; end
            [cv_inds, fid_pts, ~, ~] = PulseAnalyse10(sig, options);
            [cv_inds_new, cv_inds_names] = convert_var_to_min_struct(cv_inds);
            [fid_pts_new, fid_pts_names] = convert_var_to_min_struct(fid_pts);
            clear cv_inds fid_pts
            if options.do_plot
                close all
                plot(sig.v), hold on
                el = fid_pts_new(strcmp(fid_pts_names,'p1'));
                plot(el, sig.v(el), 'or')
                el = fid_pts_new(strcmp(fid_pts_names,'P1in'));
                plot(el, sig.v(el), '*r')
                el = fid_pts_new(strcmp(fid_pts_names,'p2'));
                plot(el, sig.v(el), 'ok')
                el = fid_pts_new(strcmp(fid_pts_names,'p2in'));
                plot(el, sig.v(el), '*k')
                close all
            end
            pulse_wave_inds(sim_no).P_pwa(domain_no_el).cv_inds(:,distance_el) = cv_inds_new;
            pulse_wave_inds(sim_no).P_pwa(domain_no_el).fid_pts(:,distance_el) = fid_pts_new;
            clear cv_inds_new fid_pts_new
            
            % For PPG
            sig.v = sim_data.output_data(domain_no_el).PPG(:,distance_el).^2;
            sig.fs = sim_data.output_data(domain_no_el).fs;
            options.do_plot = 0; sig.ht = 1.75;
            [cv_inds, fid_pts, ~, ~] = PulseAnalyse10(sig, options);
            [cv_inds_new, cv_inds_names] = convert_var_to_min_struct(cv_inds);
            [fid_pts_new, fid_pts_names] = convert_var_to_min_struct(fid_pts);
            clear cv_inds fid_pts
            if options.do_plot
                close all
            end
            pulse_wave_inds(sim_no).PPG_pwa(domain_no_el).cv_inds(:,distance_el) = cv_inds_new;
            pulse_wave_inds(sim_no).PPG_pwa(domain_no_el).fid_pts(:,distance_el) = fid_pts_new;
            
            clear fid_pts_new cv_inds_new sig
            
        end
        
        clear no_distance_els distance_el
        
    end
    clear domain_no_el
    
    % Add in names of variables
    if ~sum(strcmp(fieldnames(pulse_wave_inds(1)), 'cv_ind_names'))
        pulse_wave_inds(1).cv_ind_names = cv_inds_names;
        pulse_wave_inds(1).fid_pt_names = fid_pts_names;
    else
        a = [cv_inds_names, pulse_wave_inds(1).cv_ind_names];
        if ~isequal(a(:,1), a(:,2))
            error('Check this')
        end
        a = [fid_pts_names, pulse_wave_inds(1).fid_pt_names];
        if ~isequal(a(:,1), a(:,2))
            error('Check this')
        end 
        clear a
    end
    
    % Time to reflected wave
    rel_domain_el = sim_data.domain_nos == 15;
    distance_el = 1;
    rel_row = find(strcmp(pulse_wave_inds(1).fid_pt_names, 'p1'));
    pulse_wave_inds(sim_no).Tr = pulse_wave_inds(sim_no).P_pwa(rel_domain_el).fid_pts(rel_row,distance_el)/collated_data(sim_no).output_data(1).fs;  % as used in McEniery2005, and initially reported in Murgo1980
    
    % Max diastolic - min systolic flow at finger
    rel_domain_el = sim_data.domain_nos == 112;
    sig.v = sim_data.output_data(rel_domain_el).U(:,end);
    trs = find(sig.v(1:end-2) > sig.v(2:end-1) & sig.v(3:end) > sig.v(2:end-1));
    trs = trs(trs<collated_data(sim_no).input_data.sim_settings.lvet);
    [~,temp] = min(sig.v(trs));
    rel_tr = trs(temp); clear temp trs
    pks = find(sig.v(1:end-2) < sig.v(2:end-1) & sig.v(3:end) < sig.v(2:end-1));
    pks = pks(pks>collated_data(sim_no).input_data.sim_settings.lvet);
    [~, temp] = max(sig.v(pks));
    rel_pk = pks(temp); clear temp pks
    pulse_wave_inds(sim_no).Q_mm = sig.v(rel_pk) - sig.v(rel_tr);
    
    % Max diastolic - min systolic ppg at finger
    rel_domain_el = sim_data.domain_nos == 112;
    sig.v = sim_data.output_data(rel_domain_el).PPG(:,end);
    trs = find(sig.v(1:end-2) > sig.v(2:end-1) & sig.v(3:end) > sig.v(2:end-1));
    trs = trs(trs<collated_data(sim_no).input_data.sim_settings.lvet);
    [~,temp] = min(sig.v(trs));
    rel_tr = trs(temp); clear temp trs
    pks = find(sig.v(1:end-2) < sig.v(2:end-1) & sig.v(3:end) < sig.v(2:end-1));
    pks = pks(pks>collated_data(sim_no).input_data.sim_settings.lvet);
    [~, temp] = max(sig.v(pks));
    rel_pk = pks(temp); clear temp pks
    pulse_wave_inds(sim_no).Q_mm = sig.v(rel_pk) - sig.v(rel_tr); clear rel_tr rel_pk
    clear rel_domain_el sig sim_data
    
    % remove output data for this simulation to save memory
    collated_data(sim_no).input_data= [];
    collated_data(sim_no).output_data= [];
    collated_data(sim_no).system_chars= [];
end
clear sim_no

% save params
a = whos('pulse_wave_inds');
if a.bytes > 1.8e9
    save(PATHS.pulse_wave_inds, 'pulse_wave_inds', '-v7.3')
else
    save(PATHS.pulse_wave_inds, 'pulse_wave_inds')
end

end

function [new_var,new_names] = convert_var_to_min_struct(old_struct)

vars = fieldnames(old_struct);

% see if there are ".v" fields in each variable
eval(['temp = old_struct.' vars{1} ';'])
if isstruct(temp) & sum(strcmp(fieldnames(temp), 'v'))
    nested_v_fields = true;
else
    nested_v_fields = false;
end
clear temp

% extract data
new_var = nan(length(vars),1);
new_names = cell(length(vars),1);
for var_no = 1 : length(vars)
    curr_var = vars{var_no};
    if nested_v_fields
        eval(['curr_val = old_struct.' curr_var '.v;'])
    else
        eval(['curr_val = old_struct.' curr_var ';'])
    end
    new_var(var_no) = curr_val;
    new_names{var_no} = curr_var;    
    
end

end

function extract_haemodynamic_parameters(PATHS, props, collated_data)

fprintf('\n - Extracting haemodynamic parameters')

% load collated data
load(PATHS.pulse_wave_vels)
load(PATHS.pulse_wave_inds)

% setup
params = {'age',...
    'HR', 'SV', 'CO', 'LVET', 'dPdt', 'PFT', 'RFV', ... % cardiac
    'SBP_a', 'DBP_a', 'MBP_a', 'PP_a', 'SBP_b', 'DBP_b', 'MBP_b', 'PP_b', ... % routine BPs
    'SBP_f', 'DBP_f', 'MBP_f', 'PP_f', 'PP_amp', 'MBP_drop_finger', 'MBP_drop_ankle', 'SBP_diff', 'SMBP_a', ... % additional BPs
    'P1pk_a', 'P1pk_c', 'P1pk_b', 'P1pk_r', 'P1pk_d', 'P1in_a', 'P1in_c', 'P1in_b', 'P1in_r', 'P1in_d', ... % pulse wave points
    'P2pk_a', 'P2pk_c', 'P2pk_b', 'P2pk_r', 'P2pk_d', 'P2in_a', 'P2in_c', 'P2in_b', 'P2in_r', 'P2in_d', ... % pulse wave points
    'Ps_a', 'Ps_c', 'Ps_b', 'Ps_r', 'Ps_d', 'Pst_a', 'Pst_c', 'Pst_b', 'Pst_r', 'Pst_d', ... % pulse wave points
    'AI_a', 'AI_c', 'AI_b', 'AI_r', 'AI_d', 'AP_a', 'AP_c', 'AP_b', 'AP_r', 'AP_d', ... % augmentation indices
    'Tr_a', 'Tr_c', 'IAD', 'IADabs', 'svr', ... % non-routine parameters
    'PWV_a', 'PWV_cf', 'PWV_cr', 'PWV_ca', 'PWV_cb', 'PWV_bf', 'PWV_br', 'PWV_ba', 'PWV_rf', 'PWV_ra', 'PWV_fa',  ... % pulse wave velocities (measured)
    'PWVt_a', 'PWVt_cf', 'PWVt_cr', 'PWVt_ca', 'PWVt_cb', 'PWVt_bf', 'PWVt_br', 'PWVt_ba', 'PWVt_rf', 'PWVt_ra', 'PWVt_fa',  ... % pulse wave velocities (measured)
    'dia_asc_a', 'dia_desc_thor_a', 'dia_abd_a', 'dia_car', 'len_prox_a', ... % geometry
    'RI', 'SI', 'AGI_mod'}; % PPG indices
if exist(PATHS.system_chars, 'file')
    load(PATHS.system_chars)
    params = [params, {'pvr', 'pvc', 'pvc_iw', 'ac', 'c', 'tau'}];   % system chars
end
    %'PPamp_s_a', 'PPamp_P1_a', 'PPamp_P2_a', 'PPamp_s_c', 'PPamp_P1_c', 'PPamp_P2_c'};
for param_no = 1 : length(params)
    eval(['haemodynamic_params.' params{param_no} ' = nan(length(collated_data),1);']);
end
prop_brach = 3/4;   % distance along "brachial artery" (which is actually axillary and brachial in one long segment) to take cuff measurements

fid_pt_names = pulse_wave_inds(1).fid_pt_names;
cv_ind_names = pulse_wave_inds(1).cv_ind_names;

% extract parameters for each simulation
for sim_no = 1 : length(collated_data)
    sim_data.input_data = collated_data(sim_no).input_data;
    sim_data.output_data = collated_data(sim_no).output_data;
    sim_data.inds = pulse_wave_inds(sim_no);
    sim_data.vels = pulse_wave_vels(sim_no);
    sim_data.domain_nos = extractfield(sim_data.output_data, 'domain_no');
    
    % extract each parameter from this simulation in turn
    for param_no = 1 : length(params)
        
        % extract the parameter
        curr_param = params{param_no};
        
        if strcmp(curr_param, 'age')
            curr_val = sim_data.input_data.sim_settings.age;
        elseif (length(curr_param)>1 && sum(strcmp(curr_param(1:2), {'P1', 'P2', 'Ps', 'AI', 'AP'}))) ...
                || (length(curr_param) > 6 && strcmp(curr_param(1:6), 'PPamp_'))
            temp.sep = strfind(curr_param, '_');
            temp.quantity = curr_param(1:temp.sep(1)-1);
            if strcmp(temp.quantity(end), 't')
                temp.type = 'time';
                temp.quantity = temp.quantity(1:end-1);
            else
                temp.type = 'value';
            end
            temp.location = curr_param(temp.sep(end)+1:end);
            switch temp.location
                case 'a'
                    rel_domain_no = 1;
                    rel_dist_prop = 0;
                case 'c'
                    rel_domain_no = 15;
                    rel_dist_prop = props.carotid;
                case 'd'
                    rel_domain_no = 112;
                    rel_dist_prop = 1;
                case 'r'
                    rel_domain_no = 22;
                    rel_dist_prop = 1;
                case 'b'
                    rel_domain_no = 21;
                    rel_dist_prop = 0.75;
            end
            rel_row = find(sim_data.domain_nos == rel_domain_no);
            temp.distances = sim_data.output_data(rel_row).distances;
            [~, temp.rel_distance_el] = min(abs(temp.distances - (rel_dist_prop*sim_data.input_data.sim_settings.network_spec.length(rel_domain_no))));
            rel_sig = sim_data.output_data(rel_row).P(:,temp.rel_distance_el);
            %curr_val = 60/(length(rel_sig)/sim_data.output_data(rel_row).fs); % in bpm
            
            switch temp.quantity
                case 'P1pk'
                    rel_val_row = find(strcmp(fid_pt_names, 'p1pk'));
                    rel_pt.el = sim_data.inds.P_pwa(rel_row).fid_pts(rel_val_row, temp.rel_distance_el);
                case 'P2pk'
                    rel_val_row = find(strcmp(fid_pt_names, 'p2pk'));
                    rel_pt.el = sim_data.inds.P_pwa(rel_row).fid_pts(rel_val_row, temp.rel_distance_el);
                case 'P1in'
                    rel_val_row = find(strcmp(fid_pt_names, 'p1in'));
                    rel_pt.el = sim_data.inds.P_pwa(rel_row).fid_pts(rel_val_row, temp.rel_distance_el);
                case 'P2in'
                    rel_val_row = find(strcmp(fid_pt_names, 'p2in'));
                    rel_pt.el = sim_data.inds.P_pwa(rel_row).fid_pts(rel_val_row, temp.rel_distance_el);
                case 'Ps'
                    rel_val_row = find(strcmp(fid_pt_names, 's'));
                    rel_pt.el = sim_data.inds.P_pwa(rel_row).fid_pts(rel_val_row, temp.rel_distance_el);
                case 'AI'
                    rel_val_row1 = find(strcmp(fid_pt_names, 'p1in'));
                    rel_val_row2 = find(strcmp(fid_pt_names, 'p2pk'));
                    el1 = sim_data.inds.P_pwa(rel_row).fid_pts(rel_val_row1, temp.rel_distance_el);
                    el2 = sim_data.inds.P_pwa(rel_row).fid_pts(rel_val_row2, temp.rel_distance_el);
                    curr_wav = sim_data.output_data(rel_row).P(:,temp.rel_distance_el);
                    if isnan(el1)
                        p1 = nan;
                    else
                        p1 = curr_wav(el1);
                    end
                    if isnan(el2)
                        p2 = nan;
                    else
                        p2 = curr_wav(el2);
                    end
                    if ~isnan(p2) & ~isnan(p1)
                        rel_pt.v = 100*(p2-p1)/(max(curr_wav)-min(curr_wav)); % as percent
                    else
                        rel_pt.v = nan;
                    end
                    clear rel_val_row1 rel_val_row2 el1 el2 curr_wav p1 p2
                case 'AP'
                    rel_val_row1 = find(strcmp(fid_pt_names, 'p1in'));
                    rel_val_row2 = find(strcmp(fid_pt_names, 'p2pk'));
                    el1 = sim_data.inds.P_pwa(rel_row).fid_pts(rel_val_row1, temp.rel_distance_el);
                    el2 = sim_data.inds.P_pwa(rel_row).fid_pts(rel_val_row2, temp.rel_distance_el);
                    curr_wav = sim_data.output_data(rel_row).P(:,temp.rel_distance_el);
                    if isnan(el1)
                        p1 = nan;
                    else
                        p1 = curr_wav(el1);
                    end
                    if isnan(el2)
                        p2 = nan;
                    else
                        p2 = curr_wav(el2);
                    end
                    if ~isnan(p2) & ~isnan(p1)
                        rel_pt.v = (p2-p1)/133.33; % in mmHg
                    else
                        rel_pt.v = nan;
                    end
                    clear rel_val_row1 rel_val_row2 el1 el2 curr_wav p1 p2
                case 'PP_amp'
                    error('Check this')
            end
            if ~sum(strcmp(fieldnames(rel_pt),'v'))
                if isnan(rel_pt.el)
                    rel_pt.t = nan;
                    rel_pt.v = nan;
                else
                    rel_pt.t = (rel_pt.el-1)/sim_data.output_data(rel_row).fs;
                    rel_pt.v = rel_sig(rel_pt.el)/133.33;
                end
            end
            
            switch temp.type
                case 'time'
                    curr_val = rel_pt.t;
                case 'value'
                    curr_val = rel_pt.v;
            end
            clear temp rel_pt
            
        else
            switch curr_param
                case 'HR'
                    rel_domain_no = 1;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).P(:,1);
                    curr_val = 60/(length(rel_sig)/sim_data.output_data(rel_row).fs); % in bpm
                case 'SV'
                    rel_domain_no = 1;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).U(:,1).*sim_data.output_data(rel_row).A(:,1);
                    curr_val = sum(rel_sig)/sim_data.output_data(rel_row).fs; % in m3
                    curr_val = 1000*1000*curr_val;  % in ml
                case 'CO'
                    curr_val = haemodynamic_params(sim_no).HR * haemodynamic_params(sim_no).SV / 1000; % in l/min
                case 'LVET'
                    rel_domain_no = 1;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).U(:,1);
                    temp = find(rel_sig(11:end)<0,1); temp = temp+10-1;
                    curr_val = 1000*temp/sim_data.output_data(rel_row).fs; % in ms
                    clear temp
                case 'dPdt'
                    rel_domain_no = 1;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).P(:,1);
                    dPdt = (diff(rel_sig)/133.33).*(sim_data.output_data(rel_row).fs); % dP/dt in mmHg/s
                    curr_val = max(dPdt); clear dPdt
                case 'PFT'
                    rel_domain_no = 1;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).U(:,1).*sim_data.output_data(rel_row).A(:,1);
                    [~, temp_el] = max(rel_sig);
                    curr_val = 1000*temp_el/sim_data.output_data(rel_row).fs;  % in ms
                case 'RFV'
                    rel_domain_no = 1;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).U(:,1).*sim_data.output_data(rel_row).A(:,1);
                    rel_els = rel_sig<0;
                    curr_val = -1*sum(rel_sig(rel_els))/sim_data.output_data(rel_row).fs; % in m3
                    curr_val = 1000*1000*curr_val;  % in ml
                case 'SBP_a'
                    rel_domain_no = 1;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).P(:,1);
                    curr_val = max(rel_sig)/133.33;  % convert Pa to mmHg
                case 'DBP_a'
                    rel_domain_no = 1;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).P(:,1);
                    curr_val = min(rel_sig)/133.33;  % convert Pa to mmHg
                case 'MBP_a'
                    rel_domain_no = 1;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).P(:,1);
                    curr_val = mean(rel_sig)/133.33;  % convert Pa to mmHg
                case 'PP_a'
                    rel_domain_no = 1;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).P(:,1);
                    curr_val = range(rel_sig)/133.33;  % convert Pa to mmHg
                case 'SBP_b'
                    rel_domain_no = 21;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    temp.distances = sim_data.output_data(rel_row).distances;
                    [~, temp.rel_distance_el] = min(abs(temp.distances - (prop_brach*sim_data.input_data.sim_settings.network_spec.length(rel_domain_no))));
                    rel_sig = sim_data.output_data(rel_row).P(:,temp.rel_distance_el); clear temp
                    curr_val = max(rel_sig)/133.33;  % convert Pa to mmHg
                case 'DBP_b'
                    rel_domain_no = 21;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    temp.distances = sim_data.output_data(rel_row).distances;
                    [~, temp.rel_distance_el] = min(abs(temp.distances - (prop_brach*sim_data.input_data.sim_settings.network_spec.length(rel_domain_no))));
                    rel_sig = sim_data.output_data(rel_row).P(:,temp.rel_distance_el); clear temp
                    curr_val = min(rel_sig)/133.33;  % convert Pa to mmHg
                case 'MBP_b'
                    rel_domain_no = 21;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    temp.distances = sim_data.output_data(rel_row).distances;
                    [~, temp.rel_distance_el] = min(abs(temp.distances - (prop_brach*sim_data.input_data.sim_settings.network_spec.length(rel_domain_no))));
                    rel_sig = sim_data.output_data(rel_row).P(:,temp.rel_distance_el); clear temp
                    curr_val = mean(rel_sig)/133.33;  % convert Pa to mmHg
                case 'PP_b'
                    rel_domain_no = 21;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    temp.distances = sim_data.output_data(rel_row).distances;
                    [~, temp.rel_distance_el] = min(abs(temp.distances - (prop_brach*sim_data.input_data.sim_settings.network_spec.length(rel_domain_no))));
                    rel_sig = sim_data.output_data(rel_row).P(:,temp.rel_distance_el); clear temp
                    curr_val = range(rel_sig)/133.33;  % convert Pa to mmHg
                case 'SBP_f'
                    rel_domain_no = 112;  % digital
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).P(:,end);
                    curr_val = max(rel_sig)/133.33;  % convert Pa to mmHg
                case 'DBP_f'
                    rel_domain_no = 112;  % digital
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).P(:,end);
                    curr_val = min(rel_sig)/133.33;  % convert Pa to mmHg
                case 'MBP_f'
                    rel_domain_no = 112;  % digital
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).P(:,end);
                    curr_val = mean(rel_sig)/133.33;  % convert Pa to mmHg
                case 'PP_f'
                    rel_domain_no = 112;  % digital
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).P(:,end);
                    curr_val = range(rel_sig)/133.33;  % convert Pa to mmHg
                case 'PP_amp'
                    curr_val = haemodynamic_params(sim_no).PP_b / haemodynamic_params(sim_no).PP_a;
                case 'MBP_drop_finger'
                    rel_domain_no = 1;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).P(:,1);
                    central_mbp = mean(rel_sig)/133.33;  % convert Pa to mmHg
                    rel_domain_no = 112;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).P(:,end);
                    peripheral_mbp = mean(rel_sig)/133.33;  % convert Pa to mmHg
                    curr_val = central_mbp - peripheral_mbp;
                case 'MBP_drop_ankle'
                    rel_domain_no = 1;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).P(:,1);
                    central_mbp = mean(rel_sig)/133.33;  % convert Pa to mmHg
                    rel_domain_no = 49;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).P(:,end);
                    peripheral_mbp = mean(rel_sig)/133.33;  % convert Pa to mmHg
                    curr_val = central_mbp - peripheral_mbp;
                case 'SBP_diff'
                    rel_domain_no = 21;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    temp.distances = sim_data.output_data(rel_row).distances;
                    [~, temp.rel_distance_el] = min(abs(temp.distances - (prop_brach*sim_data.input_data.sim_settings.network_spec.length(rel_domain_no))));
                    rel_sig = sim_data.output_data(rel_row).P(:,temp.rel_distance_el); clear temp
                    temp.sbp_b = max(rel_sig)/133.33;
                    rel_domain_no = 1;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).P(:,1);
                    temp.sbp_a = max(rel_sig)/133.33;
                    curr_val = temp.sbp_b - temp.sbp_a;  % convert Pa to mmHg
                case 'SMBP_a'
                    rel_domain_no = 1;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_sig = sim_data.output_data(rel_row).P(:,1);
                    temp.sbp = max(rel_sig)/133.33;
                    temp.mbp = mean(rel_sig)/133.33;
                    curr_val = temp.sbp - temp.mbp;  % convert Pa to mmHg
                case 'IAD'
                    rel_domain_no = 7;  % right brachial
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    temp.distances = sim_data.output_data(rel_row).distances;
                    [~, distance_el] = min(abs(temp.distances - (prop_brach*sim_data.input_data.sim_settings.network_spec.length(rel_domain_no))));
                    rel_sig = sim_data.output_data(rel_row).P(:,distance_el);
                    right_arm_sbp = max(rel_sig)/133.33;  % convert Pa to mmHg
                    actual_prop_dist = temp.distances(distance_el)/sim_data.input_data.sim_settings.network_spec.length(rel_domain_no);
                    rel_domain_no = 21;  % left brachial
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    temp.distances = sim_data.output_data(rel_row).distances;
                    [~, distance_el] = min(abs(temp.distances - (actual_prop_dist*sim_data.input_data.sim_settings.network_spec.length(rel_domain_no))));
                    rel_sig = sim_data.output_data(rel_row).P(:,distance_el);
                    left_arm_sbp = max(rel_sig)/133.33;  % convert Pa to mmHg
                    curr_val = right_arm_sbp - left_arm_sbp;
                case 'IADabs'
                    rel_domain_no = 7;  % right brachial
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    temp.distances = sim_data.output_data(rel_row).distances;
                    [~, distance_el] = min(abs(temp.distances - (prop_brach*sim_data.input_data.sim_settings.network_spec.length(rel_domain_no))));
                    rel_sig = sim_data.output_data(rel_row).P(:,distance_el);
                    right_arm_sbp = max(rel_sig)/133.33;  % convert Pa to mmHg
                    actual_prop_dist = temp.distances(distance_el)/sim_data.input_data.sim_settings.network_spec.length(rel_domain_no);
                    rel_domain_no = 21;  % left brachial
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    temp.distances = sim_data.output_data(rel_row).distances;
                    [~, distance_el] = min(abs(temp.distances - (actual_prop_dist*sim_data.input_data.sim_settings.network_spec.length(rel_domain_no))));
                    rel_sig = sim_data.output_data(rel_row).P(:,distance_el);
                    left_arm_sbp = max(rel_sig)/133.33;  % convert Pa to mmHg
                    curr_val = abs(right_arm_sbp - left_arm_sbp);
                case 'Tr_c'
                    % - carotid
                    rel_domain_no = 15; % carotid
                    rel_prop = props.carotid;
                    % working
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    temp.distances = sim_data.output_data(rel_row).distances;
                    [~, distance_el] = min(abs(temp.distances - (rel_prop*sim_data.input_data.sim_settings.network_spec.length(rel_domain_no))));
                    rel_val_row = find(strcmp(fid_pt_names, 'p1in'));
                    curr_val = 1000*(sim_data.inds.P_pwa(rel_row).fid_pts(rel_val_row, distance_el)-1)/sim_data.output_data(1).fs;  % as used in McEniery2005, and initially reported in Murgo1980
                case 'Tr_a'
                    % - aortic root
                    rel_domain_no = 1;
                    rel_prop = 0;
                    % working
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    temp.distances = sim_data.output_data(rel_row).distances;
                    [~, distance_el] = min(abs(temp.distances - (rel_prop*sim_data.input_data.sim_settings.network_spec.length(rel_domain_no))));
                    rel_val_row = find(strcmp(fid_pt_names, 'p1in'));
                    curr_val = 1000*(sim_data.inds.P_pwa(rel_row).fid_pts(rel_val_row, distance_el)-1)/sim_data.output_data(1).fs;  % as used in McEniery2005, and initially reported in Murgo1980
                case 'PWV_a'
                    curr_val = sim_data.vels.pwv.aorta_iliacbif;
                case 'PWV_cf'
                    curr_val = sim_data.vels.pwv.carotid_femoral;
                case 'PWV_cr'
                    curr_val = sim_data.vels.pwv.carotid_radial;
                case 'PWV_ca'
                    curr_val = sim_data.vels.pwv.carotid_ankle;
                case 'PWV_cb'
                    curr_val = sim_data.vels.pwv.carotid_brachial;
                case 'PWV_bf'
                    curr_val = sim_data.vels.pwv.brachial_femoral;
                case 'PWV_br'
                    curr_val = sim_data.vels.pwv.brachial_radial;
                case 'PWV_ba'
                    curr_val = sim_data.vels.pwv.brachial_ankle;
                case 'PWV_rf'
                    curr_val = sim_data.vels.pwv.radial_femoral;
                case 'PWV_ra'
                    curr_val = sim_data.vels.pwv.radial_ankle;
                case 'PWV_fa'
                    curr_val = sim_data.vels.pwv.femoral_ankle;
                case 'PWVt_a'
                    curr_val = sim_data.vels.pwv_theor.aorta_iliacbif;
                case 'PWVt_cf'
                    curr_val = sim_data.vels.pwv_theor.carotid_femoral;
                case 'PWVt_cr'
                    curr_val = sim_data.vels.pwv_theor.carotid_radial;
                case 'PWVt_ca'
                    curr_val = sim_data.vels.pwv_theor.carotid_ankle;
                case 'PWVt_cb'
                    curr_val = sim_data.vels.pwv_theor.carotid_brachial;
                case 'PWVt_bf'
                    curr_val = sim_data.vels.pwv_theor.brachial_femoral;
                case 'PWVt_br'
                    curr_val = sim_data.vels.pwv_theor.brachial_radial;
                case 'PWVt_ba'
                    curr_val = sim_data.vels.pwv_theor.brachial_ankle;
                case 'PWVt_rf'
                    curr_val = sim_data.vels.pwv_theor.radial_femoral;
                case 'PWVt_ra'
                    curr_val = sim_data.vels.pwv_theor.radial_ankle;
                case 'PWVt_fa'
                    curr_val = sim_data.vels.pwv_theor.femoral_ankle;
                    
                    
                    % Diameters using values measured from the simulated waves
                case 'dia_asc_a'
                    % - asc aorta
                    rel_domain_nos = find(sim_data.input_data.sim_settings.network_spec.asc_aorta);
                    [~,domain_els,~] = intersect(sim_data.domain_nos, rel_domain_nos);
                    for s = 1 : length(domain_els)
                        curr_domain_el = domain_els(s);
                        ave_areas(s,1) = mean([max(sim_data.output_data(curr_domain_el).A(:,1)), max(sim_data.output_data(curr_domain_el).A(:,end))]);
                    end
                    lengths = sim_data.input_data.sim_settings.network_spec.length(rel_domain_nos);
                    curr_val = 1000*2*(sum(lengths.*sqrt(ave_areas./pi)))/sum(lengths);
                    clear ave_areas lengths
                case 'dia_desc_thor_a'
                    % - desc thor aorta
                    rel_domain_nos = find(sim_data.input_data.sim_settings.network_spec.desc_thor_aorta);
                    [~,domain_els,~] = intersect(sim_data.domain_nos, rel_domain_nos);
                    for s = 1 : length(domain_els)
                        curr_domain_el = domain_els(s);
                        ave_areas(s,1) = mean([max(sim_data.output_data(curr_domain_el).A(:,1)), max(sim_data.output_data(curr_domain_el).A(:,end))]);
                    end
                    lengths = sim_data.input_data.sim_settings.network_spec.length(rel_domain_nos);
                    curr_val = 1000*2*(sum(lengths.*sqrt(ave_areas./pi)))/sum(lengths);
                    clear ave_areas lengths
                case 'dia_abd_a'
                    % - abd aorta
                    rel_domain_nos = find(sim_data.input_data.sim_settings.network_spec.abd_aorta);
                    [~,domain_els,~] = intersect(sim_data.domain_nos, rel_domain_nos);
                    for s = 1 : length(domain_els)
                        curr_domain_el = domain_els(s);
                        ave_areas(s,1) = mean([max(sim_data.output_data(curr_domain_el).A(:,1)), max(sim_data.output_data(curr_domain_el).A(:,end))]);
                    end
                    lengths = sim_data.input_data.sim_settings.network_spec.length(rel_domain_nos);
                    curr_val = 1000*2*(sum(lengths.*sqrt(ave_areas./pi)))/sum(lengths);
                    clear ave_areas lengths
                case 'dia_car'
                    % - carotid
                    rel_domain_nos = find(sim_data.input_data.sim_settings.network_spec.both_carotid);
                    [~,domain_els,~] = intersect(sim_data.domain_nos, rel_domain_nos);
                    for s = 1 : length(domain_els)
                        curr_domain_el = domain_els(s);
                        ave_areas(s,1) = mean([mean(sim_data.output_data(curr_domain_el).A(:,1)), mean(sim_data.output_data(curr_domain_el).A(:,end))]);
                    end
                    lengths = sim_data.input_data.sim_settings.network_spec.length(rel_domain_nos);
                    curr_val = 1000*2*(sum(lengths.*sqrt(ave_areas./pi)))/sum(lengths);
                    clear ave_areas lengths
                case 'len_prox_a'
                    % - proximal aortic length
                    rel_els = sim_data.input_data.sim_settings.network_spec.proximal_aorta;
                    curr_val = 1000*sum( sim_data.input_data.sim_settings.network_spec.length(rel_els) );
                case 'svr'
                    % vascular
                    rel_domain_no = 1;
                    rel_row = find(sim_data.domain_nos == rel_domain_no);
                    rel_dist_el = 1;
                    curr_val = mean(sim_data.output_data(rel_row).P(:,rel_dist_el))/mean(sim_data.output_data(rel_row).U(:,rel_dist_el).*sim_data.output_data(rel_row).A(:,rel_dist_el))/(1e6);
                    % PPG
                case 'RI'
                    rel_domain_no = 112; domain_el = find(sim_data.domain_nos == rel_domain_no);
                    rel_dist_el = 3;
                    rel_val_row = find(strcmp(cv_ind_names, 'RI'));
                    curr_val = sim_data.inds.PPG_pwa(domain_el).cv_inds(rel_val_row, rel_dist_el);
                case 'SI'
                    rel_domain_no = 112; domain_el = find(sim_data.domain_nos == rel_domain_no);
                    rel_dist_el = 3;
                    rel_val_row = find(strcmp(cv_ind_names, 'delta_t'));
                    curr_val = 1.75/sim_data.inds.PPG_pwa(domain_el).cv_inds(rel_val_row, rel_dist_el);
                case 'AGI_mod'
                    rel_domain_no = 112; domain_el = find(sim_data.domain_nos == rel_domain_no);
                    rel_dist_el = 3;
                    rel_val_row = find(strcmp(cv_ind_names, 'AGI_mod'));
                    curr_val = sim_data.inds.PPG_pwa(domain_el).cv_inds(rel_val_row, rel_dist_el);
                case 'pvr'
                    curr_val = system_chars.pvr(sim_no);
                case 'pvc'
                    curr_val = system_chars.pvc(sim_no);
                case 'pvc_iw'
                    curr_val = system_chars.pvc_iw(sim_no);
                case 'ac'
                    curr_val = system_chars.ac(sim_no);
                case 'c'
                    curr_val = system_chars.c(sim_no);
                case 'tau'
                    curr_val = system_chars.tau(sim_no);
            end
            clear temp
        end
        
        % store this extracted parameter
        eval(['haemodynamic_params(sim_no).' curr_param ' = curr_val;'])
        
        clear rel_dist_el curr_val rel_domain_no domain_el rel_val_row
    end
    
end

% save params
save(PATHS.haemodynamic_params, 'haemodynamic_params')

end
