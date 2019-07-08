function PATHS = setup_paths_for_post_processing(pwdb_no)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      This line specifies the location of the storage folder      %
%         which the input and output files are copied into         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PATHS.storage_folder = '/Users/petercharlton/Documents/Data/Nektar1D/ageing_sims/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      This line specifies the location of the shared folder       %
%         containing the simulation input and output files         %
%       (this can be left alone unless reproducing the PWDB)       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PATHS.shared_folder = '/Users/petercharlton/Documents/VM-share/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   The rest of this function can be left alone   %%%%%%%%

% setup
if nargin>0
    fprintf(['  (PWDB no. ' num2str(pwdb_no), ')\n'])
end
close all
    
PATHS.slash_dirn = filesep;

% add additional functions to path
[a,~,~] = fileparts(mfilename('fullpath'));
addpath(genpath(a))

if nargin > 0
    % Specify paths
    PATHS.pwdb_filename_prefix = 'pwdb';   % for use in creating exported data files.
    PATHS.OutputFiles = [PATHS.storage_folder, 'pwdb_' num2str(pwdb_no), PATHS.slash_dirn, 'output_data', PATHS.slash_dirn];
    PATHS.InputFiles = [PATHS.storage_folder, 'pwdb_' num2str(pwdb_no), PATHS.slash_dirn, 'input_data', PATHS.slash_dirn];
    PATHS.ProcessedData = [PATHS.storage_folder, 'pwdb_' num2str(pwdb_no), PATHS.slash_dirn, 'processed_data', PATHS.slash_dirn];
    PATHS.exported_data = [PATHS.storage_folder, 'pwdb_' num2str(pwdb_no), filesep, 'exported_data', filesep];
    PATHS.CaseStudies = [PATHS.storage_folder, 'pwdb_' num2str(pwdb_no), filesep, 'case_studies', filesep];
    PATHS.Analysis = [PATHS.storage_folder, 'pwdb_' num2str(pwdb_no), filesep, 'analysis', filesep];
    PATHS.Analysis_tables = [PATHS.storage_folder, 'pwdb_' num2str(pwdb_no), filesep, 'analysis', filesep, 'tables', filesep];
    PATHS.Analysis_figures = [PATHS.storage_folder, 'pwdb_' num2str(pwdb_no), filesep, 'analysis', filesep, 'figures', filesep];
    PATHS.history_files_data = [PATHS.ProcessedData, 'history_files_data.mat'];
    PATHS.collated_data = [PATHS.ProcessedData, 'collated_data.mat'];
    PATHS.haemodynamic_params = [PATHS.ProcessedData, 'haemodynamic_params.mat'];
    PATHS.pulse_wave_params = [PATHS.ProcessedData, 'pulse_wave_params.mat'];
    PATHS.pulse_wave_vels = [PATHS.ProcessedData, 'pulse_wave_vels.mat'];
    PATHS.pulse_wave_inds = [PATHS.ProcessedData, 'pulse_wave_inds.mat'];
    PATHS.table_data = [PATHS.Analysis_tables, 'table_data.mat'];
    PATHS.baseline_waves_fig_a = [PATHS.Analysis_figures, 'baseline_waves_fig_a'];
    PATHS.baseline_waves_fig_b = [PATHS.Analysis_figures, 'baseline_waves_fig_b'];
    PATHS.pw_propagation_fig = [PATHS.Analysis_figures, 'pw_propagation_fig'];
    PATHS.pw_path_fig = [PATHS.Analysis_figures, 'pw_path_fig'];
    PATHS.wave_speed_fig = [PATHS.Analysis_figures, 'wave_speed_fig'];
    PATHS.characteristics_table = [PATHS.Analysis_tables, 'characteristics_table'];
    PATHS.characteristics_table_range = [PATHS.Analysis_tables, 'characteristics_table_range'];
    PATHS.system_chars = [PATHS.ProcessedData, 'system_chars.mat'];
    PATHS.exported_data_mat = [PATHS.exported_data, 'PWs', filesep, 'mat', filesep];
    PATHS.exported_data_csv = [PATHS.exported_data, 'PWs', filesep, 'csv', filesep];
    PATHS.exported_data_wfdb = [PATHS.exported_data, 'PWs', filesep, 'wfdb', filesep];
    PATHS.exported_data_geo = [PATHS.exported_data, 'geo', filesep];
    PATHS.exported_data_mat_pwdb_data = [PATHS.exported_data, PATHS.pwdb_filename_prefix, '_data.mat'];
    PATHS.exported_data_mat_pwdb_data_w_aorta_finger_path = [PATHS.exported_data, PATHS.pwdb_filename_prefix, '_data_w_aorta_finger_path.mat'];
    PATHS.exported_data_mat_pwdb_data_w_aorta_foot_path_p = [PATHS.exported_data, PATHS.pwdb_filename_prefix, '_data_w_aorta_foot_path_p.mat'];
    PATHS.exported_data_mat_pwdb_data_w_aorta_foot_path_u = [PATHS.exported_data, PATHS.pwdb_filename_prefix, '_data_w_aorta_foot_path_u.mat'];
    PATHS.exported_data_mat_pwdb_data_w_aorta_foot_path_a = [PATHS.exported_data, PATHS.pwdb_filename_prefix, '_data_w_aorta_foot_path_a.mat'];
    PATHS.exported_data_mat_pwdb_data_w_aorta_brain_path = [PATHS.exported_data, PATHS.pwdb_filename_prefix, '_data_w_aorta_brain_path.mat'];
    PATHS.exported_data_mat_pwdb_data_w_aorta_rsubclavian_path = [PATHS.exported_data, PATHS.pwdb_filename_prefix, '_data_w_aorta_rsubclavian_path.mat'];
    PATHS.exported_data_model_configs = [PATHS.exported_data, PATHS.pwdb_filename_prefix, '_model_configs.csv'];
    PATHS.peripheral_boundarys = [PATHS.ProcessedData, 'peripheral_boundarys.mat'];
    PATHS.exported_data_onset_times = [PATHS.exported_data, PATHS.pwdb_filename_prefix, '_onset_times.csv'];
    
    % case studies
    PATHS.collated_ASIs = [PATHS.CaseStudies, 'collated_ASIs.mat'];
    PATHS.ASI_results = [PATHS.CaseStudies, 'ASI_results.mat'];
    
    % Create required folders
    req_folders = {PATHS.storage_folder, PATHS.InputFiles, PATHS.OutputFiles, PATHS.ProcessedData, PATHS.CaseStudies};
    for folder_no = 1 : length(req_folders)
        curr_folder = req_folders{folder_no};
        if ~exist(curr_folder, 'dir')
            mkdir(curr_folder)
        end
        clear curr_folder
    end
    clear folder_no req_folders
end

end