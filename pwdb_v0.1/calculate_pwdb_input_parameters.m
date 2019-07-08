function parameters = calculate_pwdb_input_parameters
% CALCULATE_PWDB_INPUT_PARAMETERS calculates a set of Nektar1D input
% parameters for the virtual subjects in the pwdb Pulse Wave DataBase.
%
%               calculate_pwdb_input_parameters
%
%   Inputs:     none required - please see the "setup_up"
%               function below to adjust the virtual database configuration
%               to your needs.
%
%   Outputs:    - parameters: a Matlab variable containing the input
%                   parameters.
%               - a single file called "inputs.mat", which contains the
%                   parameters variable required by Nektar1D input files. This
%                   will be saved in the path specified in "setup_up", which
%                   should be adjusted before running this script.
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
% v.1.0  Peter H. Charlton, King's College London

%% Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTINGS TO CHANGE: This function specifies where to save the outputs,  %%%%
%%%%                 and which database configuration to use                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
up = setup_up;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify desired virtual database characteristics
vdb_characteristics = setup_vdb_characteristics(up);   % Specifies variations (in terms of the no. of SD) in each parameter at baseline age and other ages

%% Make plots of how parameters change with age and vary
if up.do_plot, make_param_plots(up); end

%% Calculate equations for input parameters as a function of age
eqns = calculate_equations(up);

%% Find variations in parameters for each simulation
variations = find_vdb_variations(vdb_characteristics);   % Specifies the variation of each parameter (in terms of SD) from the age-specific baseline value for each simulation

%% Find values of parameters for each simulation
parameters = calculate_vdb_parameters(variations, up);  % Calculates the desired model parameter values for each simulation

%% Add fixed Wk vascular bed parameters
parameters = add_wk_vascular_bed_parameters(parameters);

%% Set fixed simulation parameters
parameters = specify_fixed_simulation_parameters(parameters);

%% Make plots of expected haemodynamic characteristics
if up.do_plot, make_haemod_plots(parameters, up); end

%% Generate aortic inflow waveforms
inflow = generate_aortic_inflow_waveforms(parameters);

%% Save parameters to file
save([up.savefolder, 'inputs'], 'parameters', 'inflow')

end

function up = setup_up
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This line specifies the where the outputs are saved        %%%%
%%%%  - This should be changed for your computer                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

up.savefolder = '/Users/petercharlton/Documents/Data/Nektar1D/ageing_sims/';  % (needs to have a slash at the end)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% These lines specify which database configuration to use    %%%%
%%%%  - This should be changed according to your needs          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

possible_configs = {'pwdb', 'initial_simulations', 'independent_variations', 'baseline_subject', 'baseline_subjects_at_each_age'};
up.db_config = possible_configs{1};   % specifies which of the configurations to use.

% See this webpage for further details on each of the configurations:
% https://github.com/peterhcharlton/pwdb/wiki/Generating-the-Input-Files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% The rest of this function specifies some preliminary settings, which don't need to be changed %%%%

% Initial setup
fprintf(['\n\n ~~~~~ Calculating input parameters for ' upper(up.db_config) ' database ~~~~~'])
fprintf('\n - Setting up universal parameters')
close all

% add subdirectories within this file's directory to the path
[curr_folder, ~, ~] = fileparts(mfilename('fullpath'));
addpath(genpath(curr_folder));

% Specify all possible ages (only a subset of these will actually be used in the database):
up.mod_ages = 25:75;
% Specify the baseline age, which is the one the baseline model configuration corresponds to:
up.baseline_age = 25;

% Give details of the Excel file which contains the baseline model geometry:
up.network_spec_sheet = 'Haemod 116 segments';
curr_file_root = fileparts(mfilename('fullpath'));
up.network_spec_file = [curr_file_root, filesep, 'Input Data', filesep, '116_artery_model.txt'];

% Specify whether or not to make plots
up.do_plot = 0;

end

function vdb_characteristics = setup_vdb_characteristics(up)
% This function specifies each of the input parameters for the virtual
% database configuration selected in the up.db_config variable, in the
% setup_up function.

fprintf('\n - Setting up parameter variations across database subjects')

%% Ages
% Specify the age of the baseline subject
vdb_characteristics.age.baseline = up.baseline_age;

%% - Physiological Properties
%
% Specify for each property:
% - The number of values of this property to be changed:
%    (a) in combination with other properties (i.e. if both HR and SV are varied in combination, then all possible combinations of these properties will be simulated)
%    (b) independently of other properties (i.e. keep all properties except one at the mean value, whilst that property is varied on its own. Then repeat for other properties which are to be varied).
%   at:
%    (i) the baseline age (usually 25 years), and
%    (ii) other ages (usually 35, 45, 55, 65 and 75 years).
% - The range over which to vary the property, specified as the number of
%   SDs from the mean value.
%
% e.g. a value of "2" for the number of variations, and a value of "1" for
% the number of SDs, means that the parameter will be varied twice over the
% range of +/- 1 SD. i.e. it will take the values:
%     [mean - 1SD, mean, mean + 1SD]
%
% e.g. a value of "4" for the number of variations, and a value of "1" for
% the number of SDs, means that the parameter will be varied four times
% over the range of +/- 1 SD. i.e. it will take the values:
%     [mean - 1SD, mean - 0.5SD, mean, mean + 0.5SD, mean + 1SD]

switch up.db_config
    
    %% %% %% Initial Simulations %% %% %%
    % These consist of changing 10 parameters independently by +/- 1SD for the
    % 25-year old baseline model, resulting in a total of 21 simulations (a
    % baseline simulation, and then 2 simulations per parameter).
    case 'initial_simulations'
        
        vdb_characteristics.age.all = 25;
        
        %% Cardiac properties
        
        % - Heart Rate (HR)
        vdb_characteristics.hr.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.hr.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.hr.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.hr.no_independent_variations_at_baseline_age = 2;
        vdb_characteristics.hr.variation_sds = 1;
        
        % - Stroke Volume (SV)
        vdb_characteristics.sv.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.sv.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.sv.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.sv.no_independent_variations_at_baseline_age = 2;
        vdb_characteristics.sv.variation_sds = 1;
        
        % - Left Ventricular Ejection Time (LVET): duration of systole
        vdb_characteristics.lvet.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.lvet.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.lvet.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.lvet.no_independent_variations_at_baseline_age = 2;
        vdb_characteristics.lvet.variation_sds = 1;
        
        % - Time to Peak Flow (t_pf): the time from the beginning of the aortic flow wave to its peak
        vdb_characteristics.t_pf.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.t_pf.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.t_pf.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.t_pf.no_independent_variations_at_baseline_age = 2;
        vdb_characteristics.t_pf.variation_sds = 1;
        
        % - Regurgitation volume (reg_vol): the amount of backwards flow at the end of systole
        vdb_characteristics.reg_vol.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.reg_vol.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.reg_vol.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.reg_vol.no_independent_variations_at_baseline_age = 2;
        vdb_characteristics.reg_vol.variation_sds = 1;
        
        %% Vascular bed properties
        
        % - Peripheral Vascular Compliance (pvc)
        vdb_characteristics.pvc.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.pvc.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.pvc.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.pvc.no_independent_variations_at_baseline_age = 2;
        vdb_characteristics.pvc.variation_sds = 1;
        
        % - Outflow Pressure (p_out)
        vdb_characteristics.p_out.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.p_out.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.p_out.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.p_out.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.p_out.variation_sds = 0;
        
        % - Mean Blood Pressure (mbp): the peripheral vascular resistance is set to achieve the desired value of MBP
        vdb_characteristics.mbp.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.mbp.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.mbp.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.mbp.no_independent_variations_at_baseline_age = 2;
        vdb_characteristics.mbp.variation_sds = 1;
        
        % - Reflections logical: whether to include reflections
        vdb_characteristics.reflect_log.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.reflect_log.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.reflect_log.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.reflect_log.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.reflect_log.variation_sds = 0;
        
        %% Arterial properties
        
        % - Pulse Wave Velocity (pwv)
        vdb_characteristics.pwv.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.pwv.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.pwv.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.pwv.no_independent_variations_at_baseline_age = 2;
        vdb_characteristics.pwv.variation_sds = 1;
        
        % - Diastolic blood pressure (dbp): the pressure corresponding to the initial luminal areas (note that this has very little affect)
        vdb_characteristics.dbp.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.dbp.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.dbp.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.dbp.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.dbp.variation_sds = 0;
        
        % - Pressure drop (p_drop): the assumed pressure drop from the aortic root to the end of the arterial network
        vdb_characteristics.p_drop.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.p_drop.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.p_drop.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.p_drop.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.p_drop.variation_sds = 0;
        
        % - Length of the proximal aorta (len)
        vdb_characteristics.len.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.len.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.len.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.len.no_independent_variations_at_baseline_age = 2;
        vdb_characteristics.len.variation_sds = 1;
        
        % - Diameter of the larger arteries (dia, e.g. ascending, descending aorta, carotid, brachial, ...)
        vdb_characteristics.dia.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.dia.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.dia.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.dia.no_independent_variations_at_baseline_age = 2;
        vdb_characteristics.dia.variation_sds = 1;
        
        % - Wall viscosity
        vdb_characteristics.gamma.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.gamma.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.gamma.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.gamma.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.gamma.variation_sds = 0;
        
        %% Blood properties
        
        % - Blood density (rho)
        vdb_characteristics.rho.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.rho.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.rho.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.rho.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.rho.variation_sds = 0;
        
        % - Blood viscosity (mu)
        vdb_characteristics.mu.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.mu.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.mu.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.mu.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.mu.variation_sds = 0;
        
        %% - Simulation properties
        
        % - alpha: controls the shape of the velocity profile (assumed to be constant throughout the arterial tree)
        vdb_characteristics.alpha.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.alpha.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.alpha.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.alpha.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.alpha.variation_sds = 0;
        
        % - time_step: which is the numerical time step used for calculations
        vdb_characteristics.time_step.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.time_step.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.time_step.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.time_step.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.time_step.variation_sds = 0;
        
        % - visco_elastic_log: determines whether the simulations are elastic or visco-elastic
        vdb_characteristics.visco_elastic_log.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.visco_elastic_log.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.visco_elastic_log.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.visco_elastic_log.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.visco_elastic_log.variation_sds = 0;
        
        
    %% %% %% Pulse Wave Database %% %% %%
    % These consist of changing 6 parameters in combination with each other by
    % +/- 1SD for all ages (25, 35, 45, 55, 65 and 75), resulting in a total of
    % 4,374 simulations. This is the configuration used for the Pulse Wave
    % Database.
    case 'pwdb'
        
        vdb_characteristics.age.all = 25:10:75;  % Ten-year intervals
        
        %% Cardiac properties
        
        % - Heart Rate (HR)
        vdb_characteristics.hr.no_combination_variations_at_other_ages = 2;
        vdb_characteristics.hr.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.hr.no_combination_variations_at_baseline_age = 2;
        vdb_characteristics.hr.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.hr.variation_sds = 1;
        
        % - Stroke Volume (SV)
        vdb_characteristics.sv.no_combination_variations_at_other_ages = 2;
        vdb_characteristics.sv.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.sv.no_combination_variations_at_baseline_age = 2;
        vdb_characteristics.sv.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.sv.variation_sds = 1;
        
        % - Left Ventricular Ejection Time (LVET): duration of systole
        vdb_characteristics.lvet.no_combination_variations_at_other_ages = 2;
        vdb_characteristics.lvet.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.lvet.no_combination_variations_at_baseline_age = 2;
        vdb_characteristics.lvet.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.lvet.variation_sds = 1;
        
        % - Time to Peak Flow (t_pf): the time from the beginning of the aortic flow wave to its peak
        vdb_characteristics.t_pf.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.t_pf.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.t_pf.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.t_pf.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.t_pf.variation_sds = 1;
        
        % - Regurgitation volume (reg_vol): the amount of backwards flow at the end of systole
        vdb_characteristics.reg_vol.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.reg_vol.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.reg_vol.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.reg_vol.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.reg_vol.variation_sds = 1;
        
        % - Reflections logical: whether to include reflections
        vdb_characteristics.reflect_log.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.reflect_log.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.reflect_log.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.reflect_log.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.reflect_log.variation_sds = 0;
        
        %% Vascular bed properties
        
        % - Peripheral Vascular Compliance (pvc)
        vdb_characteristics.pvc.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.pvc.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.pvc.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.pvc.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.pvc.variation_sds = 1;
        
        % - Outflow Pressure (p_out)
        vdb_characteristics.p_out.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.p_out.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.p_out.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.p_out.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.p_out.variation_sds = 1;
        
        % - Mean Blood Pressure (mbp): the peripheral vascular resistance is set to achieve the desired value of MBP
        vdb_characteristics.mbp.no_combination_variations_at_other_ages = 2;
        vdb_characteristics.mbp.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.mbp.no_combination_variations_at_baseline_age = 2;
        vdb_characteristics.mbp.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.mbp.variation_sds = 1;
        
        %% Arterial properties
        
        % - Pulse Wave Velocity (pwv)
        vdb_characteristics.pwv.no_combination_variations_at_other_ages = 2;
        vdb_characteristics.pwv.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.pwv.no_combination_variations_at_baseline_age = 2;
        vdb_characteristics.pwv.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.pwv.variation_sds = 1;
        
        % - Diastolic blood pressure (dbp): the pressure corresponding to the initial luminal areas (note that this has very little affect)
        vdb_characteristics.dbp.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.dbp.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.dbp.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.dbp.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.dbp.variation_sds = 0;
        
        % - Pressure drop (p_drop): the assumed pressure drop from the aortic root to the end of the arterial network
        vdb_characteristics.p_drop.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.p_drop.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.p_drop.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.p_drop.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.p_drop.variation_sds = 0;
        
        % - Length of the proximal aorta (len)
        vdb_characteristics.len.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.len.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.len.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.len.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.len.variation_sds = 1;
        
        % - Diameter of the larger arteries (dia, e.g. ascending, descending aorta, carotid, brachial, ...)
        vdb_characteristics.dia.no_combination_variations_at_other_ages = 2;
        vdb_characteristics.dia.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.dia.no_combination_variations_at_baseline_age = 2;
        vdb_characteristics.dia.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.dia.variation_sds = 1;
        
        % - Wall viscosity
        vdb_characteristics.gamma.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.gamma.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.gamma.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.gamma.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.gamma.variation_sds = 0;
        
        %% Blood properties
        
        % - Blood density (rho)
        vdb_characteristics.rho.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.rho.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.rho.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.rho.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.rho.variation_sds = 0;
        
        % - Blood viscosity (mu)
        vdb_characteristics.mu.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.mu.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.mu.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.mu.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.mu.variation_sds = 0;
        
        %% - Simulation properties
        
        % - alpha: controls the shape of the velocity profile (assumed to be constant throughout the arterial tree)
        vdb_characteristics.alpha.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.alpha.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.alpha.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.alpha.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.alpha.variation_sds = 0;
        
        % - time_step: which is the numerical time step used for calculations
        vdb_characteristics.time_step.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.time_step.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.time_step.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.time_step.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.time_step.variation_sds = 0;
        
        % - visco_elastic_log: determines whether the simulations are elastic or visco-elastic
        vdb_characteristics.visco_elastic_log.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.visco_elastic_log.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.visco_elastic_log.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.visco_elastic_log.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.visco_elastic_log.variation_sds = 0;
        
    %% %% %% Independent variation in properties across ages %% %% %%
    % These consist of changing each parameter independently by
    % +/- 0.5 SD and +/- 1SD for all ages (25, 35, 45, 55, 65 and 75). This
    % was the approach used in the preliminary pulse wave database.
    case 'independent_variations'
        
        vdb_characteristics.age.all = 25:10:75;  % Ten-year intervals
        
        %% Cardiac properties
        
        % - Heart Rate (HR)
        vdb_characteristics.hr.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.hr.no_independent_variations_at_baseline_age = 4;
        vdb_characteristics.hr.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.hr.no_independent_variations_at_other_ages = 4;
        vdb_characteristics.hr.variation_sds = 1;
        
        % - Stroke Volume (SV)
        vdb_characteristics.sv.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.sv.no_independent_variations_at_other_ages = 4;
        vdb_characteristics.sv.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.sv.no_independent_variations_at_baseline_age = 4;
        vdb_characteristics.sv.variation_sds = 1;
        
        % - Left Ventricular Ejection Time (LVET): duration of systole
        vdb_characteristics.lvet.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.lvet.no_independent_variations_at_other_ages = 4;
        vdb_characteristics.lvet.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.lvet.no_independent_variations_at_baseline_age = 4;
        vdb_characteristics.lvet.variation_sds = 1;
        
        % - Time to Peak Flow (t_pf): the time from the beginning of the aortic flow wave to its peak
        vdb_characteristics.t_pf.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.t_pf.no_independent_variations_at_other_ages = 4;
        vdb_characteristics.t_pf.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.t_pf.no_independent_variations_at_baseline_age = 4;
        vdb_characteristics.t_pf.variation_sds = 1;
        
        % - Regurgitation volume (reg_vol): the amount of backwards flow at the end of systole
        vdb_characteristics.reg_vol.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.reg_vol.no_independent_variations_at_other_ages = 2;
        vdb_characteristics.reg_vol.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.reg_vol.no_independent_variations_at_baseline_age = 2;
        vdb_characteristics.reg_vol.variation_sds = 1;
        
        %% Vascular bed properties
        
        % - Peripheral Vascular Compliance (pvc)
        vdb_characteristics.pvc.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.pvc.no_independent_variations_at_other_ages = 2;
        vdb_characteristics.pvc.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.pvc.no_independent_variations_at_baseline_age = 2;
        vdb_characteristics.pvc.variation_sds = 1;
        
        % - Outflow Pressure (p_out)
        vdb_characteristics.p_out.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.p_out.no_independent_variations_at_other_ages = 2;
        vdb_characteristics.p_out.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.p_out.no_independent_variations_at_baseline_age = 2;
        vdb_characteristics.p_out.variation_sds = 1;
        
        % - Mean Blood Pressure (mbp): the peripheral vascular resistance is set to achieve the desired value of MBP
        vdb_characteristics.mbp.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.mbp.no_independent_variations_at_other_ages = 4;
        vdb_characteristics.mbp.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.mbp.no_independent_variations_at_baseline_age = 4;
        vdb_characteristics.mbp.variation_sds = 1;
        
        % - Reflections logical: whether to include reflections
        vdb_characteristics.reflect_log.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.reflect_log.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.reflect_log.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.reflect_log.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.reflect_log.variation_sds = 0;
        
        %% Arterial properties
        
        % - Pulse Wave Velocity (pwv)
        vdb_characteristics.pwv.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.pwv.no_independent_variations_at_other_ages = 4;
        vdb_characteristics.pwv.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.pwv.no_independent_variations_at_baseline_age = 4;
        vdb_characteristics.pwv.variation_sds = 1;
        
        % - Diastolic blood pressure (dbp): the pressure corresponding to the initial luminal areas (note that this has very little affect)
        vdb_characteristics.dbp.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.dbp.no_independent_variations_at_other_ages = 4;
        vdb_characteristics.dbp.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.dbp.no_independent_variations_at_baseline_age = 4;
        vdb_characteristics.dbp.variation_sds = 1;
        
        % - Pressure drop (p_drop): the assumed pressure drop from the aortic root to the end of the arterial network
        vdb_characteristics.p_drop.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.p_drop.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.p_drop.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.p_drop.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.p_drop.variation_sds = 0;
        
        % - Length of the proximal aorta (len)
        vdb_characteristics.len.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.len.no_independent_variations_at_other_ages = 4;
        vdb_characteristics.len.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.len.no_independent_variations_at_baseline_age = 4;
        vdb_characteristics.len.variation_sds = 1;
        
        % - Diameter of the larger arteries (dia, e.g. ascending, descending aorta, carotid, brachial, ...)
        vdb_characteristics.dia.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.dia.no_independent_variations_at_other_ages = 4;
        vdb_characteristics.dia.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.dia.no_independent_variations_at_baseline_age = 4;
        vdb_characteristics.dia.variation_sds = 1;
        
        % - Wall viscosity
        vdb_characteristics.gamma.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.gamma.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.gamma.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.gamma.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.gamma.variation_sds = 0;
        
        %% Blood properties
        
        % - Blood density (rho)
        vdb_characteristics.rho.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.rho.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.rho.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.rho.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.rho.variation_sds = 0;
        
        % - Blood viscosity (mu)
        vdb_characteristics.mu.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.mu.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.mu.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.mu.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.mu.variation_sds = 0;
        
        %% - Simulation properties
        
        % - alpha: controls the shape of the velocity profile (assumed to be constant throughout the arterial tree)
        vdb_characteristics.alpha.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.alpha.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.alpha.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.alpha.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.alpha.variation_sds = 0;
        
        % - time_step: which is the numerical time step used for calculations
        vdb_characteristics.time_step.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.time_step.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.time_step.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.time_step.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.time_step.variation_sds = 0;
        
        % - visco_elastic_log: determines whether the simulations are elastic or visco-elastic
        vdb_characteristics.visco_elastic_log.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.visco_elastic_log.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.visco_elastic_log.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.visco_elastic_log.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.visco_elastic_log.variation_sds = 0;
        
        %% %% %% Single baseline simulation %% %% %%
        % This is a single simulation of the baseline 25-year old model.
    case 'baseline_subject'
        
        vdb_characteristics.age.all = 25;
        
        %% Cardiac properties
        
        % - Heart Rate (HR)
        vdb_characteristics.hr.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.hr.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.hr.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.hr.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.hr.variation_sds = 0;
        
        % - Stroke Volume (SV)
        vdb_characteristics.sv.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.sv.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.sv.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.sv.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.sv.variation_sds = 0;
        
        % - Left Ventricular Ejection Time (LVET): duration of systole
        vdb_characteristics.lvet.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.lvet.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.lvet.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.lvet.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.lvet.variation_sds = 0;
        
        % - Time to Peak Flow (t_pf): the time from the beginning of the aortic flow wave to its peak
        vdb_characteristics.t_pf.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.t_pf.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.t_pf.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.t_pf.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.t_pf.variation_sds = 0;
        
        % - Regurgitation volume (reg_vol): the amount of backwards flow at the end of systole
        vdb_characteristics.reg_vol.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.reg_vol.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.reg_vol.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.reg_vol.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.reg_vol.variation_sds = 0;
        
        %% Vascular bed properties
        
        % - Peripheral Vascular Compliance (pvc)
        vdb_characteristics.pvc.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.pvc.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.pvc.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.pvc.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.pvc.variation_sds = 0;
        
        % - Outflow Pressure (p_out)
        vdb_characteristics.p_out.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.p_out.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.p_out.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.p_out.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.p_out.variation_sds = 0;
        
        % - Mean Blood Pressure (mbp): the peripheral vascular resistance is set to achieve the desired value of MBP
        vdb_characteristics.mbp.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.mbp.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.mbp.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.mbp.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.mbp.variation_sds = 0;
        
        % - Reflections logical: whether to include reflections
        vdb_characteristics.reflect_log.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.reflect_log.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.reflect_log.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.reflect_log.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.reflect_log.variation_sds = 0;
        
        %% Arterial properties
        
        % - Pulse Wave Velocity (pwv)
        vdb_characteristics.pwv.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.pwv.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.pwv.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.pwv.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.pwv.variation_sds = 0;
        
        % - Diastolic blood pressure (dbp): the pressure corresponding to the initial luminal areas (note that this has very little affect)
        vdb_characteristics.dbp.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.dbp.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.dbp.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.dbp.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.dbp.variation_sds = 0;
        
        % - Pressure drop (p_drop): the assumed pressure drop from the aortic root to the end of the arterial network
        vdb_characteristics.p_drop.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.p_drop.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.p_drop.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.p_drop.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.p_drop.variation_sds = 0;
        
        % - Length of the proximal aorta (len)
        vdb_characteristics.len.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.len.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.len.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.len.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.len.variation_sds = 0;
        
        % - Diameter of the larger arteries (dia, e.g. ascending, descending aorta, carotid, brachial, ...)
        vdb_characteristics.dia.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.dia.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.dia.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.dia.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.dia.variation_sds = 0;
        
        % - Wall viscosity
        vdb_characteristics.gamma.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.gamma.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.gamma.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.gamma.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.gamma.variation_sds = 0;
        
        %% Blood properties
        
        % - Blood density (rho)
        vdb_characteristics.rho.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.rho.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.rho.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.rho.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.rho.variation_sds = 0;
        
        % - Blood viscosity (mu)
        vdb_characteristics.mu.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.mu.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.mu.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.mu.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.mu.variation_sds = 0;
        
        %% - Simulation properties
        
        % - alpha: controls the shape of the velocity profile (assumed to be constant throughout the arterial tree)
        vdb_characteristics.alpha.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.alpha.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.alpha.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.alpha.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.alpha.variation_sds = 0;
        
        % - time_step: which is the numerical time step used for calculations
        vdb_characteristics.time_step.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.time_step.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.time_step.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.time_step.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.time_step.variation_sds = 0;
        
        % - visco_elastic_log: determines whether the simulations are elastic or visco-elastic
        vdb_characteristics.visco_elastic_log.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.visco_elastic_log.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.visco_elastic_log.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.visco_elastic_log.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.visco_elastic_log.variation_sds = 0;
        
        %% %% %% Baseline simulation at each age %% %% %%
        
    case 'baseline_subjects_at_each_age'
        
        vdb_characteristics.age.all = 25:10:75;
        
        %% Cardiac properties
        
        % - Heart Rate (HR)
        vdb_characteristics.hr.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.hr.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.hr.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.hr.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.hr.variation_sds = 0;
        
        % - Stroke Volume (SV)
        vdb_characteristics.sv.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.sv.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.sv.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.sv.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.sv.variation_sds = 0;
        
        % - Left Ventricular Ejection Time (LVET): duration of systole
        vdb_characteristics.lvet.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.lvet.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.lvet.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.lvet.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.lvet.variation_sds = 0;
        
        % - Time to Peak Flow (t_pf): the time from the beginning of the aortic flow wave to its peak
        vdb_characteristics.t_pf.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.t_pf.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.t_pf.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.t_pf.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.t_pf.variation_sds = 0;
        
        % - Regurgitation volume (reg_vol): the amount of backwards flow at the end of systole
        vdb_characteristics.reg_vol.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.reg_vol.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.reg_vol.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.reg_vol.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.reg_vol.variation_sds = 0;
        
        %% Vascular bed properties
        
        % - Peripheral Vascular Compliance (pvc)
        vdb_characteristics.pvc.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.pvc.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.pvc.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.pvc.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.pvc.variation_sds = 0;
        
        % - Outflow Pressure (p_out)
        vdb_characteristics.p_out.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.p_out.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.p_out.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.p_out.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.p_out.variation_sds = 0;
        
        % - Mean Blood Pressure (mbp): the peripheral vascular resistance is set to achieve the desired value of MBP
        vdb_characteristics.mbp.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.mbp.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.mbp.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.mbp.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.mbp.variation_sds = 0;
        
        % - Reflections logical: whether to include reflections
        vdb_characteristics.reflect_log.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.reflect_log.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.reflect_log.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.reflect_log.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.reflect_log.variation_sds = 0;
        
        %% Arterial properties
        
        % - Pulse Wave Velocity (pwv)
        vdb_characteristics.pwv.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.pwv.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.pwv.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.pwv.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.pwv.variation_sds = 0;
        
        % - Diastolic blood pressure (dbp): the pressure corresponding to the initial luminal areas (note that this has very little affect)
        vdb_characteristics.dbp.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.dbp.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.dbp.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.dbp.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.dbp.variation_sds = 0;
        
        % - Pressure drop (p_drop): the assumed pressure drop from the aortic root to the end of the arterial network
        vdb_characteristics.p_drop.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.p_drop.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.p_drop.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.p_drop.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.p_drop.variation_sds = 0;
        
        % - Length of the proximal aorta (len)
        vdb_characteristics.len.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.len.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.len.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.len.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.len.variation_sds = 0;
        
        % - Diameter of the larger arteries (dia, e.g. ascending, descending aorta, carotid, brachial, ...)
        vdb_characteristics.dia.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.dia.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.dia.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.dia.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.dia.variation_sds = 0;
        
        % - Wall viscosity
        vdb_characteristics.gamma.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.gamma.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.gamma.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.gamma.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.gamma.variation_sds = 0;
        
        %% Blood properties
        
        % - Blood density (rho)
        vdb_characteristics.rho.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.rho.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.rho.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.rho.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.rho.variation_sds = 0;
        
        % - Blood viscosity (mu)
        vdb_characteristics.mu.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.mu.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.mu.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.mu.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.mu.variation_sds = 0;
        
        %% - Simulation properties
        
        % - alpha: controls the shape of the velocity profile (assumed to be constant throughout the arterial tree)
        vdb_characteristics.alpha.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.alpha.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.alpha.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.alpha.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.alpha.variation_sds = 0;
        
        % - time_step: which is the numerical time step used for calculations
        vdb_characteristics.time_step.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.time_step.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.time_step.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.time_step.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.time_step.variation_sds = 0;
        
        % - visco_elastic_log: determines whether the simulations are elastic or visco-elastic
        vdb_characteristics.visco_elastic_log.no_combination_variations_at_other_ages = 0;
        vdb_characteristics.visco_elastic_log.no_independent_variations_at_other_ages = 0;
        vdb_characteristics.visco_elastic_log.no_combination_variations_at_baseline_age = 0;
        vdb_characteristics.visco_elastic_log.no_independent_variations_at_baseline_age = 0;
        vdb_characteristics.visco_elastic_log.variation_sds = 0;
        
end

% If the number of SDs of variation for any parameter which is to be varied
% is set at zero, then change it to a default value of 1.
params = fieldnames(vdb_characteristics); params = params(~strcmp(params, 'age'));
for param_no = 1 : length(params)
    curr_param = params{param_no};
    eval(['curr_data = vdb_characteristics.' curr_param ';']);
    if curr_data.variation_sds == 0 & ...
            sum([abs(curr_data.no_combination_variations_at_other_ages), ...
            abs(curr_data.no_independent_variations_at_other_ages), ...
            abs(curr_data.no_combination_variations_at_baseline_age), ...
            abs(curr_data.no_independent_variations_at_baseline_age)]) > 0
        eval(['vdb_characteristics.' curr_param '.variation_sds = 1;']);
    end
    
end

end

function variations = find_vdb_variations(vdb_characteristics)
% This function converts the database configuration (set in
% "setup_vdb_characteristics") into the number of SD to vary each parameter
% by in each simulation.

fprintf('\n - Finding variation values')

%% Find variations for each characteristic
% i.e. find out what values (in terms of SD from the mean) each parameter
% should take in the database, at what age(s), and whether to inclue these
% variations independently or in combination with each other.
varying_params = setxor(fieldnames(vdb_characteristics), 'age');
ages = vdb_characteristics.age;
for varying_param_no = 1 : length(varying_params)
    eval(['curr_varying_param = vdb_characteristics.' varying_params{varying_param_no} ';'])
    eval(['param_variations.' varying_params{varying_param_no} ' = find_param_values(curr_varying_param, ages);']);
end

%% Assemble list of variations for each simulation

[variations.params, variations.ages, variations.age_nos] = deal([]);
variations.param_names = varying_params;
for age_no = 1 : length(ages.all)
    curr_age = ages.all(age_no);
    
    % identify variations in each parameter (in terms of SD from the mean)
    % for this age, and categorise them according to whether these
    % variations are to be used independently, or in combination with each other.
    [ind_variations, comb_variations] = deal(cell(length(varying_params),1));
    for varying_param_no = 1 : length(varying_params)
        eval(['curr_variations = param_variations.' varying_params{varying_param_no} ';']);
        rel_variation_els = curr_variations.ages == curr_age;
        rel_variations = curr_variations.variations(rel_variation_els); rel_variations = rel_variations(:)';
        rel_comb_log = curr_variations.comb_log(rel_variation_els);
        
        ind_variations{varying_param_no} = rel_variations(~rel_comb_log);
        comb_variations{varying_param_no} = rel_variations(rel_comb_log);
        
        clear rel_comb_log rel_variations rel_variation_els curr_variations
    end
    clear varying_param_no
    
    % find independent variations
    all_ind_variations = [];
    for varying_param_no = 1 : length(ind_variations)
        rel_ind_variations = ind_variations{varying_param_no};
        temp = zeros(length(rel_ind_variations),length(varying_params));
        temp(:,varying_param_no) = rel_ind_variations(:);
        all_ind_variations = [all_ind_variations; temp];
        clear temp rel_ind_variations
    end
    clear varying_param_no
    % NB: all_ind_variations is a matrix specifying the variation (in terms
    % of SD) of each parameter (columns correspond to different parameters)
    % for each simulation (rows correspond to simulations).
    
    % find combination variations
    combvec_text = '';
    for s = 1 : length(comb_variations)
        combvec_text = [combvec_text, 'comb_variations{' num2str(s) '}, '];
    end
    combvec_text = combvec_text(1:end-2);
    eval(['all_comb_variations = combvec(' combvec_text ');']);
    all_comb_variations = all_comb_variations';
    % NB: all_comb_variations is a matrix specifying the variation (in terms
    % of SD) of each parameter (columns correspond to different parameters)
    % for each simulation (rows correspond to simulations).
    
    % extract baseline simulation
    baseline_el = sum(abs(all_comb_variations')) == 0;
    baseline_sim = all_comb_variations(baseline_el,:);
    all_comb_variations = all_comb_variations(~baseline_el,:);
    
    % sort variations to put independent ones immediately after baseline ones
    additional_ind_variations = sum(all_comb_variations==0,2) == size(all_comb_variations,2)-1;
    if ~isempty(all_comb_variations)
        all_ind_variations = [all_ind_variations; all_comb_variations(additional_ind_variations,:)];
    end
    all_ind_variations = sortrows(all_ind_variations);
    all_comb_variations = all_comb_variations(~additional_ind_variations,:);
    
    % list variations
    variations.params = [variations.params; baseline_sim; all_ind_variations; all_comb_variations];
    variations.ages = [variations.ages; curr_age*ones(1+size(all_comb_variations,1)+size(all_ind_variations,1),1)];
    variations.age_nos = [variations.age_nos; age_no*ones(1+size(all_comb_variations,1)+size(all_ind_variations,1),1)];
    
end

%% Re-order variations to put the priority ones first:
% Prioritise:
% (i) baseline for age
priority1 = sum(variations.params~=0,2) == 0;
% (ii) single parameter changed
priority2 = sum(variations.params~=0,2) == 1 & ~priority1;
% (iii) all parameters changed
priority3 = sum(variations.params~=0,2) == max(sum(variations.params~=0,2))  & ~priority1 & ~priority2;
% (iv) remainder
remainder = ~priority1 & ~priority2 & ~priority3;
% Re-order
new_order = [find(priority1); find(priority2); find(priority3); find(remainder)];
clear priority1 priority2 priority3 remainder
variations.params = variations.params(new_order,:);
variations.ages = variations.ages(new_order);
variations.age_nos = variations.age_nos(new_order);
clear new_order

%% Eliminate non-sensical variations of categorical options
% This are non-numerical parameters - i.e. logicals which indicate whether
% or not to model something (such as whether or not to model
% visco-elasticity in arterial walls).
cat_options = {'visco_elastic_log', 'reflect_log'};

for param_no = 1 : length(variations.param_names)
    
    % skip if this parameter isn't a categorical option
    if ~sum(strcmp(variations.param_names{param_no}, cat_options))
        continue
    end
    
    % eliminate non-sensical (i.e. negative) variations
    variations_to_keep = variations.params(:,param_no) == 0 | variations.params(:,param_no) == -1;
    variations.params = variations.params(variations_to_keep,:);
    variations.ages = variations.ages(variations_to_keep);
    variations.age_nos = variations.age_nos(variations_to_keep);
    clear variations_to_keep
    
end

end

function param_variations = find_param_values(curr_param, ages)
% This function works out what values (in terms of SD from the mean) each parameter
% should take in the database, at what age(s), and whether to inclue these
% variations independently or in combination with each other.

% cycle through each age
param_variations.ages = [];
param_variations.variations = [];
param_variations.comb_log = [];
param_variations.baseline_log = [];
    
for age_no = 1 : length(ages.all)
    curr_age = ages.all(age_no);
    
    % extract relevant info on desired variations
    if curr_age == ages.baseline
        no_combination_variations = curr_param.no_combination_variations_at_baseline_age;
        no_independent_variations = curr_param.no_independent_variations_at_baseline_age;
    else
        no_combination_variations = curr_param.no_combination_variations_at_other_ages;
        no_independent_variations = curr_param.no_independent_variations_at_other_ages;
    end
        
    % Find variations
    total_variations = no_combination_variations + no_independent_variations;
    if total_variations > 0
    all_variations = linspace(-1*curr_param.variation_sds,curr_param.variation_sds, total_variations+1);
    else
        all_variations = 0;
    end
    clear total_variations
    
    % Find which of these are in combination, and which are independent
    if no_combination_variations > 0
        comb_variations = linspace(-1*curr_param.variation_sds,curr_param.variation_sds, no_combination_variations+1);
    else
        comb_variations = [];
    end
    comb_variation_log = false(length(all_variations),1);
    for s = 1 : length(all_variations)
        if sum(all_variations(s) == comb_variations)
            comb_variation_log(s) = true;
        end
    end
    clear s comb_variations no_combination_variations no_independent_variations
    
    % find which are baseline
    baseline_variation_log = logical(all_variations' == 0);
    comb_variation_log = baseline_variation_log | comb_variation_log;
        
    % assemble matrix of desired variations
    param_variations.ages = [param_variations.ages; curr_age*ones(length(all_variations),1)]; clear curr_age
    param_variations.variations = [param_variations.variations; all_variations(:)]; clear all_variations
    param_variations.comb_log = logical([param_variations.comb_log; comb_variation_log]); clear comb_variation_log
    param_variations.baseline_log = logical([param_variations.baseline_log; baseline_variation_log]); clear baseline_variation_log
    
end

end

function parameters = calculate_vdb_parameters(variations, up)
% This function converts the list of variations of each parameter (in terms
% of SD from the mean) into absolute values which can be prescribed to the
% model. For instance, the "variations" variable might say that heart rate
% should be varied by [-1, 0, 1] SD from the mean at age 25. This function
% converts that information into actual heart rate values (such as [60, 70,
% 80] bpm, if the mean +/- SD HR at age 25 was 70+/-10 bpm).

fprintf('\n - Calculating model input parameters\n\n')

% setup
required_age = unique(variations.ages);
parameters.age = variations.ages;
parameters.variations = rmfield(variations, {'ages', 'age_nos'});

%% Calculate mean and SD values for parameters which either: (i) are independent of age, or (ii) are solely dependent on age
ind_params = {'t_pf', 'reg_vol', 'p_out', 'rho', 'mu', 'alpha', 'time_step', 'p_drop', 'visco_elastic_log', 'reflect_log'};
age_dep_params = {'dbp', 'pvc', 'dbp', 'mbp', 'hr', 'sv', 'len', 'dia_asc', 'dia_desc_thor', 'dia_abd', 'dia_carotid'};
params = [age_dep_params, ind_params];
eqns = struct;

% Calculate values at each age ...
for curr_age_no = 1 : length(required_age)
    curr_age = required_age(curr_age_no);
    
    % ... for each parameter
    for param_no = 1 : length(params)
        curr_param = params{param_no};
        
        % Find the mean and sd values for this parameter for this age
        [baseline_val, sd, eqn] = extract_values_from_literature(curr_age, curr_param, up);
        
        % store results
        eval(['phys_vals.' curr_param '.val(curr_age_no,1) = baseline_val;'])
        eval(['phys_vals.' curr_param '.sd(curr_age_no,1) = sd;'])
        
        % store equation
        eval(['eqns.' curr_param ' = eqn;'])
        
        clear baseline_val sd curr_param abbr eqn
    end
    clear param_no curr_age
end
clear curr_age_no params age_dep_params ind_params

%% Calculate values of parameters which are either: (i) independent of, or (ii) dependent on, the model baseline for each simulation
ind_params = {'hr', 'sv', 't_pf', 'reg_vol', 'dbp', 'mbp', 'mu', 'alpha', 'time_step', 'p_drop', 'visco_elastic_log', 'reflect_log'};
dep_params = {'pvc', 'p_out', 'rho'};
params = [ind_params, dep_params]; clear ind_params dep_params

% cycle through each parameter
for param_no = 1 : length(params)
    curr_param = params{param_no};
    curr_param_col = find(strcmp(variations.param_names, curr_param));
    
    % extract variations for this parameter for each simulation 
    param_variations.no_sd = variations.params(:,curr_param_col); clear curr_param_col
    param_variations.age_ind = variations.age_nos;
    
    % calculate values for this parameter for each simulation 
    eval(['curr_param_phys_vals = phys_vals.' curr_param ';'])
    param_values = curr_param_phys_vals.val(param_variations.age_ind) + (param_variations.no_sd .* curr_param_phys_vals.sd(param_variations.age_ind));
    clear param_variations
    
    % store these values
    eval(['parameters.' curr_param ' = param_values;'])
    
    clear curr_param
end
clear param_no params

%% Correct time step for visco-elastic simulations
% This ensures that the time step is short enough for visco-elastic simulations to run

max_visco_elastic_time_step  = 1e-5;
% identify visco-elastic simulations which have a time step above the maximum value
time_steps_to_change = parameters.time_step > max_visco_elastic_time_step & ...
    parameters.visco_elastic_log == 1;
% reset the time step of these simulations to 1e-5 (the max for visco-elastic)
parameters.time_step(time_steps_to_change) = max_visco_elastic_time_step;


%% Calculate parameters which are dependent on several initial model parameters

%% -- len --
% length of proximal aorta
curr_param = 'len';

% import geometry of arterial tree from data file
baseline_network_spec = obtain_network_spec(eqns, up);

% extract param variations for this parameter
curr_param_col = find(strcmp(variations.param_names, curr_param));
param_variations.no_sd = variations.params(:,curr_param_col);
param_variations.age_ind = variations.age_nos;

% Find desired proximal aortic lengths for each simulation
eval(['curr_param_phys_vals = phys_vals.' curr_param ';'])
temp_param_values = curr_param_phys_vals.val(param_variations.age_ind) + (param_variations.no_sd .* curr_param_phys_vals.sd(param_variations.age_ind));  % values in mm

% Scale these desired lengths so that the desired length at the baseline
% age is equal to the length in the initial model configuration.
param_baseline_val = 1000*sum(baseline_network_spec.length(baseline_network_spec.proximal_aorta)); % length of proximal aorta in mm
model_scale_factor = param_baseline_val/curr_param_phys_vals.val(required_age == up.baseline_age);
model_param_values = temp_param_values*model_scale_factor/1000; % the required lengths of the proximal aorta in each of the simulations, in m

% For each simulation, find the lengths of the individual proximal aortic
% segments which sum to give these scaled desired lengths.
baseline_proximal_aorta_lengths = baseline_network_spec.length(baseline_network_spec.proximal_aorta);
baseline_proximal_aorta_lengths = baseline_proximal_aorta_lengths(:)';
proximal_aorta_param_values = repmat(baseline_proximal_aorta_lengths, [length(model_param_values),1]) .* (model_param_values/sum(baseline_proximal_aorta_lengths));  % in m

% store these values in the network spec for each simulation
for sim_no = 1 : length(variations.ages)
    parameters.network_spec{sim_no} = baseline_network_spec;
    parameters.network_spec{sim_no}.length(baseline_network_spec.proximal_aorta) = proximal_aorta_param_values(sim_no,:);
end

%% -- dia --
% diameter of asc aorta, desc aorta, abd aorta, and common carotid artery
% (more or less - it's actually all segments with a radius of >= 0.003462 m) 
curr_param = 'dia';

% extract variations for this parameter
curr_param_col = find(strcmp(variations.param_names, curr_param));
param_variations.no_sd = variations.params(:,curr_param_col);
param_variations.age_ind = variations.age_nos;

% find overall percentage change in aortic diameter across all four sites for each age
baseline_dia_age_no = 1;
curr_artery_dia = phys_vals.dia_asc.val; perc_change(:,1) = 100*curr_artery_dia/curr_artery_dia(baseline_dia_age_no);   % the values are for each baseline simulation for each age
curr_artery_dia = phys_vals.dia_desc_thor.val; perc_change(:,2) = 100*curr_artery_dia/curr_artery_dia(baseline_dia_age_no);
curr_artery_dia = phys_vals.dia_abd.val; perc_change(:,3) = 100*curr_artery_dia/curr_artery_dia(baseline_dia_age_no);
curr_artery_dia = phys_vals.dia_carotid.val; perc_change(:,4) = 100*curr_artery_dia/curr_artery_dia(baseline_dia_age_no);
curr_param_phys_vals.val = mean(perc_change,2);
clear curr_artery_dia perc_change

% To find eqn: perc_change = curr_param_phys_vals.val; age = 25:10:75; age = age(:); tbl = table(perc_change, age); mdl = fitlm(tbl, 'perc_change ~ age')

% Find the SD (as a proportion of the mean diameter)
sd_div_mean = mean([2.4/25.7,2.3/26.5,2.1/28.8,3.1/29.0]); % From Agmon 2003
% Find the SD (absolute value)
curr_param_phys_vals.sd = sd_div_mean*curr_param_phys_vals.val;

% find the desired diameter as a proportion of the baseline diameter for each simulation
prop_of_baseline = curr_param_phys_vals.val(param_variations.age_ind) + (param_variations.no_sd .* curr_param_phys_vals.sd(param_variations.age_ind));
prop_of_baseline = prop_of_baseline/100; % converts from percent to proportion

% store these values in the network spec for each simulation
av_radius = mean([parameters.network_spec{1}.inlet_radius, parameters.network_spec{1}.outlet_radius],2);
segs_to_increase = av_radius >= 0.003462*sqrt(1.5); % all segments with radius greater than 4.24mm (3.46 mm scaled by sqrt(1.5)) are increased in diameter
outlet_nodes_to_increase = parameters.network_spec{1}.outlet_node(segs_to_increase); % the outlet diameters for all these segments should be increased.
additional_inlets_to_increase = false(length(parameters.network_spec{1}.inlet_node),1); % the inlets of segments attached to the ends of the segments whose outlet diameters are increased should also be increased.
for s = 1 : length(parameters.network_spec{1}.inlet_node)
    curr_inlet_node = parameters.network_spec{1}.inlet_node(s);
    if sum(outlet_nodes_to_increase == curr_inlet_node) > 0
        additional_inlets_to_increase(s) = true;
    end
end
inlets_to_increase = segs_to_increase | additional_inlets_to_increase;
for sim_no = 1 : length(variations.ages)
    parameters.network_spec{sim_no}.inlet_radius(inlets_to_increase) = parameters.network_spec{sim_no}.inlet_radius(inlets_to_increase)*prop_of_baseline(sim_no);
    parameters.network_spec{sim_no}.outlet_radius(segs_to_increase) = parameters.network_spec{sim_no}.outlet_radius(segs_to_increase)*prop_of_baseline(sim_no);
end

%% -- lvet -- 
curr_param = 'lvet';

% extract variations for this parameter
curr_param_col = find(strcmp(variations.param_names, curr_param));
param_variations.no_sd = variations.params(:,curr_param_col);
param_variations.age_ind = variations.age_nos;

% setup model for lvet (based on data from Weissler 1961)
hr = [68,68,85,72,49,120,56,60,66,65,56,72,62,68,57,65,80,   93,84,80,74,100]'; sv = [82,81,62,95,98,44,106,102, 96,82,92,74,107,93,101,109,79,    49,71,59,76,43]'; lvet = 1000*[0.27,0.275,0.235,0.265,0.290,0.2,0.315,0.32,0.31,0.285,0.305,0.260,0.300,0.270,0.305,0.285,0.250,    0.2,0.235,0.24,0.24,0.19]';
tbl = table(hr,sv,lvet);
modelspec = 'lvet ~ 1 + hr + sv';
mdl = fitlm(tbl,modelspec);  % model of LVET as a function of HR and SV
[lvet_val,lvet_ci] = predict(mdl,[parameters.hr, parameters.sv], 'Alpha', 0.3173);  % I think a coefficient of 0.3173 gives one SD
lvet_sd = lvet_val - lvet_ci(:,1);  % SD of LVET

% calculate values of this parameter for each simulation
param_values = lvet_val + (param_variations.no_sd .* lvet_sd);

% normalise LVETs for each age to have the desired mean and SD
baseline_vals.mean = 282; % Using Mynard 2015's value for the mean
baseline_vals.sd = 282*(24.43/295.29); % Use the values from Gerstenblith 1977 (295.29 +/- 24.43 ms) to find the SD

ages = unique(variations.ages);
for age_no = 1 : length(ages)
    
    % identify data corresponding to this age
    rel_els = variations.ages == ages(age_no);
    rel_lvet = param_values(rel_els);
    
    % normalise lvet data
    rel_lvet = baseline_vals.mean * rel_lvet/mean(rel_lvet);
    if std(rel_lvet) > 0  % i.e. if there's more than one simulation at this age.
        rel_lvet = mean(rel_lvet) + (baseline_vals.sd * (rel_lvet-mean(rel_lvet))/std(rel_lvet));
    end
    
    % re-insert normalised data
    param_values(rel_els) = rel_lvet;
    
end

% store these values
eval(['parameters.' curr_param ' = param_values;'])

%% -- pvr --
% scale the peripheral vascular resistance Windkessel values to achieve the desired MBP
sims_co = parameters.hr.*parameters.sv/1000;        % cardiac output for each simulation in l/min
mean_flow = sims_co/(1000*60);                      % cardiac output converted to m3/sec
network_r = (133.33*parameters.p_drop)./mean_flow;  % resistance of network
parameters.pvr = (((133.33*parameters.mbp) - parameters.p_out)./mean_flow) - network_r; % where PVR is in Pa s/m3

%% -- pwv --
% Find the appropriate constants (k1, k2, k3) for the equation describing
% the relationship between arterial stiffness and arterial radius for each
% simulation. This is done in four steps.

% STEP 1: Find reference PWV values for different ages and different MBPs

% original data from Mattace-Raso2010, table 6: PWV = A*age + B*age^2 + C specified for different BP ranges 
A = [0.000, 0.000, 0.000, 0.000, 0.044];
B = [0.83, 0.99, 1.05, 1.18, 0.85]*1e-3;
C = [5.55, 5.69, 5.91, 6.17, 5.73];
sbp = [115, 125, 135, 150, 170];
dbp = [77.5, 82.5, 87.5, 95, 105];
init_mbp_vals = dbp+(0.4*(sbp-dbp));  % calculate approximate MBP corresponding to these SBP and DBP values using 0.4 constant from article; similar to 0.412 constant from: "Formula and nomogram for the sphygmomanometric calculation of the mean arterial pressure", http://heart.bmj.com/content/heartjnl/84/1/64.full.pdf

% extend range of MBP values by extrapolating at:
% - lower end
no_to_interpolate = 5;
interp_els = (1-no_to_interpolate):0;
new_sbp_vals = interp1(1:3, sbp(1:3), interp_els, 'linear', 'extrap');
new_dbp_vals = interp1(1:3, dbp(1:3), interp_els, 'linear', 'extrap');
new_A_vals = interp1(1:3, A(1:3), interp_els, 'linear', 'extrap');
new_B_vals = interp1(1:3, B(1:3), interp_els, 'linear', 'extrap');
new_C_vals = interp1(1:3, C(1:3), interp_els, 'linear', 'extrap');
% store these new values
sbp = [new_sbp_vals, sbp];
dbp = [new_dbp_vals, dbp];
A = [new_A_vals, A];
B = [new_B_vals, B];
C = [new_C_vals, C];

% - upper end
no_to_interpolate = 3;
interp_els = 2+(1:no_to_interpolate);
new_sbp_vals = interp1(1:2, sbp([end-1, end]), interp_els, 'linear', 'extrap');
new_dbp_vals = interp1(1:2, dbp([end-1, end]), interp_els, 'linear', 'extrap');
new_A_vals = interp1(1:2, A([end-1, end]), interp_els, 'linear', 'extrap');
new_B_vals = interp1(1:2, B([end-1, end]), interp_els, 'linear', 'extrap');
new_C_vals = interp1(1:2, C([end-1, end]), interp_els, 'linear', 'extrap');
% store these new values
sbp = [sbp, new_sbp_vals];
dbp = [dbp, new_dbp_vals];
A = [A, new_A_vals];
B = [B, new_B_vals];
C = [C, new_C_vals];

% calculate mbp
mbps.vals = dbp+(0.4*(sbp-dbp));  % approximation as above: 0.4 constant from article; 0.412 constant from: http://heart.bmj.com/content/heartjnl/84/1/64.full.pdf
mbps.inds = 1 : length(mbps.vals);

% data from Mattace-Raso2010 table 5 (median and 10 and 90 pc, therefore these represent +/- 40% in each direction)
no_sds = range(norminv([0.5, 0.9])); % in this case, for 80% confidence interval

sd_percentile.lower_v = (1/no_sds)*100* [(6.0-5.2)/6.0, (6.4-5.7)/6.4, (6.7-5.8)/6.7, (7.2-5.7)/7.2, (7.6-5.9)/7.6; ...
    (6.5-5.4)/6.5, (6.7-5.3)/6.7, (7.0-5.5)/7.0, (7.2-5.5)/7.2, (7.6-5.8)/7.6; ...
    (6.8-5.8)/6.8, (7.4-6.2)/7.4, (7.7-6.5)/7.7, (8.1-6.8)/8.1, (9.2-7.1)/9.2; ...
    (7.5-6.2)/7.5, (8.1-6.7)/8.1, (8.4-7.0)/8.4, (9.2-7.2)/9.2, (9.7-7.4)/9.7; ...
    (8.7-7.0)/8.7, (9.3-7.6)/9.3, (9.8-7.9)/9.8, (10.7-8.4)/10.7, (12.0-8.5)/12.0; ...
    (10.1-7.6)/10.1, (11.1-8.6)/11.1, (11.2-8.6)/11.2, (12.7-9.3)/12.7, (13.5-10.3)/13.5 ...
    ];
sd_percentile.upper_v = (1/no_sds)*100* [(7.0-6.0)/6.0, (7.5-6.4)/6.4, (7.9-6.7)/6.7, (9.3-7.2)/7.2, (9.9-7.6)/7.6; ...
    (7.9-6.5)/6.5, (8.2-6.7)/6.7, (8.8-7.0)/7.0, (9.3-7.2)/7.2, (11.2-7.6)/7.6; ...
    (8.5-6.8)/6.8, (9.0-7.4)/7.4, (9.5-7.7)/7.7, (10.8-8.1)/8.1, (13.2-9.2)/9.2; ...
    (9.2-7.5)/7.5, (10.4-8.1)/8.1, (11.3-8.4)/8.4, (12.5-9.2)/9.2, (14.9-9.7)/9.7; ...
    (11.4-8.7)/8.7, (12.2-9.3)/9.3, (13.2-9.8)/9.8, (14.1-10.7)/10.7, (16.5-12.0)/12.0; ...
    (13.8-10.1)/10.1, (15.5-11.1)/11.1, (15.8-11.2)/11.2, (16.7-12.7)/12.7, (18.2-13.5)/13.5 ...
    ];
sd_percentile.age = 25:10:75;

% add on extrapolated values (assuming a constant percentage value for the SD outside of the specified data range)
sd_percentile.mbp = [0, init_mbp_vals, 500];
sd_percentile.lower_v = [sd_percentile.lower_v(:,1), sd_percentile.lower_v, sd_percentile.lower_v(:,end)];
sd_percentile.upper_v = [sd_percentile.upper_v(:,1), sd_percentile.upper_v, sd_percentile.upper_v(:,end)];

% find reference PWV values for different ages and different MBPs
age = 20:80;
for n = 1 : length(A)
    pwv_age(:,n) = A(n)*age + B(n)*(age.^2) + C(n);
end

% STEP 2: Calculate mean and SD values for carotid-femoral PWV for each simulation according to age and MBP

% calculate expected MBP for each simulation
sims_co = parameters.hr.*parameters.sv/1000; % in l/min
mean_flow = sims_co/(1000*60); % in m3/sec
sims_mbp = ((mean_flow.*(parameters.pvr+network_r)) + parameters.p_out)/133.33;   % where PVR is in Pa s/m3

% setup variables
[initial_pwv.cf.v, initial_pwv.cf.sd_perc_lower, initial_pwv.cf.sd_perc_upper, parameters.desired_pwv_aorta, parameters.desired_pwv_leg, parameters.desired_pwv_arm, parameters.expected_pwv_aorta, parameters.expected_pwv_leg, parameters.expected_pwv_arm] = deal(nan(length(parameters.age),1));
for sim_no = 1 : length(parameters.age)
    
    % find values for this age
    curr_age = parameters.age(sim_no);
    [~, rel_age_el] = min(abs(age-curr_age));
    rel_pwv_age_vals = pwv_age(rel_age_el,:);
    
    % find values for this MBP from within those for this age
    rel_mbp_ind = interp1(mbps.vals, mbps.inds, sims_mbp(sim_no), 'linear', 'extrap');
    if rel_mbp_ind > length(mbps.vals)
        rel_mbp_ind = length(mbps.vals);
    end
    if rel_mbp_ind < 1
        rel_mbp_ind = 1;
    end
    
    % find mean value for cf PWV for this age and mbp
    initial_pwv.cf.v(sim_no) = interp1(mbps.inds, rel_pwv_age_vals, rel_mbp_ind);
    parameters.desired_pwv_leg(sim_no) = 0.056.*((initial_pwv.cf.v(sim_no)-6.15)/0.092)+7.91; % obtained from Avolio 1983 by combining eqns for PWVao and PWVleg
    parameters.desired_pwv_arm(sim_no) = 0.048.*((initial_pwv.cf.v(sim_no)-6.15)/0.092)+9.98; % obtained from Avolio 1983 by combining eqns for PWVao and PWVarm
    
    % find SD value for cf PWV for this age and mbp
    rel_sd_percentiles_lower = sd_percentile.lower_v(sd_percentile.age == curr_age, :);
    rel_sd_percentiles_upper = sd_percentile.upper_v(sd_percentile.age == curr_age, :);
    initial_pwv.cf.sd_perc_lower(sim_no,1) = interp1(sd_percentile.mbp, rel_sd_percentiles_lower, sims_mbp(sim_no));
    initial_pwv.cf.sd_perc_upper(sim_no,1) = interp1(sd_percentile.mbp, rel_sd_percentiles_upper, sims_mbp(sim_no));
    
end

% STEP 3: Calculate the desired cf-PWV value for each simulation based on:
% (i) the mean and SD values for this age and MBP, and
% (ii) the variation in PWV from the mean for this simulation.

% calculate actual cf-pwv values, taking into account variation from mean values for the age group
curr_param = 'pwv';
curr_param_col = find(strcmp(variations.param_names, curr_param));
param_variations.no_sd = variations.params(:,curr_param_col);
lower_els = param_variations.no_sd < 0;
parameters.desired_pwv_aorta = nan(size(lower_els));
parameters.desired_pwv_aorta(~lower_els) = initial_pwv.cf.v(~lower_els) .* (1 + (param_variations.no_sd(~lower_els).*(initial_pwv.cf.sd_perc_upper(~lower_els)/100)));
parameters.desired_pwv_aorta(lower_els) = initial_pwv.cf.v(lower_els) .* (1 + (param_variations.no_sd(lower_els).*(initial_pwv.cf.sd_perc_lower(lower_els)/100)));

% STEP 4: Find appropriate values of k1, k2 and k3 to achieve these desired cf-PWVs

k1 = 3.0e6;  % [g/s^2/cm]   % (Mynard's)
k2 = -9*1.5;  % [/cm]       % This is increased by a factor of 1.5 from Mynard 2015's value to achieve the desired pulse wave shapes (and amplification)
k3 = 3.37e5;  % [g/s^2/cm]   % (Mynard's)

% Prepare figure to plot wave speeds
if up.do_plot && strcmp(up.db_config, 'baseline_subjects_at_each_age')
    max_vals = [0,0]; min_val = inf;
    paper_size = [1200,400]; figure('Position', [20,20,paper_size])
end
for sim_no = 1 : length(parameters.age)

    % extract desired cf PWVs for this simulation
    desired_wave_speeds.aorta = parameters.desired_pwv_aorta(sim_no);
    desired_wave_speeds.arm = parameters.desired_pwv_arm(sim_no);
    desired_wave_speeds.leg = parameters.desired_pwv_leg(sim_no);

    % extract network parameters for this simulation
    rel_network_spec = parameters.network_spec{sim_no};

    % extract density of blood
    rho = parameters.rho(sim_no);

    % determine optimal value for k3 using a search algorithm (which
    % minimises the difference between the desired and expected cf PWV).
    options = optimset('MaxFunEvals',200, 'display', 'off');
    optimal_k1_and_2 = [k1,k2];
    k3_cost_function = @(k3)calculate_k3_cost_function(rel_network_spec, optimal_k1_and_2, k3, rho, desired_wave_speeds);
    [optimal_k3, ~] = fminsearch(k3_cost_function,k3, options);
    optimal_k = [k1, k2, optimal_k3];
    
    % store the optimal values of the constants
    store_k(sim_no,:) = optimal_k;
    parameters.network_spec{sim_no}.k = optimal_k;

    optimal_wave_speeds = calculate_theoretical_wave_speeds(rel_network_spec, optimal_k, rho);

    errors(sim_no,:) = [desired_wave_speeds.aorta - optimal_wave_speeds.aorta, ...
        desired_wave_speeds.arm - optimal_wave_speeds.arm, ...
        desired_wave_speeds.leg - optimal_wave_speeds.leg];

    parameters.expected_pwv_aorta(sim_no) = optimal_wave_speeds.aorta;
    parameters.expected_pwv_arm(sim_no) = optimal_wave_speeds.arm;
    parameters.expected_pwv_leg(sim_no) = optimal_wave_speeds.leg;

    if up.do_plot && strcmp(up.db_config, 'baseline_subjects_at_each_age')
        
        % Find expected wave speeds for all arterial segments
        ave_radius = mean([rel_network_spec.inlet_radius(:), rel_network_spec.outlet_radius(:)], 2);
        ave_radius_cm = 100*ave_radius;
        Eh_D0 = (optimal_k(1)*exp(optimal_k(2)*ave_radius_cm))+optimal_k(3); % Eh/D0 (from Mynard's 2015 paper, eqn 3)
        c0_squared = (2/3)*(Eh_D0/(rho/1000)); % from Mynard's 2015 paper, eqn 3, converts rho from kg/m3 to g/cm3
        wave_speed = sqrt(c0_squared)/100;  % converts from cm/s to m/s.
        [~, order] = sort(ave_radius_cm);
        
        % plot
        color_intensity = 0.2+0.5*(sim_no/length(parameters.age));
        curr_color = color_intensity*[1,1,1];
        leg_h(sim_no,1) = plot(2*10*ave_radius_cm(order), wave_speed(order), 'Color', curr_color, 'LineWidth', 2);
        hold on
        
        % tidy up
        ftsize = 28;
        set(gca, 'FontSize', ftsize)
        ylabel('Wave Speed [m/s]', 'FontSize', ftsize);
        xlabel('Arterial Diameter [mm]', 'FontSize', ftsize);
        title('PWV vs Diameter', 'FontSize', ftsize);
        labels{sim_no,1} = num2str(parameters.age(sim_no));
        box off
        max_vals(1) = max([max_vals(1), max(2*10*ave_radius_cm)]);
        min_val = min([min_val, min(wave_speed)]);
        max_vals(2) = max([max_vals(2), max(wave_speed)]);
    end

end

if up.do_plot && strcmp(up.db_config, 'baseline_subjects_at_each_age')
    xlim([0 1.15*max_vals(1)])
    ylim([0.95*min_val max_vals(2)+0.05*min_val])
    lgd = legend(flipud(leg_h), flipud(labels), 'FontSize', ftsize, 'Location', 'SouthEast');
    title(lgd,'Age')
    PrintFigs(gcf, paper_size/70, [up.savefolder, 'Wave_speeds'])
    
end

close all
% Make plot
if up.do_plot
    
    paper_size = [900, 270];
    figure('Position', [20,20,paper_size])
    ftsize = 16;
    lims = [5.5, 12.5];
    subplot(1,3,1)
    plot(parameters.desired_pwv_aorta, parameters.expected_pwv_aorta, 'bx', 'MarkerSize', 8), hold on, 
    plot(lims, lims, 'k')
    xlabel('Desired', 'FontSize', ftsize), ylabel('Prescribed', 'FontSize', ftsize)
    title('Carotid-femoral PWV', 'FontSize', ftsize)
    set(gca, 'FontSize', ftsize)
    xlim(lims), ylim(lims)
    
    subplot(1,3,2)
    plot(parameters.desired_pwv_arm, parameters.expected_pwv_arm, 'bx', 'MarkerSize', 8), hold on, 
    plot(lims, lims, 'k')
    xlabel('Desired', 'FontSize', ftsize), ylabel('Prescribed', 'FontSize', ftsize)
    title('Brachial-radial PWV', 'FontSize', ftsize)
    set(gca, 'FontSize', ftsize)
    xlim(lims), ylim(lims)
    
    subplot(1,3,3)
    plot(parameters.desired_pwv_leg, parameters.expected_pwv_leg, 'bx', 'MarkerSize', 8), hold on, 
    plot(lims, lims, 'k')
    xlabel('Desired', 'FontSize', ftsize), ylabel('Prescribed', 'FontSize', ftsize)
    title('Femoral-ankle PWV', 'FontSize', ftsize)
    set(gca, 'FontSize', ftsize)
    xlim(lims), ylim(lims)
    
    % Save figure
    PrintFigs(gcf, paper_size/70, [up.savefolder, 'Prescribed_vs_desired_PWVs'])
    
end

%% --- gamma coefficients ---
baseline_Gamma_b0 = 400*1.5; % a factor of 1.5 times larger than 400 g/s in Mynard 2015
baseline_Gamma_b1 = 100*1.5; % a factor of 1.5 times larger than 100 g cm/s in Maynard 2015
[parameters.gamma_b1, parameters.gamma_b0] = deal(nan(size(parameters.age)));

% find out how many SD to vary gamma by
curr_param = 'gamma';
curr_param_col = strcmp(variations.param_names, curr_param);
param_variations.no_sd = variations.params(:,curr_param_col);
parameters.gamma_b0 = baseline_Gamma_b0 .* (1 + param_variations.no_sd);
parameters.gamma_b1 = baseline_Gamma_b1.* (1 + param_variations.no_sd); 
    
%% Save table of equations
make_table_equations(up.savefolder, eqns);

end

function make_param_plots(up)

% setup
required_age = 25:1:75;
params = {'pwv_cf', 'pwv_fa', 'pwv_br', 't_pf', 'lvet', 'reg_vol', 'pvc', 'dbp', 'mbp', 'hr', 'sv', 'len', 'dia_asc', 'dia_desc_thor', 'dia_abd', 'dia_carotid'};
ftsize = 32;
lwidth = 2;
paper_size = [500,400];

% Extract data from literature for these params
for curr_age_no = 1 : length(required_age)
    curr_age = required_age(curr_age_no);
    
    for param_no = 1 : length(params)
        curr_param = params{param_no};
        
        if strcmp(curr_param, 'lvet')
            baseline_val = 295.29;
            sd = 24.43;
        elseif strcmp(curr_param, 'pwv_cf')
            raw_vals = [5.94618071428571;6.53537900000000;7.39549485714286;8.28360714285714;9.23771257142857;10.3004642857143];
            baseline_val = interp1(25:10:75, raw_vals, curr_age);
            raw_sd_vals = [10.4040552809651;13.2051470873787;11.8934986269314;13.5115098873697;15.0291697484146;19.1777415260153].*raw_vals/100;
            sd_min = interp1(25:10:75,raw_sd_vals, curr_age);
            raw_sd_vals = [13.0050691012063;16.8065508384820;18.5716246860543;19.2832834122860;24.2420005849521;28.7696946479979].*raw_vals/100;
            sd_max = interp1(25:10:75,raw_sd_vals, curr_age);
            clear raw_vals raw_sd_vals
        elseif strcmp(curr_param, 'pwv_fa')
            baseline_val = (5.6*curr_age+791)/100;
            sd_min = zeros(size(curr_age));
            sd_max = zeros(size(curr_age));
        elseif strcmp(curr_param, 'pwv_br')
            baseline_val = (4.8*curr_age+998)/100;
            sd_min = zeros(size(curr_age));
            sd_max = zeros(size(curr_age));
        else
            % calculate this parameter's baseline val and SD for this age
            [baseline_val, sd, eqn] = extract_values_from_literature(curr_age, curr_param, up);
            %[baseline_val, sd, eqn] = calculate_val(curr_age, curr_param, value_source, up);
        end
        
        % store results
        eval(['phys_vals.' curr_param '.val(curr_age_no,1) = baseline_val;'])
        if exist('sd', 'var')
            eval(['phys_vals.' curr_param '.sd_min(curr_age_no,1) = sd;'])
            eval(['phys_vals.' curr_param '.sd_max(curr_age_no,1) = sd;'])
        else
            eval(['phys_vals.' curr_param '.sd_min(curr_age_no,1) = sd_min;'])
            eval(['phys_vals.' curr_param '.sd_max(curr_age_no,1) = sd_max;'])
        end
        
        clear baseline_val sd curr_param abbr eqn
        
    end
    clear param_no curr_age
    
end
clear curr_age_no

% Make a plot of each of these params in turn
for param_no = 1 : length(params)
    curr_param = params{param_no};
    
    % get relevant data
    eval(['rel_data = phys_vals.' curr_param ';'])
    if strcmp(curr_param, 'dbp')
        rel_data.val = rel_data.val/133.33;
        rel_data.sd_min = rel_data.sd_min/133.33;
        rel_data.sd_max = rel_data.sd_max/133.33;
    end
    if strcmp(curr_param, 'pwv_fa') || strcmp(curr_param, 'pwv_br')
        rel_data = phys_vals.pwv_cf;
    end
    rel_data.sd_min_1 = rel_data.val - rel_data.sd_min;
    rel_data.sd_min_2 = rel_data.val - 2*rel_data.sd_min;
    rel_data.sd_plus_1 = rel_data.val + rel_data.sd_max;
    rel_data.sd_plus_2 = rel_data.val + 2*rel_data.sd_max;
    
    % change peripheral PWVs
    if strcmp(curr_param, 'pwv_fa')
        rel_data.val = 0.056.*((rel_data.val-6.15)/0.092)+7.91; % obtained from \cite{Avolio1983} by combining eqns for PWVao and PWVleg
        rel_data.sd_min_1 = 0.056.*((rel_data.sd_min_1-6.15)/0.092)+7.91; 
        rel_data.sd_min_2 = 0.056.*((rel_data.sd_min_2-6.15)/0.092)+7.91; 
        rel_data.sd_plus_1 = 0.056.*((rel_data.sd_plus_1-6.15)/0.092)+7.91; 
        rel_data.sd_plus_2 = 0.056.*((rel_data.sd_plus_2-6.15)/0.092)+7.91; 
    elseif strcmp(curr_param, 'pwv_br')
        rel_data.val = 0.048.*((rel_data.val-6.15)/0.092)+9.98; % obtained from \cite{Avolio1983} by combining eqns for PWVao and PWVarm
        rel_data.sd_min_1 = 0.048.*((rel_data.sd_min_1-6.15)/0.092)+9.98; 
        rel_data.sd_min_2 = 0.048.*((rel_data.sd_min_2-6.15)/0.092)+9.98; 
        rel_data.sd_plus_1 = 0.048.*((rel_data.sd_plus_1-6.15)/0.092)+9.98; 
        rel_data.sd_plus_2 = 0.048.*((rel_data.sd_plus_2-6.15)/0.092)+9.98; 
    end
    
    % Make figure
    figure('Position', [20, 20, paper_size])
    subplot('Position', [0.19,0.21,0.80,0.70])
    fprintf(['\n - Making plot for ' curr_param])
    
    % plot data
    plot(required_age, rel_data.val, 'k', 'LineWidth', lwidth),
    hold on
    plot(required_age, rel_data.sd_min_1, '--k', 'LineWidth', lwidth)
    plot(required_age, rel_data.sd_min_2, '--k', 'LineWidth', lwidth)
    plot(required_age, rel_data.sd_plus_1, '--k', 'LineWidth', lwidth)
    plot(required_age, rel_data.sd_plus_2, '--k', 'LineWidth', lwidth)
    
    % tidy up
    xlim([min(required_age), max(required_age)]);
    temp = [min(rel_data.sd_min_2), max(rel_data.sd_plus_2)];
    ylim([min(temp)-0.1*range(temp), max(temp)+0.1*range(temp)])
    set(gca, 'FontSize', ftsize-2)
    xlabel('Age [years]', 'FontSize', ftsize)
    [label, units, abbr, graph_title, graph_title_no_units] = make_param_label(curr_param);
    if strcmp(curr_param, 'len')
        units = 'mm';
        abbr = 'Length';
    end
    if length(curr_param)>3 & strcmp(curr_param(1:3), 'dia')
        units = 'mm';
        abbr = 'Diameter';
    end
    temp = strfind(abbr, '_');
    if ~isempty(temp)
        abbr = abbr(1:temp-1);
    end
    ylabel([abbr, ' [' units ']'], 'FontSize', ftsize)
    title(graph_title_no_units, 'FontSize', ftsize)
    box off
    grid on
    
    % save
    savepath = [up.savefolder, 'changes_w_age_', curr_param];
    PrintFigs(gcf, paper_size/60, savepath)
    
end

end

function make_haemod_plots(parameters, up)

% setup
params = fieldnames(parameters);
params = params(~strcmp(params,'variations') & ~strcmp(params,'fixed') & ~strcmp(params,'wk_params')...
    & ~strcmp(params,'desired_pwv_aorta') & ~strcmp(params,'desired_pwv_arm') & ~strcmp(params,'desired_pwv_leg')...
    & ~strcmp(params,'network_spec') & ~strcmp(params,'age'));
ftsize = 32;
lwidth = 2;
paper_size = [500,400];
age = parameters.age;
required_age = unique(age);

% Make a plot of each of these params in turn
for param_no = 1 : length(params)
    curr_param = params{param_no};
    
    [rel_data.val, sds] = deal(nan(length(required_age),1));
    for age_no = 1 : length(required_age)
        
        rel_els = age == required_age(age_no);
        eval(['curr_data = parameters.' curr_param '(rel_els);'])
        if strcmp(curr_param, 'dbp')
            curr_data = curr_data/133.33;  % convert from Pa to mmHg
        end
        rel_data.val(age_no) = mean(curr_data);
        sds(age_no) = std(curr_data);
    end
    rel_data.sd_min_1 = rel_data.val - sds;
    rel_data.sd_min_2 = rel_data.val - 2*sds;
    rel_data.sd_plus_1 = rel_data.val + sds;
    rel_data.sd_plus_2 = rel_data.val + 2*sds;
    
    % Make figure
    figure('Position', [20, 20, paper_size])
    subplot('Position', [0.18,0.21,0.80,0.70])
    fprintf(['\n - Making plot for ' curr_param])
    
    % plot data
    plot(required_age, rel_data.val, 'k', 'LineWidth', lwidth),
    hold on
    plot(required_age, rel_data.sd_min_1, '--k', 'LineWidth', lwidth)
    plot(required_age, rel_data.sd_min_2, '--k', 'LineWidth', lwidth)
    plot(required_age, rel_data.sd_plus_1, '--k', 'LineWidth', lwidth)
    plot(required_age, rel_data.sd_plus_2, '--k', 'LineWidth', lwidth)
    
    % tidy up
    xlim([min(required_age), max(required_age)]);
    temp = [min(rel_data.sd_min_2), max(rel_data.sd_plus_2)];
    if range(temp) == 0
        temp = [temp(1)-1, temp(1)+1];
    end
    ylim([min(temp)-0.1*range(temp), max(temp)+0.1*range(temp)])
    set(gca, 'FontSize', ftsize-2)
    xlabel('Age [years]', 'FontSize', ftsize)
    [label, units, abbr, graph_title, graph_title_no_units] = make_param_label(curr_param);
    temp = strfind(abbr, '_');
    if ~isempty(temp)
        abbr = abbr(1:temp-1);
    end
    ylabel([abbr, ' [' units ']'], 'FontSize', ftsize)
    title(graph_title_no_units, 'FontSize', ftsize)
    box off
    grid on
    
    % save
    savepath = [up.savefolder, 'expected_changes_w_age_', curr_param];
    PrintFigs(gcf, paper_size/60, savepath)
    
end

end

function PrintFigs(h, paper_size, savepath)
set(h,'PaperUnits','inches');
set(h,'PaperSize', [paper_size(1), paper_size(2)]);
set(h,'PaperPosition',[0 0 paper_size(1) paper_size(2)]);
set(gcf,'color','w');
print(h,'-depsc',savepath)
%print(h,'-dpdf',savepath)
close all;

% save 
fid = fopen([savepath, '.txt'], 'w');
p = mfilename('fullpath');
p = strrep(p, '\', '\\');
fprintf(fid, ['Figures generated by:\n\n ' p '.m \n\n on ' datestr(today)]);
fclose all;

end

function wave_speeds = calculate_theoretical_wave_speeds(rel_network_spec, k, rho)

% find current wave speeds of each segment (based on average radius)
ave_radius = mean([rel_network_spec.inlet_radius(:), rel_network_spec.outlet_radius(:)], 2);
ave_radius_cm = 100*ave_radius;
wave_speed = empirical_wave_speed(ave_radius_cm, k, rho);

% find current wave speeds along certain paths
% - aorta:carotid-femoral (noting that carotid and femoral measurements are taken half way along carotid and femoral arteries)
lens_aorta_carotid = [rel_network_spec.length([1,2]); rel_network_spec.length(15)/2];
carotid_radius_cm = 100* (rel_network_spec.inlet_radius(15) - (0.25*(rel_network_spec.inlet_radius(15) - rel_network_spec.outlet_radius(15))));
path_speeds_aorta_carotid = [wave_speed([1,2]); empirical_wave_speed(carotid_radius_cm, k, rho)];

lens_aorta_femoral = [rel_network_spec.length([1,2,14,18,27,28,35,37,39,41,42,44]); rel_network_spec.length(46)/2];
femoral_radius_cm = 100* (rel_network_spec.inlet_radius(46) - (0.25*(rel_network_spec.inlet_radius(46) - rel_network_spec.outlet_radius(46))));
path_speeds_aorta_femoral = [wave_speed([1,2,14,18,27,28,35,37,39,41,42,44]); empirical_wave_speed(femoral_radius_cm, k, rho)];

path_len = sum(lens_aorta_femoral)-sum(lens_aorta_carotid);
time_taken = sum(lens_aorta_femoral./path_speeds_aorta_femoral) - sum(lens_aorta_carotid./path_speeds_aorta_carotid);
wave_speeds.aorta = path_len/time_taken;

% - arm:brachial-radial (noting that brachial measurements are taken half way along brachial)
lens_aorta_radial = rel_network_spec.length([1,2,14,19,21,22]);
path_speeds_aorta_radial = wave_speed([1,2,14,19,21,22]);
lens_aorta_brachial = [rel_network_spec.length([1,2,14,19]); rel_network_spec.length(21)*0.75];
lens_aorta_brachial = [rel_network_spec.length([1,2,14,19]); rel_network_spec.length(21)/2];   %%%%% CHANGE
brachial_radius_cm = 100* (rel_network_spec.inlet_radius(21) - (0.375*(rel_network_spec.inlet_radius(21) - rel_network_spec.outlet_radius(21))));
brachial_radius_cm = 100* (rel_network_spec.inlet_radius(21) - (0.5*(rel_network_spec.inlet_radius(21) - rel_network_spec.outlet_radius(21))));   %%%%%% CHANGE
path_speeds_aorta_brachial = [wave_speed([1,2,14,19]); empirical_wave_speed(brachial_radius_cm, k, rho)];

path_len = sum(lens_aorta_radial)-sum(lens_aorta_brachial);
time_taken = sum(lens_aorta_radial./path_speeds_aorta_radial) - sum(lens_aorta_brachial./path_speeds_aorta_brachial);
wave_speeds.arm = path_len/time_taken;

% - leg:femoral-ankle
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

function wave_speed_cost_function = calculate_k3_cost_function(rel_network_spec, k, k3, rho, desired_wave_speeds)

wave_speeds = calculate_theoretical_wave_speeds(rel_network_spec, [k(1), k(2), k3], rho);

wave_speed_diffs = [wave_speeds.aorta - desired_wave_speeds.aorta, ...
    wave_speeds.arm - desired_wave_speeds.arm, ...
    wave_speeds.leg - desired_wave_speeds.leg];

% - cost function is solely based on aortic wave speed
temp = wave_speed_diffs(1);
wave_speed_cost_function = sqrt(mean(temp.^2));

end

function [baseline_val, sd, eqn] = extract_values_from_literature(required_age, param, up)

% extract data from the most reliable article(s) for this parameter
article_data = extract_data_from_articles(param, up);

% interpolate article data (using best fit) to provide values at each age of interest
abs_data = fit_article_data(article_data, param, up);

% Calculate baseline value and SD
eqn = abs_data;
if sum(strcmp(fieldnames(eqn), 'mean_f'))
    % for linear changes where you can simply use the linear regression equation
    baseline_val = eqn.mean_f.p1 * required_age  +  eqn.mean_f.p2;
    sd = eqn.sd_f.p1 * required_age  +  eqn.sd_f.p2;
else
    % for non-linear changes where you have to query the value at this particular age
    rel_el = eqn.age == required_age;
    baseline_val = eqn.mean(rel_el);
    sd = eqn.sd(rel_el);
end

end

function article_data = extract_data_from_articles(param, up)

switch param
    
    case 'sv'
        
        % All data obtained from male european data given in Table A2a and A2b of Poppe 2015
        article_data.val.age = 20:80;
        sv_u = 102.35 - 0.3856*article_data.val.age;
        sv_l = 42.94 - 0.1193*article_data.val.age;
        article_data.val.v = (sv_u + sv_l)/2; % approx median value
        article_data.sd.age = 20:80;
        article_data.sd.mean_val = article_data.val.v;
        article_data.sd.v= (sv_u-sv_l)/range(norminv([0.05,0.95])); % the 90% range specified by the upper and lower reference values, divided by the number of SDs from the mean (measuring both directions) which contain 90% of normally distributed values.
        
    case 'hr'
        
        % Mean data obtained from male data in Fig.1 of Yashin2006
        article_data.val.age = 32.5:5:92.5;
        temp = [80,95,97.5,96,94,93,90,84,76,66,57,51,51];
        article_data.val.v = 66 + (12*(temp-5)/(99-5));
        
        % SD data calculated from Table 1 of Petersen2017 (assumed to remain constant with age) 
        
        means = [67,69,70]; sds = [10,12,11]; n = [240,333,231];
        overall_n = sum(n);
        overall_mean = sum(means.*n)/overall_n;
        overall_sd_squared = 0;
        for dataset_no = 1 : length(n)
            curr_n = n(dataset_no);
            curr_mean = means(dataset_no);
            curr_sd = sds(dataset_no);
            curr_contribution_to_sd_squared = ((curr_sd^2) * (curr_n-1)) + ...
                (curr_n*curr_mean^2) + (curr_n*(overall_mean^2)) - ...
                (2*curr_n*curr_mean*overall_mean);
            overall_sd_squared = overall_sd_squared + (curr_contribution_to_sd_squared/(overall_n-1));
        end
        overall_SD = sqrt(overall_sd_squared);
        
        article_data.sd.age = [51, 59, 68];
        article_data.sd.mean_val = ones(1,3)*overall_mean;
        article_data.sd.v = ones(1,3)*overall_SD;  
        
    case 'dia_asc' % 'Dia - Aorta Asc'
        
        % Mean data obtained from Fig.4 of Hickson2010 (mm) - L1
        article_data.val.age = [24,34,45,57,63,73];
        article_data.val.v = 10 + (25*[35.5,38.5,41.75,42.0,44.0,44.5]./47);
        
        % SD data obtained from male data in Table 2 of Agmon2003 (assumed to remain constant with age and distance along the aorta)
        % - obtained 8.9 % variation by considering male data at all three relevant aortic sites
        article_data.sd.age = article_data.val.age; % these are from Hickson2010
        article_data.sd.mean_val = article_data.val.v; % these are from Hickson2010
        article_data.sd.v = (8.9/100)*article_data.val.v;
        
    case 'dia_desc_thor' % 'Dia - Aorta D Thor'
        
        % Mean data obtained from Fig.4 of Hickson2010 (mm) - L2
        article_data.val.age = [24,34,45,57,63,73];
        article_data.val.v = 10 + (25*[48,53.5,58.5,58,65,64.5]./108);
        
        % SD data obtained from male data in Table 2 of Agmon2003 (assumed to remain constant with age and distance along the aorta)
        % - obtained 8.9 % variation by considering male data at all three relevant aortic sites
        article_data.sd.age = article_data.val.age; % these are from Hickson2010
        article_data.sd.mean_val = article_data.val.v;  % these are from Hickson2010
        article_data.sd.v = (8.9/100)*article_data.val.v;
        
    case 'dia_abd' % 'Dia - Aorta Abd'
        
        % Mean data obtained from Fig. 4 of Hickson2010 (mm) - L4
        article_data.val.age = [24,34,45,57,63,73];
        article_data.val.v = 10 + (25*[29.5,35,38,37,42.5,44]./108);
        
        % SD data obtained from male data in Table 2 of Agmon2003 (assumed to remain constant with age and distance along the aorta)
        % - obtained 8.9 % variation by considering male data at all three relevant aortic sites
        article_data.sd.age = article_data.val.age; % these are from Hickson2010
        article_data.sd.mean_val = article_data.val.v;  % these are from Hickson2010
        article_data.sd.v = (8.9/100)*article_data.val.v;
        
    case 'dia_carotid' % 'Dia - Carotid'
        
        % All data taken from the male data in Table 1 of Hansen1995 (mm)
        article_data.val.age = [14,25,46,60,71];
        article_data.val.v = [6.9, 7.7, 8.3, 8.2, 9.0];
        article_data.sd.age = [14,25,46,60,71];
        article_data.sd.mean_val = [6.9, 7.7, 8.3, 8.2, 9.0];
        article_data.sd.v = [0.7, 0.4, 0.8, 0.8, 0.7];
        
    case 'len' % 'Length - Aorta Asc'
        
        % Mean data from Table 1 and Fig. 4 of Hickson2010 (mm) - ages in table 1, lengths in fig 4
        article_data.val.age = [24,34,45,57,63,73];  % From Hickson
        article_data.val.v = [90, 100, 100, 108, 118, 128];  % NEEDS MEASURING
        
        % SD data obtaine from Table 1 of Bensalah2014
        article_data.sd.age = [24,34,45,57,63,73];  % From Hickson
        article_data.sd.mean_val = [90, 100, 100, 108, 118, 128];  % NEEDS MEASURING
        article_data.sd.v = mean([14.0/105.5, 17.2/126.9])*article_data.sd.mean_val;  % from the two groups in Table1 of Bensalah 2014
        
    case 't_pf'
        
        % All data obtained from Table 1 of Bensalah2014 (ms)
        article_data.val.v = 104;
        article_data.sd.mean_val = 104;
        article_data.sd.v = 22;
        
        % All data obtained from PEF control group in Table 1 of Kamimura2016 (ms)
        article_data.val.v = 79;
        article_data.sd.mean_val = 79;
        article_data.sd.v = 11;
        
    case 'reg_vol'
        
        % All data obtained from Table 1 of Bensalah2014 (ml)
        article_data.val.v = 0.73;
        article_data.sd.mean_val = 0.73;
        article_data.sd.v = 0.63;
        
    case 'pvr'
                
        % Mean data obtained from Fig. 2 of McVeigh1999 (dyne s cm-5)
        article_data.val.age = 20:80;
        article_data.val.v = (8.1*article_data.val.age) + 926.9;  % invasive (Fig. 2)
        %article_data.val.v = (12.3*article_data.val.age) + 691.7;  % non-invasive (Fig. 3)
        article_data.sd.age = article_data.val.age;
        
        % SD data obtained from text in Bertand1985 (dyne s cm-5)
        article_data.sd.mean_val = article_data.val.v;
        article_data.sd.v = (180/1170)*article_data.sd.mean_val; % Data from Bertand1985 (dyne s cm-5)
        
        
    case 'pvc'
        
        % Mean data obtained from Fig.2 of McVeigh1999 (ml/mmHg)
        article_data.val.age = 20:80;
        article_data.val.v = (-0.001*article_data.val.age) + 0.113;  % invasive (Fig.2)
        %article_data.val.v = (-0.002*article_data.val.age) + 0.151;  % non-invasive (Fig. 3)
        article_data.sd.age = article_data.val.age;
        
        % SD data obtained from the results section (p.1245) of Resnick2000 (ml/mmHg)
        article_data.sd.mean_val = article_data.val.v;
        article_data.sd.v = (0.02/0.073)*article_data.sd.mean_val; % Data from Resnick2000, results p.1245 (ml/mmHg)
        
        % scale to make pvc a scaling factor
        scaling_factor = 1/article_data.val.v(article_data.val.age == 25);
        article_data.val.v = article_data.val.v*scaling_factor;
        article_data.sd.mean_val = article_data.val.v;
        article_data.sd.v = article_data.sd.v*scaling_factor;
        
    case 'p_out'
        
        % Data taken from DOI: 10.1152/jappl.1993.74.2.946 Parazynski 1993
        article_data.val.v = 33.2*133.33;
        article_data.sd.mean_val = article_data.val.v;
        article_data.sd.v = 1.5*sqrt(13)*133.3;
        
    case 'rho'
        
        article_data.val.v = 1060;  % kg/m3
        article_data.sd.mean_val = article_data.val.v;
        article_data.sd.v = 0.1*article_data.sd.mean_val; % no particular reason for using this value
        
    case 'mbp'
        
        % All data obtained from male data in Table 1 of McEniery2005
        article_data.val.age = 15.5:10:85.5;  % assumed mid-points of age categories
        article_data.val.v = [88, 89, 92, 95, 95, 94, 93, 92];        % i think these are obtained by integrating the central BP waveform, which was obtained from radial waveform using Sphygmocor
        article_data.sd.age = article_data.val.age;
        article_data.sd.mean_val = article_data.val.v;
        article_data.sd.v = [8, 8, 8, 7, 7, 7, 7, 8];
        
    case 'dbp'
        
        % All data obtained from male data in Table 1 of McEniery2005
        % although they are peripheral, looking at CSBP and CPP, they seem to be similar to central DBPs.
        article_data.val.age = 15.5:10:85.5;  % assumed mid-points of age categories
        article_data.val.v = 133.33*[73, 75, 77, 79, 79, 78, 76, 75];        % i think these are obtained by integrating the central BP waveform, which was obtained from radial waveform using Sphygmocor
        article_data.sd.age = article_data.val.age;
        article_data.sd.mean_val = article_data.val.v;
        article_data.sd.v = 133.33*[8, 8, 7, 6, 6, 6, 6, 9];
        
    case 'p_drop'
        
        article_data.val.v = 3.0;     
        article_data.val.v = 0;     
        article_data.sd.mean_val = article_data.val.v;
        article_data.sd.v = 10; % temporary data (no particular reason for using this value)
        
    case 'mu'
        
        article_data.val.v = 2.5e-3;  % Pa s (2e-3 provides a lower pressure drop between aorta and periphery than 3e-3).
        article_data.sd.mean_val = article_data.val.v;
        article_data.sd.v = 1e-3; % temporary data (no particular reason for using this value)
        
    case 'alpha'
        
        article_data.val.v = 4/3;
        article_data.sd.mean_val = article_data.val.v;
        article_data.sd.v = 0.16; % temporary data (no particular reason for using this value)
        
    case 'time_step'
        
        article_data.val.v = 5e-5;  % s
        article_data.sd.mean_val = article_data.val.v;
        article_data.sd.v = 2e-5; % temporary data (no particular reason for using this value)
        
    case 'visco_elastic_log'
        
        article_data.val.v = 1;  % logical: 1 indicates visco-elastic
        article_data.sd.mean_val = article_data.val.v;
        article_data.sd.v = 1;
        
    case 'reflect_log'
        
        article_data.val.v = 1;  % logical: 1 indicates normal reflections
        article_data.sd.mean_val = article_data.val.v;
        article_data.sd.v = 1;
end


end

function abs_data = fit_article_data(article_data, param, up)

change_type = identify_chosen_change(param);

abs_data.age = up.mod_ages;
thresh = 1e-10;

if strcmp(change_type, 'increase') || strcmp(change_type, 'decrease')
    
    % fit mean values (using linear fit)
    abs_data.mean_f=fit(article_data.val.age(:),article_data.val.v(:),'poly1');
    abs_data.mean = feval(abs_data.mean_f,abs_data.age); abs_data.mean = abs_data.mean(:)';
    % fit SDs (using linear fit)
    abs_data.sd_f=fit(article_data.sd.age(:),article_data.sd.v(:),'poly1');
    abs_data.sd = feval(abs_data.sd_f,abs_data.age); abs_data.sd = abs_data.sd(:)';
    
elseif strcmp(change_type, 'complex')
    
    % fit mean values (using non-linear fit)
    abs_data.mean = interp1(article_data.val.age, article_data.val.v, abs_data.age, 'pchip', 'extrap');
    % fit SDs (using linear fit)
    abs_data.sd_f=fit(article_data.sd.age(:),article_data.sd.v(:),'poly1');
    abs_data.sd = feval(abs_data.sd_f,abs_data.age); abs_data.sd = abs_data.sd(:)';
    
elseif strcmp(change_type, 'none')
    
    % copy across mean
    abs_data.mean = mean(article_data.val.v)*ones(length(abs_data.age),1);
    abs_data.mean_f.p1 = 0;
    abs_data.mean_f.p2 = mean(article_data.val.v);
    % copy across SD
    abs_data.sd = mean(article_data.sd.v)*ones(length(abs_data.age),1);
    abs_data.sd_f.p1 = 0;
    abs_data.sd_f.p2 = mean(article_data.sd.v);
        
end

end

function [change_type, sd_type] = identify_chosen_change(param)

switch param
    case 'sv'
        change_type = 'decrease';
    case 'hr'
        change_type = 'complex';
    case 'dia_asc'
        change_type = 'increase';
    case 'dia_desc_thor'
        change_type = 'increase';
    case 'dia_abd'
        change_type = 'increase';
    case 'dia_carotid'
        change_type = 'increase';
    case 'len'
        change_type = 'increase';
    case 'Length - Aorta D Thor'
        change_type = 'none';
    case 'Length - Aorta Abd'
        change_type = 'none';
    case 'Length - Carotid'
        change_type = 'none';
    case 'lvet'
        change_type = 'none';
    case 't_pf'
        change_type = 'none';
    case 'reg_vol'
        change_type = 'none';
    case 'pvr'
        change_type = 'increase';
    case 'pvc'
        change_type = 'decrease';
        sd_type = 'relative';
    case 'Stiffness - Radial'
        change_type = 'none';
    case 'Stiffness - Iliac'
        change_type = 'none';
    case 'Stiffness - Femoral'
        change_type = 'none';
    case 'Stiffness - Carotid'
        change_type = 'increase';
    case 'Stiffness - Brachial'
        change_type = 'increase';
    case 'Stiffness - Aorta D Thor'
        change_type = 'increase';
    case 'Stiffness - Aorta Asc'
        change_type = 'increase';
    case 'Stiffness - Aorta Abd'
        change_type = 'increase';
    case 'PWV - Aorta'
        change_type = 'increase';
    case 'PWV - Arm'
        change_type = 'increase';
    case 'PWV - Leg'
        change_type = 'increase';
    case 'p_out'
        change_type = 'none';
    case 'rho'
        change_type = 'none';
    case 'mu'
        change_type = 'none';
    case 'alpha'
        change_type = 'none';
    case 'time_step'
        change_type = 'none';
    case 'mbp'
        change_type = 'complex';
    case 'dbp'
        change_type = 'complex';
    case 'p_drop'
        change_type = 'none';
    case 'visco_elastic_log'
        change_type = 'none';
    case 'reflect_log'
        change_type = 'none';
end

end

function network_spec = obtain_network_spec(eqns, up)
% This function loads the network spec from the file specifying the
% arterial geometry.

% Load data
temp = tdfread(up.network_spec_file);

% Store relevant fields
network_spec.seg_no = temp.Segment_No;
network_spec.inlet_node = temp.Inlet_node;
network_spec.outlet_node = temp.Outlet_node;
network_spec.length = temp.Length_0x5Bm0x5D;
network_spec.inlet_radius = temp.Inlet_Radius_0x5Bm0x5D;
network_spec.outlet_radius = temp.Outlet_Radius_0x5Bm0x5D;
network_spec.segment_name = temp.Name;

% identify various properties
[network_spec.proximal_aorta, network_spec.aorta, network_spec.asc_aorta, network_spec.desc_thor_aorta, network_spec.abd_aorta, network_spec.both_carotid, network_spec.femoral_dorsalis, network_spec.brachial_radial] = deal(false(length(network_spec.seg_no),1));
network_spec.proximal_aorta([1,2]) = true;  % excluded 14 because it seemed to be beyond the usual definition of "aortic arch"
network_spec.aorta([1,2,14,18,27,28,35,37,39,41]) = true;
network_spec.asc_aorta([1,2]) = true;
network_spec.desc_thor_aorta([14,18,27]) = true; % needs checking
network_spec.abd_aorta([28,35,37,39,41]) = true; % needs checking
network_spec.both_carotid([5,15]) = true;
network_spec.femoral_dorsalis([46,49]) = true; % the dorsalis is "a continuation of the left anterior artery": https://en.wikipedia.org/wiki/Dorsalis_pedis_artery
network_spec.brachial_radial([21,22]) = true;  % left side
network_spec.asc_aorta_femoral([1,2,14,18,27,28,35,37,39,41,42,44]) = true;  % left side

end

function parameters = add_wk_vascular_bed_parameters(parameters)

fprintf('\n - Adding in Wk vascular bed parameters')

%% Load Windkessel data
parameters.wk_params = [6.9040,0.021726;20.755,0.026499;16.852,0.032637;16.852,0.032637;10.824,0.050811;8.368,0.065727;8.368,0.065727;21.252,0.02606;27.832,0.019761;54.084,0.010169;10.374,0.053018;26.119,0.021057;32.29,0.017033;32.29,0.017033;14.786,0.037197;14.824,0.037102];

end

function parameters = specify_fixed_simulation_parameters(parameters)

fprintf('\n - Specifying fixed simulation parameters')

%% ~~~~~~~~~~~~~~~~~~ General Model Settings ~~~~~~~~~~~~~~~~~~ 

%- Terminal conditions
%parameters.fixed.types.terminal = 'reflect';  % 'reflect' or 'absorb'

%- Arterial tree
parameters.fixed.arterial_tree = 'M116art'; % 'M55art', 'M96_art', or 'M116_art'

% specify which segments are to be modelled as visco-elastic (if visco-elastic modelling is being used):
parameters.fixed.visco_elastic_segments = 'some'; % 'all', 'hand', or 'some'
switch parameters.fixed.visco_elastic_segments
    % if element is less than or equal to this length then make elastic rather than visco-elastic
    case 'all'   % results in flow not being subsonic in simulation
        parameters.fixed.elmt_shortL_vw = 0.000; %Short element [m]
    case 'hand'  % ensures that all the arteries in the path of the digital artery are visco-elastic
        parameters.fixed.elmt_shortL_vw = 0.009; %Short element [m]
    case 'some'  % original setting from Marie's work  (NB: it doesn't seem to run if it's less than 0.011)
        parameters.fixed.elmt_shortL_vw = 0.011; %Short element [m]
end

%% ~~~~~~~~~~~~~~~~~~ Numerical Simulation Settings ~~~~~~~~~~~~~~~~~~ 

%- Simulation settings
parameters.fixed.duration = 6.5;     %Duration of simulation [s]  (originally 8.5 s)
parameters.fixed.output_dt = 1e-3;   %Time step of .his output history files [s]
parameters.fixed.iostep_t = 0.1;     %Time step of .out output files (flag -O) [s]
parameters.fixed.LinearForm = 0;     %1 if linear, 0 if nonlinear
parameters.fixed.p_order = 3;        %Polynomial order
parameters.fixed.q_order = 3;        %Quadrature order
parameters.fixed.elmt_shortL_reduced = 0.016; % reduces complexity if element is less than this length [m] % original setting from Marie's work
parameters.fixed.p_reduced = 2;      %Reduced order for short segments
parameters.fixed.q_reduced = 2;
parameters.fixed.int_time = 2;       %Time integration order (INTTYPE)
parameters.fixed.elmt = 0.02;        %Length of a single element in spatial discretisation [m]
parameters.fixed.R1_FIXED = 0;       %=1 if windkessel R1 is fixed via input dataFile; =0 if R1 = characteristic impedance Z0

%- Simulation outputs
parameters.fixed.artery_output = [1 2 5 7 8 14 15 18 21 22 27 38 39 41 46 49 55 56 60 61 62 63 68 71 74 75 76 85 87 102 103 104 105 112 113 114 115]; %Vessel number where output is created, at proximal, mid, and distal points
parameters.fixed.artery_output = [1 2 5 7 8 14 15 18 19 21 22 27 28 35 37 39 41 42 44 46 49 87 108 112]; % Adapted by PC to give only the segments he's interested in
parameters.fixed.artery_output = 1:116;  % trying all of them.

%% ~~~~~~~~~~~~~~~~~~ Post-processing Settings ~~~~~~~~~~~~~~~~~~ 

%- Arterial segments and paths of interest (used for post-processing)
parameters.fixed.aorta_path = [1 2 14 18 27 28 35 37 39 41];
% these go: segment number; relevant element of that segment; path of segments leading up to that segment
parameters.fixed.Art_aorta = 1;  parameters.fixed.Pt_aorta = 1; parameters.fixed.Path_aorta = [];
parameters.fixed.Art_asc   = 1;  parameters.fixed.Pt_asc = 2;   parameters.fixed.Path_asc = []; % +mid length 1
parameters.fixed.Art_desc  = 18; parameters.fixed.Pt_desc = 2;  parameters.fixed.Path_desc = [1 2 14]; % +mid length 18
parameters.fixed.Art_thoAo = 27; parameters.fixed.Pt_thoAo = 2; parameters.fixed.Path_thoAo = [1 2 14 18]; %+ mid length 27
parameters.fixed.Art_bif   = 41; parameters.fixed.Pt_bif = 3;   parameters.fixed.Path_bif = parameters.fixed.aorta_path;
parameters.fixed.Art_car   = 15; parameters.fixed.Pt_car = 2;   parameters.fixed.Path_car = [1 2]; % +mid length 15
parameters.fixed.Art_ren   = 38; parameters.fixed.Pt_ren = 2;   parameters.fixed.Path_ren = [1 2 14 18 27 28 35 37]; % +mid length 38
%parameters.fixed.Art_iliac = 50; parameters.fixed.Pt_iliac = 2; parameters.fixed.Path_iliac = [aorta_path, 42]; % +mid length 50
parameters.fixed.Art_brach = 21; parameters.fixed.Pt_brach = 3; parameters.fixed.Path_brach = [1 2 14 19 21]; % up to the end of brachial
parameters.fixed.Art_rad   = 22; parameters.fixed.Pt_rad = 2;   parameters.fixed.Path_rad = [1 2 14 19 21]; % +mid length 22
parameters.fixed.Art_dig   = 112;parameters.fixed.Pt_dig = 3;   parameters.fixed.Path_dig = [1 2 14 19 21 22 108 112]; % up to the end of digital
parameters.fixed.Art_ankle = 49; parameters.fixed.Pt_ankle = 3; parameters.fixed.Path_ankle = [parameters.fixed.aorta_path, 42, 44, 46, 49]; % up to the end of ankle
parameters.fixed.Art_fem   = 46; parameters.fixed.Pt_fem = 1;   parameters.fixed.Path_fem = [parameters.fixed.aorta_path, 42 44]; %
parameters.fixed.Path_brain = [1, 2, 15, 16, 79, 65, 96, 71, 72];
%parameters.fixed.Art_tibPost=54; parameters.fixed.Pt_tibPost = 3; parameters.fixed.Path_tibPost = [parameters.fixed.aorta_path, 42, 50, 52, 54];
parameters.fixed.arteriograph_path = [18 27 28 35 37 39 41]; % added arteriograph path length, from: DOI 10.1007/s10439-010-9945-1
parameters.fixed.relevant_paths = {'dig', 'ankle'};

end

function inflow = generate_aortic_inflow_waveforms(parameters)

fprintf('\n - Generating aortic inflow waves')

for sim_no = 1 : length(parameters.age)
    
    % setup input parameters for this subject
    input_params.wave_type = 'Mynard';   % use the template wave shape from Mynard 2015.
    input_params.HR = parameters.hr(sim_no);
    input_params.SV = parameters.sv(sim_no);
    input_params.T_Peak_Flow = parameters.t_pf(sim_no);
    input_params.Reg_Vol = parameters.reg_vol(sim_no);
    input_params.LVET = parameters.lvet(sim_no);
    input_params.plot_multiple = true;
    input_params.do_plot = false;
    
    % save space
    temp = AorticFlowWave_orig(input_params);
    temp = rmfield(temp, 't');
    temp.v = single(temp.v);
    inflow{sim_no} = temp; clear temp
    
    if rem(sim_no,100) == 0
        fprintf(['\n ' num2str(sim_no)])
    end
end

end

function make_table_equations(savefolder, eqns)

%% Find equation text

cardiac_params = {'HR', 'SV', 't_PF', 'reg_vol'};
arterial_params = {'MBP'};
vascular_params = {'PVC'};
params = [cardiac_params, arterial_params, vascular_params];

no_sf = 3;

for param_no = 1 : length(params)
    curr_param = params{param_no};
    eval(['curr_eqn = eqns.' lower(curr_param) ';'])
    
    % Formulate equation for mean
    if sum(strcmp(fieldnames(curr_eqn), 'mean_f'))
        curr_mean_eqn = num2str(round(curr_eqn.mean_f.p2, no_sf, 'significant'));
        if curr_eqn.mean_f.p1 ~= 0 && abs(curr_eqn.mean_f.p1) > 1e-10
            temp = num2str(round(curr_eqn.mean_f.p1, no_sf, 'significant'));
            if strcmp(temp(1), '-')
                temp = [' -- ', temp(2:end)];
            else
                temp = [' + ', temp];
            end
            curr_mean_eqn = [curr_mean_eqn, temp ' * age'];
        end
        
    else
        curr_mean_eqn = 'non-linear change';
    end
    
    % Formulate equation for SD
    curr_sd_eqn = num2str(round(curr_eqn.sd_f.p2, no_sf, 'significant'));
    if curr_eqn.sd_f.p1 ~= 0 && abs(curr_eqn.sd_f.p1) > 1e-10
        temp = num2str(round(curr_eqn.sd_f.p1, no_sf, 'significant'));
        if strcmp(temp(1), '-')
            temp = [' -- ', temp(2:end)];
        else
            temp = [' + ', temp];
        end
        curr_sd_eqn = [curr_sd_eqn, temp ' * age'];
    end
    
    
    % store equations
    eval(['eqns_text.' curr_param '.mean = curr_mean_eqn;'])
    eval(['eqns_text.' curr_param '.sd = curr_sd_eqn;'])
end


%% Create Table

table_text = ['\\textbf{Cardiac} & & \\\\', newline];
table_text = add_params_table_text(table_text, cardiac_params, eqns_text);
table_text = [table_text, '\\hline', newline];
table_text = [table_text, '\\textbf{Arterial} & & \\\\', newline];
table_text = add_params_table_text(table_text, arterial_params, eqns_text);
table_text = [table_text, '\\hline', newline];
table_text = [table_text, '\\textbf{Vascular Beds} & & \\\\', newline];
table_text = add_params_table_text(table_text, vascular_params, eqns_text);

fid = fopen([savefolder, 'eqns_table.txt'], 'w');
fprintf(fid, table_text);
fclose(fid);

end

function table_text = add_params_table_text(table_text, curr_params, eqns_text)

for param_no = 1 : length(curr_params)
    
    curr_param = curr_params{param_no};
    
    % Make label
    [label, units, abbr, graph_title] = make_param_label(curr_param);
    
    % Extract relevant text
    eval(['curr_eqns_text = eqns_text.' curr_param ';'])
    
    % Insert this line of the table
    table_line = ['- ', label];
    table_line = [table_line, ' & ' curr_eqns_text.mean ' & '  curr_eqns_text.sd ];
    table_line = [table_line, '\\\\', newline];
    table_text = [table_text, table_line];
    clear table_line
    
end

end

function eqns = calculate_equations(up)

cardiac_params = {'sv', 'hr', 't_pf', 'reg_vol'};
arterial_params = {'mbp', 'len'};
vascular_params = {'pvc', 'p_out'};
blood_params = {'rho'};
params = [cardiac_params, arterial_params, vascular_params, blood_params];

for param_no = 1 : length(params)
    
    param = params{param_no};
    
    %% Calculating values from raw data extracted from the literature
    
    % extract data from the most reliable article(s) for this parameter
    article_data = extract_data_from_articles(param, up);
    
    % interpolate article data (using best fit) to provide values at each age of interest
    eqn = fit_article_data(article_data, param, up);
    
    % Store equation for this parameter
    if sum(strcmp(fieldnames(eqn), 'mean_f'))
        eval(['eqns.' param '.mean_f = eqn.mean_f;']);
    end
    eval(['eqns.' param '.sd_f = eqn.sd_f;']);
    
end

% % Correct for length
% curr_init_length = eqns.len.mean_f.p1*up.baseline_age+eqns.len.mean_f.p2;
% scale_factor = 100/curr_init_length;
% eqns.len.mean_f.p1 = eqns.len.mean_f.p1*scale_factor;
% eqns.len.mean_f.p2 = eqns.len.mean_f.p2*scale_factor;
% eqns.len.sd_f.p1 = eqns.len.sd_f.p1*scale_factor;
% eqns.len.sd_f.p2 = eqns.len.sd_f.p2*scale_factor;

%% Make table for publication
make_table_equations(up.savefolder, eqns);

end