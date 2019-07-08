function create_pwdb_input_files
% CALCULATE_PWDB_INPUT_FILES creates a set of Nektar1D input files for each
% virtual subject in the Pulse Wave DataBase.
%
%               create_pwdb_input_files
%
%   Inputs:     a single file called "inputs.mat", which contains the
%                   parameters variable required by Nektar1D input files. This
%                   file can be created using the
%                   "calculate_pwdb_input_parameters.m" script.
%           
%
%   Outputs:    - input files for Nektar simulations
%               - shell script files to make it easier to run a batch of
%                   simulations
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
% v.1.0   Contributed to by: Marie Willemet, Jordi Alastruey, and Peter H. Charlton

fprintf('\n --- Creating Input Files ---')

%% Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTINGS TO CHANGE: This function specifies where to save the outputs   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
up = setup_up;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create input files
up = WriteInput_DB_pc(up);

%% Copy input files from temporary folder to Ubuntu shared folder
CopyInputFiles(up);

end

function up = setup_up

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This line specifies the location of the input parameters file  %%
% This should be the same as in calculate_pwdb_input_parameters.m  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

up.paths.savefolder = '/Users/petercharlton/Documents/Data/Nektar1D/ageing_sims/';  % (needs to have a slash at the end)
%up.paths.savefolder = '/home/pc13/Documents/Data/Nektar1D/ageing_sims/';
up.paths.input_parameters_filepath = [up.paths.savefolder, 'inputs.mat'];           % input parameters file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   These lines specify where to store temporary files created    %%
%                         to run the simulations                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (i) Computer running Matlab
up.paths.Matlab_computer_PathRoot         = '/Users/petercharlton/Desktop/temp/pwdb_files/';   % This is where the input files will be stored temporarily before being copied to the shared folder
up.paths.Matlab_computer_shared_folder    = '/Users/petercharlton/Documents/VM-share/';        % This is a shared folder which the input files are copied into (from the perspective of the Matlab computer)
%up.paths.Matlab_computer_PathRoot         = '/home/pc13/Documents/Data/Nektar1D/temp/pwdb_files/';   % This is where the input files will be stored temporarily before being copied to the shared folder
%up.paths.Matlab_computer_shared_folder    = '/home/pc13/Documents/Data/Nektar1D/VM-share/';        % This is a shared folder which the input files are copied into (from the perspective of the Matlab computer)
% (ii) Computer running Nektar1D (Ubuntu)
up.paths.Nektar_computer_shared_folder    = '/media/sf_VM-share/';                             % This is the shared folder (from the perspective of the computer used to run the simulations).
up.paths.Nektar_computer_sim_folder       = '/home/pc13/Documents/sim/';                       % This is the folder in which the output files from the simulations will be stored temporarily before being copied back to the shared folder.
%up.paths.Nektar_computer_shared_folder    = '/home/pc13/Documents/Data/Nektar1D/VM-share/';    % This is the shared folder (from the perspective of the computer used to run the simulations).
%up.paths.Nektar_computer_sim_folder       = '/media/pc13/6050B1A350B18078/Users/pc13/Simulation_data/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This line specifies how many command files to create for       %
%   running the simulations. It should be at least one, and (to    % 
%   run the simulations more quickly) the number of cores which    %
%               will be used to run the simulations.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

up.no_launch_files = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   The rest of this function can be left alone   %%%%%%%%

fprintf('\n - Setting up universal parameters')

%- file names
% used to run shell scripts to launch the simulations
up.paths.command_file       = 'command_list';
up.paths.launch_file        = 'launch_sims';
% used to record which simulations ran successfully
up.paths.FileNameResults_success = 'Nektar_Success';    % Where successful simulation names are listed during simulations
up.paths.FileNameResults_error = 'Nektar_Error';        % Where unsuccessful simulation names are listed during simulations

% - make folders if they don't already exist
if ~exist(up.paths.Matlab_computer_PathRoot)
    mkdir(up.paths.Matlab_computer_PathRoot)
end

end

function up = WriteInput_DB_pc(up)

%% Setup
fclose all;

% add directory of this file (and subdirectories) to path
filepath = fileparts(mfilename('fullpath'));
addpath(genpath(filepath)); clear filepath

% load virtual database parameters
load(up.paths.input_parameters_filepath);

% clear temporary folder
clear_temp_folder(up);

%% Create the input files for each simulation
fprintf('\n - Creating Input Files')

% Open command cloud file ready to write instructions (in temporary folder)
for file_no = 1 : up.no_launch_files
    command_file_names{file_no} = [up.paths.Matlab_computer_PathRoot, up.paths.command_file, num2str(file_no), '.txt'];
    fidCommand(file_no)  = fopen(command_file_names{file_no},'w');
end

for sim_no = 1 : length(parameters.age)
        
    fprintf(['\n Simulation: ' num2str(sim_no)])
    
    % Identify settings for this simulation
    curr_sim_setting = identify_sim_settings(parameters, inflow, sim_no);
    
    % Create Nektar input files (in temporary folder)
    FileName = create_input_files_for_a_simulation(curr_sim_setting, up);
    
    % Write commandLine into file
    current_command_file_no = rem(sim_no, up.no_launch_files);
    if current_command_file_no == 0
        current_command_file_no = up.no_launch_files;
    end
    curr_error_filename = [up.paths.FileNameResults_error, num2str(current_command_file_no)];
    curr_success_filename = [up.paths.FileNameResults_success, num2str(current_command_file_no)];
    %commandLine = ['oneDbio -L -R ', FileName, '.in && echo ''', FileName,''' >> ', curr_success_filename,'.txt || echo ''', FileName,''' >> ', curr_error_filename, '.txt']; clear curr_sim_setting
    % extra output files
    commandLine = ['oneDbio -O -L -d -t -R -s -m ', FileName, '.in && echo ''', FileName,''' >> ', curr_success_filename,'.txt || echo ''', FileName,''' >> ', curr_error_filename, '.txt']; clear curr_sim_setting
    fprintf(fidCommand(current_command_file_no), [commandLine,'\n']);
    
    clear FileName curr_sim_setting commandLine
    
end

% Finish command cloud file
for file_no = 1 : length(fidCommand)
    fclose(fidCommand(file_no));
end
clear fidCommand
fclose all;

fprintf(['\n All input files copied to ',strrep(up.paths.Matlab_computer_PathRoot, '\', '\\')])

%% Create launch file
WriteLaunchFile(up, command_file_names)

up.vdb_up = parameters.fixed;

end

function clear_temp_folder(up)

fprintf('\n - Clearing temporary folder ')

% delete any files in the temporary VDB folder
files = dir(up.paths.Matlab_computer_shared_folder);
if ~isempty(files)
    file_names = extractfield(files, 'name');
    file_dirs = cell2mat(extractfield(files, 'isdir')); clear files
    files_to_delete = file_names(~file_dirs); clear file_dirs file_names
    for file_no = 1 : length(files_to_delete)
        curr_filename = [up.paths.Matlab_computer_shared_folder, files_to_delete{file_no}];
        delete(curr_filename);
    end
    clear file_no files_to_delete
end

% Make folders to store files in
folder_path = up.paths.Matlab_computer_shared_folder;
if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

end

function sim_settings = identify_sim_settings(parameters, inflow, sim_no)

% copy across fixed parameters
sim_settings.fixed = parameters.fixed;

% copy across network spec for this simulation
sim_settings.network_spec = parameters.network_spec{sim_no};

% copy across wk params
sim_settings.wk_params = parameters.wk_params;

% copy across variations information
sim_settings.variations.params = parameters.variations.params(sim_no,:);
sim_settings.variations.param_names = parameters.variations.param_names;

% copy across remaining parameters
parameters = rmfield(parameters, {'fixed', 'network_spec', 'variations', 'wk_params'});

model_params.names = fieldnames(parameters);
for param_no = 1 : length(model_params.names)
    curr_param = model_params.names{param_no};
    eval(['sim_settings.' curr_param ' = parameters.' curr_param '(sim_no);']);
end

% copy across inflow waveform
sim_settings.inflow = inflow{sim_no}; 

sim_settings.filename = ['sim_' num2str(sim_no)];

end

function WriteLaunchFile(up, command_file_names)

fprintf('\n - Writing Launch Files')

if nargin == 1
    launch_file_names{1} = up.paths.launch_file;
    filepath{1} = [up.paths.Matlab_computer_PathRoot, launch_file_names{file_no}];
else
    for file_no = 1 : length(command_file_names)
        launch_file_names{file_no} = [up.paths.launch_file, num2str(file_no)];
        filepath{file_no} = [up.paths.Matlab_computer_PathRoot, launch_file_names{file_no}];
    end
end
newline_char = '\n';

% create starter file
starterfile_path = [up.paths.Matlab_computer_PathRoot, 'launch_starter'];
fid = fopen(starterfile_path, 'wt');
file_text = ['mkdir -p ' up.paths.Nektar_computer_sim_folder, newline_char, ...
            'cd ' up.paths.Nektar_computer_sim_folder, newline_char, ...
            'rm *.*', newline_char, ...
            'cp ' up.paths.Nektar_computer_shared_folder '*.* .', newline_char, ...
            'echo Finished operations'];
fprintf(fid, file_text);
fclose(fid);

% create finisher file
finisherfile_path = [up.paths.Matlab_computer_PathRoot, 'launch_finisher'];
fid = fopen(finisherfile_path, 'wt');
file_text = ['cd ' up.paths.Nektar_computer_sim_folder, newline_char, ...
            'cp *.* /media/sf_VM-share/.', newline_char, ...
            'echo Finished operations'];
fprintf(fid, file_text);
fclose(fid);

% create individual run files

for file_no = 1 : length(launch_file_names)
    
    fid = fopen(filepath{file_no}, 'wt');
    
    file_text = ['echo Launching file ', num2str(file_no), newline_char, ...
            'cd ', up.paths.Nektar_computer_sim_folder, newline_char, ...
            ['chmod u+x ' up.paths.command_file, num2str(file_no), '.txt'], newline_char, ...
            ['./' up.paths.command_file, num2str(file_no) '.txt'], newline_char, ...
            'echo Finished operations'];
    
    fprintf(fid, file_text);
    
    fclose(fid);
    
end

end

function CopyInputFiles(up)

fprintf('\n - Copying Input Files to Virtual Box Shared Folder')

% Identify input files to be copied
current_folder = up.paths.Matlab_computer_PathRoot;
input_files = dir(current_folder);
exc = cell2mat(extractfield(input_files, 'isdir'));
input_files = extractfield(input_files, 'name');
input_files = input_files(~exc); clear exc

% Copy each input file
for input_file_no = 1 : length(input_files)
    curr_file_name = input_files{input_file_no};
    copyfile([current_folder, curr_file_name], [up.paths.Matlab_computer_shared_folder, curr_file_name]);
end

end