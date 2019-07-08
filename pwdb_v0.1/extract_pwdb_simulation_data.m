function extract_pwdb_simulation_data(pwdb_no)
% EXTRACT_PWDB_SIMULATION_DATA extracts simulation data from the Pulse Wave
% DataBase simulation input and output files, and saves them in Matlab
% format.
%
%               extract_pwdb_simulation_data
%
%   Inputs:     the input and output files for the simulations, which
%                   should be stored in the shared folder.
%
%   Outputs:    - the files are copied into a new folder for storage
%               - a single Matlab file called 'collated_data' is generated
%                   containing the simulation data for the database.
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

fprintf('\n --- Extracting PWDB simulation data ---')

%% Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTINGS TO CHANGE: This function specifies where to save the outputs   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PATHS = setup_paths_for_post_processing;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copy input files and output files to a new folder for storage
if nargin==0
    pwdb_no = CopyFiles(PATHS);
end

% Setup paths with current simulation paths
PATHS = setup_paths_for_post_processing(pwdb_no);

% Convert output History files into Matlab format
up.all_data = false; up.all_beats = false;
up.dir = PATHS.OutputFiles; up.save_dir = PATHS.ProcessedData;
up.ds_factor = 2;
up.find_pw_els = 1;
if ~exist(PATHS.history_files_data, 'file')
    PATHS.history_files_data = ConvertHistoryFiles(up);
end

% Import arterial system characteristics
if ~exist(PATHS.system_chars, 'file')
   import_system_chars(PATHS);
end

% Import peripheral boundary conditions
if ~exist(PATHS.peripheral_boundarys, 'file')
   import_peripheral_boundarys(PATHS);
end

% Collate input and output data into single variable, and calculate PPGs
if ~exist(PATHS.collated_data, 'file')
    collate_input_and_output_data(PATHS);
end

fprintf('\n --- Finished extracting PWDB simulation data ---\n')

end

function pwdb_no = CopyFiles(PATHS)

% Identify the number of this pwdb
prev_pwdbs = dir(PATHS.storage_folder);
prev_pwdbs = extractfield(prev_pwdbs, 'name');
prev_pwdbs = prev_pwdbs(~cellfun(@isempty, strfind(prev_pwdbs, 'pwdb')));
if isempty(prev_pwdbs)
    pwdb_no = 1;
else
    pwdb_nos = nan(length(prev_pwdbs),1);
    for s = 1 : length(prev_pwdbs)
        pwdb_nos(s) = cell2mat(textscan(prev_pwdbs{s}, 'pwdb_%n'));
    end
    clear s
    pwdb_no = max(pwdb_nos)+1; clear pwdb_nos
end
clear prev_sims

% Check to see whether there are any files to be copied across:
temp = dir(PATHS.shared_folder);
possible_files = extractfield(temp, 'name');
possible_files = possible_files(~cell2mat(extractfield(temp, 'isdir')));
possible_files = possible_files(~strcmp(possible_files, '.DS_Store'));
if isempty(possible_files)
    pwdb_no = pwdb_no - 1;  % take most recent simulation
    return
end

fprintf('\n --- Copying Files ---')

% Identify files to be copied
files_to_copy = dir(PATHS.shared_folder);
exc = cell2mat(extractfield(files_to_copy, 'isdir'));
files_to_copy = extractfield(files_to_copy, 'name');
files_to_copy =files_to_copy(~exc); clear exc
files_to_copy =files_to_copy(~strcmp(files_to_copy, '.DS_Store'));

% Determine whether each file is an input or an output
file_type = nan(size(files_to_copy));
for file_no = 1 : length(files_to_copy)
    if ~isempty(strfind(files_to_copy{file_no}, 'command')) || ...
            ~isempty(strfind(files_to_copy{file_no}, 'launch')) || ...
            ~isempty(strfind(files_to_copy{file_no}, '.bcs')) || ...
            ~isempty(strfind(files_to_copy{file_no}, '.in')) || ...
            ~isempty(strfind(files_to_copy{file_no}, '.mat'))
        file_type(file_no) = 1;
    else
        file_type(file_no) = 2;
    end
end

% Copy each input file
for file_no = 1 : length(files_to_copy)
    curr_file_path = [PATHS.shared_folder, files_to_copy{file_no}];
    if file_type(file_no) == 1
        new_file_path = [input_folder, files_to_copy{file_no}];
    else
        new_file_path = [output_folder, files_to_copy{file_no}];        
    end
    if ~exist(new_file_path, 'file')
        copyfile(curr_file_path, new_file_path);
    end
    delete(curr_file_path)
end

end

function import_system_chars(PATHS)

% Check to see whether the required "period.tex" files have been outputted
% from the simulations:

rel_files = dir([PATHS.OutputFiles, '*_period.tex']);
if isempty(rel_files)
    return
end

fprintf('\n Importing arterial system characteristics')

% load output data
load(PATHS.history_files_data);
history_files_data = data; clear data
sims = fieldnames(history_files_data);

% setup input variables
[pvr, pvc, pvc_iw, ac, c, tau] = deal(nan(length(sims),1));

% extract data for each simulation
for sim_no = 1 : length(sims)
    % identify file for this simulation
    curr_sim = sims{sim_no};
    filepath = [PATHS.OutputFiles, curr_sim '_period.tex'];
    fid = fopen(filepath);
    % extract this file's data
    finished = false;
    while ~finished
        line_text = fgetl(fid);
        if ~isempty(strfind(line_text, ' Peripheral vascular resistance:'))
            % Extract PVR
            temp1 = strfind(line_text, ':');
            temp2 = strfind(line_text, 'Pa');
            rel_text = line_text(temp1+4:temp2-2);
            pvr(sim_no) = str2double(rel_text);
            
            % Extract PVC (on next line)
            line_text = fgetl(fid);
            temp1 = strfind(line_text, ':');
            temp2 = strfind(line_text, 'm$');
            rel_text = line_text(temp1+4:temp2-2);
            pvc(sim_no) = str2double(rel_text);
            
            % Extract impedance-weighted PVC (on next line)
            line_text = fgetl(fid);
            temp1 = strfind(line_text, ':');
            temp2 = strfind(line_text, 'm$');
            rel_text = line_text(temp1+5:temp2-2);
            pvc_iw(sim_no) = str2double(rel_text);
            
            % Extract total arterial compliance (on next line)
            line_text = fgetl(fid);
            temp1 = strfind(line_text, ':');
            temp2 = strfind(line_text, 'm$');
            rel_text = line_text(temp1+5:temp2-2);
            ac(sim_no) = str2double(rel_text);
            
            % Extract total compliance (on next line)
            line_text = fgetl(fid);
            temp1 = strfind(line_text, ':');
            temp2 = strfind(line_text, 'm$');
            rel_text = line_text(temp1+14:temp2-2);
            c(sim_no) = str2double(rel_text);
            
            % Extract time constant (on next line)
            line_text = fgetl(fid);
            temp1 = strfind(line_text, ':');
            temp2 = strfind(line_text, 's \\');
            rel_text = line_text(temp1+11:temp2-2);
            tau(sim_no) = str2double(rel_text);
            
            finished = true;
        end
    end
end

% save arterial system characteristics data
system_chars = table(pvr, pvc, pvc_iw, ac, c, tau);
save(PATHS.system_chars, 'system_chars')

end

function import_peripheral_boundarys(PATHS)

fprintf('\n Importing peripheral Windkessel boundary conditions')

% load output data
load(PATHS.history_files_data);
history_files_data = data; clear data
sims = fieldnames(history_files_data);

% extract data for each simulation
[c, r_sum] = deal(nan(length(sims),116));
for sim_no = 1 : length(sims)
    
    % identify file for this simulation
    curr_sim = sims{sim_no};
    filepath = [PATHS.InputFiles, curr_sim '.mat'];
    input_data = load(filepath);
    
    % Extract required Windkessel data
    if sum(strcmp(fieldnames(input_data.sim_settings), 'reflect_log')) && input_data.sim_settings.reflect_log && ~input_data.sim_settings.fixed.R1_FIXED
        temp.c_vals = input_data.params.outlet.C_WK;            %WK Compliance
        temp.c_vals(temp.c_vals == 0) = nan;
        c(sim_no,1:length(temp.c_vals)) = temp.c_vals;
        temp.r_sum_vals = input_data.params.outlet.R_WK_tot;    %WK Resistance (sum of R1 and R2)
        temp.r_sum_vals(temp.r_sum_vals == 0) = nan;
        r_sum(sim_no,1:length(temp.r_sum_vals)) = temp.r_sum_vals;
        clear temp
    else
        warning('\nHaven''t written this part yet')
        temp.c_vals = input_data.params.outlet.C_WK;            %WK Compliance
        temp.c_vals(temp.c_vals == 0) = nan;
        c(sim_no,1:length(temp.c_vals)) = temp.c_vals;
        temp.r_sum_vals = input_data.params.outlet.R_WK_tot;    %WK Resistance (sum of R1 and R2)
        temp.r_sum_vals(temp.r_sum_vals == 0) = nan;
        r_sum(sim_no,1:length(temp.r_sum_vals)) = temp.r_sum_vals;
        clear temp
    end
    
end

% save peripheral Windkessel characteristics data
peripheral_chars.c = c;
peripheral_chars.r_sum = r_sum;
save(PATHS.peripheral_boundarys, 'peripheral_chars')

end

function collate_input_and_output_data(PATHS)

fprintf('\n Collating input and output data')

% load output data
load(PATHS.history_files_data);
history_files_data = data; clear data

% Load system chars (if available)
if exist(PATHS.system_chars, 'file')
    load(PATHS.system_chars);
    do_system_chars = 1;
else
    do_system_chars = 0;
end

% import input, output and system chars data for each simulation
sims = fieldnames(history_files_data);
for sim_no = 1 : length(sims)
    % import input data
    curr_sim_input_data = [PATHS.InputFiles, sims{sim_no}];
    input_data = load(curr_sim_input_data); clear curr_sim_input_data
    input_data = rmfield(input_data, 'sim_up');
    collated_data(sim_no).input_data = input_data; clear input_data
    % import output data
    eval(['output_data = history_files_data.' sims{sim_no} ';'])
    collated_data(sim_no).output_data = output_data; clear output_data
    % import system chars
    if do_system_chars
        temp_fields = system_chars.Properties.VariableNames;
        for field_no = 1 : length(temp_fields)
            eval(['collated_data(sim_no).system_chars.' temp_fields{field_no} ' = system_chars.' temp_fields{field_no} '(sim_no);']);
        end
        clear temp_fields
    end
end

% Calculate WK PPG data
for sim_no = 1 : length(sims)
    
    for domain_el = 1 : length(collated_data(sim_no).output_data)
        
        no_distance_els = length(collated_data(sim_no).output_data(domain_el).distances);
        for distance_el = 1 : no_distance_els
            
            % Extract required data
            Q.v = collated_data(sim_no).output_data(domain_el).A(:,distance_el).*collated_data(sim_no).output_data(domain_el).U(:,distance_el);
            Q.fs = collated_data(sim_no).output_data(domain_el).fs;
            P.v = collated_data(sim_no).output_data(domain_el).P(:,distance_el);
            P.fs = collated_data(sim_no).output_data(domain_el).fs;
            P_out.v = collated_data(sim_no).input_data.sim_settings.p_out*ones(size(P.v));
            P_out.fs = collated_data(sim_no).output_data(domain_el).fs;
            Q1D.fs = collated_data(sim_no).output_data(domain_el).fs;
            Q1D.v = collated_data(sim_no).output_data(domain_el).Q1D;
            Q_out.fs = collated_data(sim_no).output_data(domain_el).fs;
            Q_out.v = collated_data(sim_no).output_data(domain_el).Q_out;
            
            % Estimate PPG using Windkessel approach
            curr_domain_no = collated_data(sim_no).output_data(domain_el).domain_no;
            curr_length = collated_data(sim_no).input_data.sim_settings.network_spec.length(curr_domain_no);
            end_of_segment_log = collated_data(sim_no).output_data(domain_el).distances(distance_el) == curr_length;
            clear curr_length curr_domain_no
            if ~isempty(Q1D.v) && end_of_segment_log
                curr_PPG_WK.v = cumsum(Q1D.v-Q_out.v);
%                 curr_PPG_WK.v = (temp-min(temp))./range(temp); clear temp
                curr_PPG_WK.v = calc_static_wave(curr_PPG_WK.v);
            else
                curr_PPG_WK = estimate_ppg_using_windkessel(P, Q, P_out);
            end
            curr_PPG_WK.v = curr_PPG_WK.v(:);
            [~,rel_el] = min(curr_PPG_WK.v);
            curr_PPG_WK.v = [curr_PPG_WK.v(rel_el:end); curr_PPG_WK.v(1:rel_el-1)];
            
            % find normalisation scaling factor
            norm_factor = 1./range(curr_PPG_WK.v);
            
            % normalise
            curr_PPG_WK.v = (curr_PPG_WK.v-min(curr_PPG_WK.v)).*norm_factor;
            
%             curr_PPG_WK.v = sqrt(curr_PPG_WK.v);  % Added to transform volume PPG into distance PPG
            collated_data(sim_no).output_data(domain_el).PPG(:,distance_el) = curr_PPG_WK.v;
            collated_data(sim_no).output_data(domain_el).PPG_start_sample(1,distance_el) = collated_data(sim_no).output_data(domain_el).start_sample(1,distance_el) + rel_el -1;
            collated_data(sim_no).output_data(domain_el).PPG_scaling_factor(1,distance_el) = norm_factor;
            clear rel_el curr_PPG_WK end_of_segment_log curr_length curr_domain_no norm_factor
            
        end
        collated_data(sim_no).output_data(domain_el).units.PPG = 'au';
        
    end
    
end

% Calculate forward and backward waves
do_extras = 0;
if do_extras
    for sim_no = 1 : length(sims)
        
        for domain_el = 1:length(collated_data(sim_no).output_data)
            
            no_distance_els = length(collated_data(sim_no).output_data(domain_el).distances);
            for distance_el = 1 : no_distance_els
                
                % Extract required data
                U.v = collated_data(sim_no).output_data(domain_el).U(:,distance_el);
                U.fs = collated_data(sim_no).output_data(domain_el).fs;
                P.v = collated_data(sim_no).output_data(domain_el).P(:,distance_el);
                P.fs = collated_data(sim_no).output_data(domain_el).fs;
                D.v = 2*sqrt(collated_data(sim_no).output_data(domain_el).A(:,distance_el)/pi);
                D.fs = collated_data(sim_no).output_data(domain_el).fs;
                PPG.v = collated_data(sim_no).output_data(domain_el).PPG(:,distance_el);
                
                % Find reservoir pressure
                Q.v = collated_data(sim_no).output_data(domain_el).Q(:,distance_el);
                R = mean(P.v)/mean(Q.v);
                
                C = 1e-8;
                tdia = collated_data(sim_no).input_data.sim_settings.lvet/1000; % ms to secs
                tdia = tdia+0.06;
                pout = collated_data(sim_no).input_data.sim_settings.p_out;
                pout = 0;
                options = optimset('MaxFunEvals',1000, 'display', 'off');
                input_vars = [C,pout];
                pres_cost_function = @(input_vars)calculate_pres_cost_function(P, Q, pout, R, tdia, input_vars);
                [optimal_input_vars, cost_func] = fminsearch(pres_cost_function,input_vars, options);
                %optimal_input_vars(2)
                Pres = calc_pres(P,Q,optimal_input_vars(2),R, optimal_input_vars(1));
%                 pres_cost_function = @(C)calculate_pres_cost_function(P, Q, pout, R, tdia, C);
%                 [optimal_C, cost_func] = fminsearch(pres_cost_function,C, options);
%                 Pres = calc_pres(P,Q,pout,R, optimal_C);
                Pexc = P.v - Pres.v;
                collated_data(sim_no).output_data(domain_el).Pexc(:,distance_el) = Pexc;
                
%                 rel_els = 309:length(P.v);
%                 rel_t = [rel_els]/P.fs; rel_t = rel_t(:);
%                 rel_p = P.v(rel_els);
%                 
%                 P_inf = collated_data(sim_no).input_data.sim_settings.p_out;
%                 %g = fittype('a*exp(b*x)');
%                 %g = fittype([num2str(P_inf), ' + ' num2str(P.v(1) - P_inf) '*exp(b*x)']);
%                 g = fittype(['a + ' num2str(P.v(1) - P_inf) '*exp(b*x)']);
%                 %options = fitoptions('Method', 'NonlinearLeastSquares','StartPoint', -0.5);
%                 myfit = fit(rel_t(:),rel_p,g); %, options);
%                 plot(myfit,rel_t,rel_p)
%                 RC = -1/myfit.b;
%                 %RC = 1.57*RC;
%                 %R = 1.57*R;
%                 C = RC/R;
%                 %P_inf = P.v(1)-myfit.a;
%                 %P_inf = collated_data(sim_no).input_data.sim_settings.p_out;
%                 P_0 = P.v(1);
%                 % New Method
%                 t = [0:length(P.v)-1]/P.fs; t = t(:);
%                 dt=1/P.fs;
%                 %RC = RC*1.57;
%                 qse=cumtrapz(Q.v.*exp(t/RC))*dt;
%                 P_wk = (exp(-t/RC).*(qse/C)) + (exp(-t/RC)*(P.v(1) - P_inf)) + P_inf;
%                 P_wk = P_wk/133.33;
%                 plot(P_wk, 'r')
%                 clear P_wk
%                 
%                 % Old Method
%                 for t_no = 1 : length(P.v)
%                     curr_t = (t_no-1)/P.fs;
%                     term2 = (P_0 - P_inf)*exp(-curr_t/(R*C));
%                     term3a = Q.v(1:t_no); term3a = term3a(:);
%                     term3b = exp((((1:t_no)-1)/P.fs)./(R*C)); term3b = term3b(:);
%                     term3 = exp(-1*curr_t/(R*C)) * (1/C) * sum(term3a.*term3b) * (1/P.fs);
%                     temp_P_wk = P_inf + term2 + term3;
%                     P_wk(t_no,1) = temp_P_wk;
%                 end
%                 P_wk = P_wk/133.33;
%                 plot(P.v/133.33), hold on
%                 plot(P_wk, '--k')
                
                [Pr,A,B,Pinf,Tn,Pn]=kreservoir_v10_pc(P.v,length(P.v)-1,282); %collated_data(sim_no).input_data.sim_settings.lvet);
                
                PexcKP = P.v(:) - Pr(:);
                collated_data(sim_no).output_data(domain_el).PexcKP(:,distance_el) = PexcKP;
                
                do_plot = 0;
                if do_plot && domain_el == 1
                    plot(P.v), hold on, plot(P.v(1)+Pexc), plot(P.v(1)+PexcKP)
                    waitfor(gcf)
                end
                
                % Estimate wave speed, c
                ave_radius_cm = 100*mean(D.v/2);
                rho = collated_data(1).input_data.sim_settings.rho;
                k = collated_data(1).input_data.sim_settings.network_spec.k;
                Eh_D0 = (k(1)*exp(k(2)*ave_radius_cm))+k(3); % Eh/D0 (from Mynard's 2015 paper, eqn 3)
                c0_squared = (2/3)*(Eh_D0/(rho/1000)); % from Mynard's 2015 paper, eqn 3, converts rho from kg/m3 to g/cm3
                c = sqrt(c0_squared)/100;  % converts from cm/s to m/s.
                
                % Perform wave separation to obtain forward and backward P and U waves
                dP_f = (diff(P.v) + (rho*c*diff(U.v)))/2;
                dP_b = (diff(P.v) - (rho*c*diff(U.v)))/2;
                dU_f = (diff(U.v) + (diff(P.v)/(rho*c)))/2;
                dU_b = (diff(U.v) - (diff(P.v)/(rho*c)))/2;
                dU_f = dP_f/(rho*c);
                dU_b = -1*dP_b/(rho*c);
                P_f = cumsum(dP_f);
                P_b = cumsum(dP_b);
                U_f = U.v(1) + cumsum(dU_f);
                U_b = U.v(1) + cumsum(dU_b);
                
                collated_data(sim_no).output_data(domain_el).P_f(:,distance_el) = P_f;
                collated_data(sim_no).output_data(domain_el).P_b(:,distance_el) = P_b;
                collated_data(sim_no).output_data(domain_el).U_f(:,distance_el) = U_f;
                collated_data(sim_no).output_data(domain_el).U_b(:,distance_el) = U_b;
                
                do_plot = 0; rel_dom_el = 7;
                if do_plot && domain_el == rel_dom_el && distance_el == no_distance_els
                    figure('Position', [20,20,1000,600])
                    subplot(2,2,1)
                    plot(P.v-P.v(1), 'k'), hold on, plot(P_f, 'b'), plot(P_b, 'r')
                    subplot(2,2,3)
                    plot(U.v, 'k'), hold on, plot(U_f, 'b'), plot(U_b, 'r')
                end
                
                % Perform wave separation to obtain forward and backward P and U waves
                dP_f = (diff(Pexc) + (rho*c*diff(U.v)))/2;
                dP_b = (diff(Pexc) - (rho*c*diff(U.v)))/2;
                %dP_b = (diff(P.v-PexcKP) - (rho*c*diff(U.v)))/2;
                dU_f = (diff(U.v) + (diff(Pexc)/(rho*c)))/2;
                dU_b = (diff(U.v) - (diff(Pexc)/(rho*c)))/2;
                P_f = cumsum(dP_f);
                P_b = cumsum(dP_b);
                U_f = U.v(1) + cumsum(dU_f);
                U_b = U.v(1) + cumsum(dU_b);
                
                if do_plot && domain_el == rel_dom_el && distance_el == no_distance_els
                    subplot(2,2,2)
                    plot(P.v-P.v(1), 'k'), hold on, plot(P_f, 'b'), plot(P_b, 'r'), plot(P.v-P.v(1)-PexcKP, 'g')
                    subplot(2,2,4)
                    plot(U.v, 'k'), hold on, plot(U_f, 'b'), plot(U_b, 'r')
                    waitfor(gcf)
                end
                
                do_plot = 0;
                if do_plot
                    figure('Position', [100,100,900,400])
                    subplot(1,2,1), plot(U.v), hold on, plot(U_f), plot(U_b)
                    subplot(1,2,2), plot(PPG.v)
                    close all
                end
            end
            
        end
        
    end
end


% save collated data
a = whos('collated_data');
if a.bytes > 1.8e9
    save(PATHS.collated_data, 'collated_data', '-v7.3')
else
    save(PATHS.collated_data, 'collated_data')
end

end

function ppg = estimate_ppg_using_windkessel(P, Q, Pout)

%% Calculate PPG

% Calculate time vectors for input signals
P.t = [0:length(P.v)-1]/P.fs;
Q.t = [0:length(Q.v)-1]/Q.fs;

% Find the resistance to flow further down the arterial tree
temp = P.v-Pout.v; % pressure drop between this segment and end of arterial tree
R = sum(temp)/sum(Q.v); % resistance is mean pressure drop over mean flow

% Find the flow into the more distal part of the arterial tree from this segment
Qout.v = (P.v-Pout.v)./R;  % I = V/R (electrical circuit)
Qout.t = P.t;

% Find the volume stored in the arterial segment
Volvb.t = Q.t;
const = 0;
Volvb.v = const + cumsum(Q.v) - cumsum(Qout.v);  % volume stored is the difference between inflow and outflow

ppg.t = Volvb.t;
% ppg.v = normwave(Volvb.v);
ppg.v = Volvb.v;

end

function norm_wave = normwave(orig_wave)

norm_wave = orig_wave;

norm_wave = (norm_wave-min(norm_wave))./range(norm_wave);

end

function static_wave = calc_static_wave(orig_wave)

orig_wave = orig_wave(:);

% Calculate expected position of next point
mean_diff = mean([diff(orig_wave(end-1:end)), diff(orig_wave(1:2))]);
next_point = orig_wave(end) + mean_diff;

% If the wave was static, then this next point would be equal to the first point
% So, we can make the wave static by making this next point equal to the first
temp = linspace(0,orig_wave(1)-next_point, length(orig_wave)); temp = temp(:);
static_wave = orig_wave + temp;

end

function pres_cost_function = calculate_pres_cost_function(P, Q, pout, R, tdia, input_vars)

% setup
C = input_vars(1);
pout = input_vars(2);

% calculate reservoir pressure
Pres = calc_pres(P,Q,pout,R, C);

% assess performance
t = [0:length(Pres.v)-1]/Pres.fs;
duration_dia = t(end) - tdia;
rel_els = t>= (tdia+(1/3)*duration_dia);
%rel_els = t>= tdia;
pres_diffs = P.v(rel_els) - Pres.v(rel_els);

% plot(t, P.v, 'b'), hold on
% plot(t(rel_els), P.v(rel_els), 'r'), hold on, plot(t(rel_els), Pres.v(rel_els), 'k'), title(num2str(C))
% waitfor(gcf)

pres_cost_function = sqrt(mean(pres_diffs.^2));

end

function Pres = calc_pres(P,Q,pout,R, C)

RC = R*C;
t = [0:length(P.v)-1]/P.fs; t = t(:);
dt=1/P.fs;
qse=cumtrapz(Q.v.*exp(t/RC))*dt;
Pres.v = (exp(-t/RC).*(qse/C)) + (exp(-t/RC)*(P.v(1) - pout)) + pout;
Pres.fs = P.fs;
end
