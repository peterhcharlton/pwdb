function inflow = AorticFlowWave(params)
% AORTICFLOWWAVE generates an exemplary aortic flow wave, with optional
% input parameters specifying its properties
%
%               AorticFlowWave
%
%   Inputs:     - params, an optional structure of input parameters
%
%   Outputs:    - inflow, a structure containing the waveform and its
%   properties.
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
%  Copyright (C) 2018  King's College London
%
% Contributed to by: Marie Willemet, Jordi Alastruey, and Peter H. Charlton
% v.1.0

%% ==== Setup input params (adding default and derived values where necessary)
if nargin < 1
    params = struct;
end
inflow = setup_params(params);

%% ==== Find Template Flow Wave
template_flow = create_template_inflow_waveform(inflow);

%% ==== Find Adjusted Flow Wave
inflow = adjust_template_flow(inflow, template_flow);

%% ==== Adjust Sampling Frequency
inflow = adjust_sampling_frequency(inflow);

%% ==== Calculate Characteristics of Inflow wave
inflow = calc_characteristics_inflow_wave(inflow);

%% ==== Plot Figures
plot_figures(inflow);    

end

function inflow = setup_params(params)

%% Setup

% store the user-specified input params
inflow.input_params = params; clear params

% shut figures unless plotting multiple flow waveforms
if ~sum(strcmp(fieldnames(inflow.input_params),'plot_multiple'))
    close all
end

%% Input parameters

% if CO and HR have been specified, but not SV, then calculate SV
curr_params = fieldnames(inflow.input_params);
if sum(strcmp(curr_params, 'HR')) && sum(strcmp(curr_params, 'CO')) && ~sum(strcmp(curr_params, 'SV'))
    inflow.input_params.SV = 1000*inflow.input_params.CO / inflow.input_params.HR;
end

% If either HR or T has been specified, but not the other, then specify it
if sum(strcmp(curr_params, 'HR')) && ~sum(strcmp(curr_params, 'T'))
    inflow.input_params.T = 60/inflow.input_params.HR;
elseif ~sum(strcmp(curr_params, 'HR')) && sum(strcmp(curr_params, 'T'))
    inflow.input_params.HR = 60/inflow.input_params.T;
end

% If LVET has been specified then make sure it is in secs
if sum(strcmp(curr_params, 'LVET')) && inflow.input_params.LVET > 2
    inflow.input_params.LVET = inflow.input_params.LVET/1000;   % convert from ms to secs
end

% If T_Peak_Flow has been specified then make sure it is in secs
if sum(strcmp(curr_params, 'T_Peak_Flow')) && inflow.input_params.T_Peak_Flow > 2
    inflow.input_params.T_Peak_Flow = inflow.input_params.T_Peak_Flow/1000;   % convert from ms to secs
end

% identify any missing parameters
curr_params = fieldnames(inflow.input_params);
req_params = {'fs', 'wave_type', 'HR', 'SV', 'CO', 'T', 'LVET', 'T_Peak_Flow', 'T_sys_upslope_ip', 'Q_sys_upslope_ip', 'T_Min_Flow', 'T_dia', 'Reg_Vol', 'rev_flow', 'contractility_const'};
curr_req_params = intersect(curr_params, req_params);
missing_params = setxor(curr_req_params, req_params); clear curr_params curr_req_params
% moves CO, LVET, T and T_sys_upslope_ip to the end, as they are dependent on other parameters.
for rel_param = req_params
    temp = find(strcmp(missing_params, rel_param{1,1}));
    if ~isempty(temp)
        missing_params = missing_params([1:temp-1, temp+1 : length(missing_params)]);
        missing_params(end+1) = {rel_param{1,1}};
    end
end


% If any parameters are missing then insert the baseline values for them
for param_no = 1 : length(missing_params)
    
    % identify the current missing parameter
    curr_missing_param = missing_params{param_no};
    
    % specify the baseline value of this parameter
    switch curr_missing_param
        case 'HR'
            val = 75;       % in bpm
            val = 75;       % Mynard's, in bpm
        case 'SV'
            val = 83;       % im ml
            val = (6.2*1000/inflow.input_params.HR);  % Mynard's, in ml
        case 'CO'
            val = inflow.input_params.HR*inflow.input_params.SV/1000;       % im l/min
            %val = 6.2;      % Mynard's, in l/min
        case 'wave_type'
            val = 'elderly';  % either 'young' or 'avg' or 'elderly'
            val = 'Mynard'; % Use Mynard2015's aortic root wave
        case 'rev_flow'
            val = 1;        % between 0 and 1 - the proportion of reverse flow to include
        case 'fs'
            val = 1000;     % in Hz
        case 'contractility_const'
            val = 1;        % a multiplication factor applied to the time of the systolic peak
        case 'T'
            val = 60/inflow.input_params.HR;
        case 'T_sys_upslope_ip'
            val = (0.020/0.085)*inflow.input_params.T_Peak_Flow;
            val = (0.010/0.080)*inflow.input_params.T_Peak_Flow; % Mynard's
        case 'Q_sys_upslope_ip'
            val = 0.32;
            % Haven't measured Mynard's - don't think it's used
        case 'LVET'
            val = calc_lvet(inflow);
            val = 0.282;   % Mynard's, s
        case 'T_Peak_Flow'
            val = (0.085/0.290)*inflow.input_params.LVET;
            val = (0.080/0.282)*inflow.input_params.LVET;  % Mynard's
        case 'T_Min_Flow'
            val = (0.310/0.290)*inflow.input_params.LVET;
            val = (0.298/0.282)*inflow.input_params.LVET;  % Mynard's
        case 'T_dia'
            val = (0.330/0.290)*inflow.input_params.LVET;
            val = (0.309/0.282)*inflow.input_params.LVET;  % Mynard's
        case 'Reg_Vol'
            val = 0.73;
            val = 1.2775; % Mynard's, ml
    end
            
    % insert this baseline value
    eval(['inflow.input_params.' curr_missing_param ' = val;'])
    
    clear val
    
end

%% Scale contractility and reverse flow magnitude
contractility_times = {'T_sys_upslope_ip', 'T_Peak_Flow'};
for s = 1 : length(contractility_times)
    eval(['inflow.input_params.' contractility_times{s} ' = inflow.input_params.contractility_const * inflow.input_params.' contractility_times{s} ';']);
end
inflow.input_params.Reg_Vol = inflow.input_params.rev_flow * inflow.input_params.Reg_Vol;

%% Input settings

% identify any missing settings
curr_params = fieldnames(inflow.input_params);
req_params = {'do_plot','save_plot', 'plot_multiple', 'plot_name', 'file_path'};
curr_req_params = intersect(curr_params, req_params);
missing_params = setxor(curr_req_params, req_params); clear curr_params req_params

% If any parameters are missing then insert the baseline values for them
for param_no = 1 : length(missing_params)
    
    % identify the current missing parameter
    curr_missing_param = missing_params{param_no};
    
    % specify the baseline value of this parameter
    switch curr_missing_param
        case 'do_plot'
            val = true;     % whether or not to plot the generated wave
        case 'save_plot'
            val = false;
        case 'plot_multiple'
            val = false;
        case 'plot_name'
            val = 'AorticFlowWave_plot';
        case 'file_path'
            val = '/Users/petercharlton/Google Drive/Work/Code/AorticFlowWave/AorticFlowWave manual/Figures/';
            if ~exist(val, 'dir')
                val = uigetdir;
            end
    end
            
    % insert this baseline value
    eval(['inflow.input_settings.' curr_missing_param ' = val;'])
    
    clear val
    
end

% move user specified settings to this new structure
for param_no = 1 : length(curr_req_params)
    
    % identify the current missing parameter
    curr_setting = curr_req_params{param_no};
            
    % insert this setting
    eval(['inflow.input_settings.' curr_setting ' = inflow.input_params.' curr_setting ';'])
    
    % remove from the params structure
    inflow.input_params = rmfield(inflow.input_params, curr_setting);
    
end

%% Put fields in alphabetical order
inflow.input_params = orderfields(inflow.input_params);
inflow.input_settings = orderfields(inflow.input_settings);

end

function LVET = calc_lvet(inflow)

% use the function provided in the following article to derive LVET:
%
%    Weissler et al, "Relationships between left ventricular ejection time,
%    stroke volume, and heart rate in normal individuals and patients with
%    cardiovascular disease", American Heart Journal, vol. 62, 1961.
%    DOI: 10.1016/0002-8703(61)90403-3 

% Function:
%    LVET = 0.266 + 0.0011*(SV - 82) - 0.0009*(HR - 73)
% where LVET is in secs, SV [ml] is stroke volume, and HR [bpm] is heart rate.

% If a LVET has been specified as an input, then take this
if ~sum(strcmp(fieldnames(inflow.input_params), 'LVET'))
   
    % otherwise, calculate using this formula
    LVET = 0.266 + 0.0011*(inflow.input_params.SV - 82) - 0.0009*(inflow.input_params.HR - 73);   % in secs

else
    LVET = inflow.input_params.LVET;
end


end

function template_flow = create_template_inflow_waveform(inflow)

%% === Setup constants
do_change = 1;
if do_change
    prop = 1.2;
    prop2 = 1.1;
    prop3 = 1.1;
    prop4 = 1.04;
    prop5 = 1.05;
    prop6 = 1.03;
else
    [prop, prop2,prop3,prop4,prop5, prop6] = deal(1);
end

% timings
template_flow.T = 1;
template_flow.ti_0   = (1/prop)*0.024*template_flow.T; %Time of inflection point in systolic increase
template_flow.tmax_0 = 0.08*template_flow.T; %Time of systolic peak
template_flow.ti2_0  = (1/prop5)*prop2*0.17*template_flow.T;  %0.150; Time of inflection point in systolic decrease (old only)
template_flow.ti3_0  = prop4*prop2*0.21*template_flow.T;  %0.150; Time of inflection point in systolic decrease (old only)
template_flow.ts_0   = 0.290*template_flow.T; %Time of start of dicrotic notch
template_flow.tmin_0 = 0.310*template_flow.T; %Time of minimum flow during dicrotic notch  %%% THIS HAS BEEN ADJUSTED from 0.30
template_flow.td_0   = 0.330*template_flow.T; %Time of start of diastole

% flows
template_flow.Qmax = 1;    %Systolic peak flow
template_flow.Qi   = 0.367; %Flow at inflection point on systolic upslope
template_flow.Qi2  = prop5*prop3*(1/prop2)*0.647; %Flow at inflection point in systolic decrease (old)
template_flow.Qi3  = prop6*(1/prop3)*(1/prop2)*0.496; %Flow at second inflection point in systolic decrease (old)

% find approximate Qmin (during reverse flow)
template_flow.Qmin = -0.1; %Minimum flow during dicrotic notch
% % Adjust if the regurgitation volume is provided
% if sum(strcmp(fieldnames(params),'Reg_Vol'))
%     desired_ratio_reverse_to_forward_flow = params.Reg_Vol/(params.SV+params.Reg_Vol);
%     curr_reverse_flow = abs(0.5*(template_flow.td_0-template_flow.ts_0)*template_flow.Qmin);
%     curr_forward_flow = 0.5*template_flow.ts_0*template_flow.Qmax;
%     approx_curr_ratio_reverse_to_forward_flow = curr_reverse_flow/curr_forward_flow;
%     scale_factor = desired_ratio_reverse_to_forward_flow/approx_curr_ratio_reverse_to_forward_flow;
%     template_flow.Qmin = template_flow.Qmin*scale_factor;
% end

% round timings to  right sizes of time vector
ti   = round(template_flow.ti_0  *1000)./1000;
tmax = round(template_flow.tmax_0*1000)./1000;
ti2  = round(template_flow.ti2_0 *1000)./1000;
ti3  = round(template_flow.ti3_0 *1000)./1000;
ts   = round(template_flow.ts_0  *1000)./1000;
tmin = round(template_flow.tmin_0*1000)./1000;
td   = round(template_flow.td_0  *1000)./1000;

% time step
dt = 1/1000;

%% === Calculate template using Marie's piecewise polynomial

% switch params.wave_type   % young, avg, or elderly
%     case 'young'
%         template_flow.age = 'young';
%         template_flow.v = Get_parametric_young(dt, template_flow.Qi, template_flow.Qmax, template_flow.Qmin, ti, tmax, tmin, ts, td);
%     case 'avg'
%         template_flow.age = 'avg';
%         young = Get_parametric_young(dt, template_flow.Qi, template_flow.Qmax, template_flow.Qmin, ti, tmax, tmin, ts, td);
%         elderly = Get_parametric_old(dt, template_flow.Qi, template_flow.Qmax, template_flow.Qmin, template_flow.Qi2, ti, tmax, tmin, ti2, ts, td);
%         template_flow.v = (young+elderly)./2; clear young elderly
%     case 'elderly'
%         template_flow.age = 'elderly';
%         template_flow.v = Get_parametric_old(dt, template_flow.Qi, template_flow.Qmax, template_flow.Qmin, template_flow.Qi2, ti, tmax, tmin, ti2, ts, td);
%     case 'impulse'
%         template_flow.age = 'impulse';
%         durn_of_impulse = 0.05*template_flow.T;
%         durn_of_impulse_samples = round(durn_of_impulse*(1/dt));
%         durn_of_beat = template_flow.T;
%         durn_of_beat_samples = round(durn_of_beat*(1/dt));
%         template_flow.v = [0, ones(1,durn_of_impulse_samples), zeros(1,durn_of_beat_samples - durn_of_impulse_samples-1)];
% end

if strcmp(inflow.input_params.wave_type, 'elderly')
    template_flow.age = 'elderly';
    template_flow.v = Get_parametric_pete(dt, template_flow.Qi, template_flow.Qmax, template_flow.Qmin, template_flow.Qi2, ti, tmax, tmin, ti2, ts, td, template_flow.Qi3, ti3);
elseif strcmp(inflow.input_params.wave_type, 'Mynard')
    template_flow.age = 'young';
    template_flow.v = Get_Mynard2015_flow_wave;
end
template_flow.t = [0:(length(template_flow.v)-1)]*(1/1000);

end

function inflow = adjust_template_flow(inflow, template_flow)
%% === Adjust template to give desired timings

% setup
dt = 1/inflow.input_params.fs;

% Find sytolic upslope
deb = 1; [~, fin] = max(template_flow.v);
sys_upslope.v = template_flow.v(deb:fin);
sys_upslope.t = linspace(0, inflow.input_params.T_Peak_Flow, length(sys_upslope.v));

% Find systolic downslope
[~, temp] = max(template_flow.v); deb = temp+1; clear temp
fin = find(template_flow.v(1:end-1) > 0 & template_flow.v(2:end) <= 0); fin = fin(1);
sys_downslope.v = template_flow.v(deb:fin);
sys_downslope.t = linspace(sys_upslope.t(end)+dt, inflow.input_params.LVET, length(sys_downslope.v));

% Find reverse flow
deb = find(template_flow.v(1:end-1) > 0 & template_flow.v(2:end) <= 0); deb = deb(1)+1;
fin = find(template_flow.v ~= 0, 1, 'last');
reverse.v = template_flow.v(deb:fin);
reverse.t = linspace(sys_downslope.t(end)+dt, inflow.input_params.T_dia, length(reverse.v));

% Find diastolic flat line
no_els = round((inflow.input_params.T - (inflow.input_params.T_dia+dt))/dt);
diastolic.t = linspace(inflow.input_params.T_dia+dt, inflow.input_params.T, no_els);
diastolic.v = zeros(size(diastolic.t));

% concatenate to give waveform
mod_flow.t = [sys_upslope.t, sys_downslope.t, reverse.t, diastolic.t];
mod_flow.v = [sys_upslope.v, sys_downslope.v, reverse.v, diastolic.v];

% resample to give constant fs, without irregular spacing at joins
inflow.t = (0 : floor(mod_flow.t(end)*inflow.input_params.fs))/inflow.input_params.fs;
inflow.v = interp1(mod_flow.t, mod_flow.v, inflow.t);

if sum(strcmp(fieldnames(inflow.input_params), 'Reg_Vol'))
    
    % Scale to give desired reverse flow volume
    scale_factor = inflow.input_params.Reg_Vol/abs(sum(inflow.v(inflow.v<0)*dt));
    inflow.v(inflow.v<0) = inflow.v(inflow.v<0)*scale_factor;  % flow is now in ml/s
    
    % Scale to give desired stroke volume
    scale_factor = (inflow.input_params.SV+inflow.input_params.Reg_Vol)/(sum(inflow.v(inflow.v>0))*dt);
    inflow.v(inflow.v>0) = inflow.v(inflow.v>0)*scale_factor;
    
else
    % Scale to give desired stroke volume
    scale_factor = inflow.input_params.SV/(sum(inflow.v)*dt);
    inflow.v = inflow.v*scale_factor;  % flow is now in ml/s
    
end

% Convert to required units
CO_in_m3_per_sec = inflow.input_params.CO/(60*1000); %(convert from l/min)
inflow.v = inflow.v./mean(inflow.v).*CO_in_m3_per_sec;

% Scale to give desired Qmax (if specified). This overrides the CO value
if sum(strcmp(fieldnames(inflow.input_params), 'Qmax'))
    scale_factor = inflow.input_params.Qmax/max(inflow.v);
    inflow.v(inflow.v>0) = inflow.v(inflow.v>0)*scale_factor;
end

% Specify dQdt (if specified)
if sum(strcmp(fieldnames(inflow.input_params), 'dQdt_max'))
    sf = 1;
    options = optimset('MaxFunEvals',200, 'display', 'off');
    dQdt_cost_function = @(sf)calculate_dQdt_cost_function(inflow, inflow.input_params.dQdt_max, sf);
    [optimal_sf, ~] = fminsearch(dQdt_cost_function, sf, options);
    inflow = calculate_new_inflow_from_sf(inflow, optimal_sf);
    current_dQdt_max = calc_dQ_dt_max(inflow);
end

end

function dQdt_cost_function = calculate_dQdt_cost_function(inflow, desired_dQdt_max, sf)

new_inflow = calculate_new_inflow_from_sf(inflow, sf);

current_dQdt_max = calc_dQ_dt_max(new_inflow);

% fprintf("\n - SF: %f,    dQdt_max: %f", sf, current_dQdt_max);

dQdt_cost_function = abs(current_dQdt_max-desired_dQdt_max);

end

function current_dQ_dt_max = calc_dQ_dt_max(inflow)
current_dQ_dt_max = max(diff(inflow.v))/(1/inflow.input_params.fs);
end

function new_inflow = calculate_new_inflow_from_sf(inflow, sf)

new_inflow = inflow;

deb = 1; [~, fin] = max(inflow.v); % identify systolic upslope
old_deriv = [0, diff(inflow.v(deb:fin))];
[~,ms] = max(old_deriv);
old_deriv_sum = sum(old_deriv);
% old way
N1 = ms;
N2 = fin-deb+1-N1;
scale_vector = ones(1,N1+N2).*[linspace(1,sf,N1), linspace(sf,1,N2)];
% % new way
% N1 = ms;
% N2 = round((0.5*(fin-deb))-ms);
% N3 = fin-deb+1-N1-N2;
% scale_vector = ones(1,N1+N2+N3).*[linspace(1,sf,N1), sf*ones(1,N2), linspace(sf,1,N3)];
% applying scaling to derivative
new_deriv = old_deriv.*scale_vector;

% to avoid shifting the point of max dQdt:
new_deriv(new_deriv>new_deriv(ms)) = new_deriv(ms);

new_deriv = old_deriv_sum*new_deriv/sum(new_deriv);
new_inflow.v = [inflow.v(1)+cumsum(new_deriv), inflow.v(fin+1:end)];

end

function Q = Get_parametric_pete(dt, Qi, Qmax, Qmin,Qi2, ti,tmax,tmin,ti2, ts, td, Qi3, ti3)

%% Find Q1 (systolic uplslope, fourth order)

t_matrix = [
    0           0           0           0       1;  % Q1(0)     = 0
    ti^4        ti^3        ti^2        ti      1;  % Q1(ti)    = Qi
    4*3*ti^2    3*2*ti      2*1         0       0;  % Q1''(ti)  = 0
    4*tmax^3    3*tmax^2    2*tmax      1       0;  % Q1'(tmax) = 0
    tmax^4      tmax^3      tmax^2      tmax    1   % Q1(tmax)  = Qmax
    ];

q_matrix = [
    0;
    Qi;
    0;
    0;
    Qmax
    ];

%Solve for coefficients t_matrix*coeffs=q_matrix: coeffs = t_matrix\q_matrix;
Q1_coeffs = t_matrix\q_matrix;

% %% Find Q2 (systolic downslope, second order)
% 
% t_matrix = [
%     tmax^2      tmax    1;  % Q2(tmax)  = Qmax
%     2*tmax      1       0;  % Q2'(tmax) = 0
%     ts^2        ts      1;  % Q2(ts)    = 0
%     ];
% 
% q_matrix = [
%     Qmax;
%     0;
%     0;
%     ];
% 
% %Solve for coefficients t_matrix*coeffs=q_matrix: coeffs = t_matrix\q_matrix;
% Q2_coeffs = t_matrix\q_matrix;
% 
% %% Find Q2 (systolic downslope with inflection point, fourth order)
% 
% t_matrix = [
%     tmax^4      tmax^3      tmax^2      tmax    1;  % Q2(tmax)  = Qmax
%     4*tmax^3    3*tmax^2    2*tmax      1       0;  % Q2'(tmax) = 0
%     ti2^4       ti2^3       ti2^2       ti2     1;  % Q2(ti2)   = Qi2
%     4*3*ti2^2   3*2*ti2     2*1         0       0;  % Q2''(ti2) = 0
%     ts^4        ts^3        ts^2        ts      1;  % Q2(ts)    = 0
%     ];
% 
% q_matrix = [
%     Qmax;
%     0;
%     Qi2;
%     0;
%     0;
%     ];
% 
% %Solve for coefficients t_matrix*coeffs=q_matrix: coeffs = t_matrix\q_matrix;
% Q2_coeffs = t_matrix\q_matrix;

%% Find Q2 (systolic downslope with two inflection points, sixth order)

t_matrix = [
    tmax^6      tmax^5      tmax^4      tmax^3      tmax^2      tmax    1;  % Q2(tmax)  = Qmax
    6*tmax^5    5*tmax^4    4*tmax^3    3*tmax^2    2*tmax      1       0;  % Q2'(tmax) = 0
    ti2^6       ti2^5       ti2^4       ti2^3       ti2^2       ti2     1;  % Q2(ti2)   = Qi2
    6*5*ti2^4   5*4*ti2^3   4*3*ti2^2   3*2*ti2     2*1         0       0;  % Q2''(ti2) = 0
    ti3^6       ti3^5       ti3^4       ti3^3       ti3^2       ti3     1;  % Q2(ti3)   = Qi3
    6*5*ti3^4   5*4*ti3^3   4*3*ti3^2   3*2*ti3     2*1         0       0;  % Q2''(ti3) = 0
    ts^6        ts^5        ts^4        ts^3        ts^2        ts      1;  % Q2(ts)    = 0
    ];

q_matrix = [
    Qmax;
    0;
    Qi2;
    0;
    Qi3;
    0;
    0;
    ];

%Solve for coefficients t_matrix*coeffs=q_matrix: coeffs = t_matrix\q_matrix;
Q2_coeffs = t_matrix\q_matrix;

%% Find Q3 (reverse flow, third order)

%Q2_deriv_ts = 2*Q2_coeffs(1).*ts + Q2_coeffs(2);

t_matrix = [
    ts^3        ts^2        ts      1;  % Q3(ts) = 0
    % 4*ts^3      3*ts^2      2*ts        1       0;  % Q3'(ts) = Q2'(ts)
    3*tmin^2    2*tmin      1       0;  % Q3'(tmin) = 0
    tmin^3      tmin^2      tmin^1  1;  % Q3(tmin)  = Qmin
    td^3        td^2        td      1;  % Q3(td)    = 0
    %3*td^2      2*td        1       0;  % Q3'(td)   = 0
    ];

q_matrix = [
    0;
    %Q2_deriv_ts;
    0;
    Qmin
    0;
    ];

%Solve for coefficients t_matrix*coeffs=q_matrix: coeffs = t_matrix\q_matrix;
Q3_coeffs = t_matrix\q_matrix;

%% Calculate template flow curve

dt=1e-3;

% Q1, t1
t1=[0:dt:tmax];
Q1 = Q1_coeffs(1).*t1.^4 + Q1_coeffs(2).*t1.^3 + Q1_coeffs(3).*t1.^2 +Q1_coeffs(4).*t1 + Q1_coeffs(5);

% Q2, t2
t2=[tmax+dt:dt:ts];
Q2 = Q2_coeffs(1)*t2.^6 + Q2_coeffs(2)*t2.^5 + Q2_coeffs(3).*t2.^4 + Q2_coeffs(4).*t2.^3 + Q2_coeffs(5).*t2.^2 +Q2_coeffs(6).*t2 + Q2_coeffs(7);
% Q2 = Q2_coeffs(1).*t2.^4 + Q2_coeffs(2).*t2.^3 + Q2_coeffs(3).*t2.^2 +Q2_coeffs(4).*t2 + Q2_coeffs(5);
%Q2 = Q2_coeffs(1).*t2.^2 + Q2_coeffs(2).*t2 + Q2_coeffs(3);

% Q3, t3
t3 = [ts+dt:dt:td];
% Q3 = Q3_coeffs(1).*t3.^4 +Q3_coeffs(2).*t3.^3 + Q3_coeffs(3).*t3.^2 + Q3_coeffs(4).*t3 + Q3_coeffs(5);
Q3 = Q3_coeffs(1).*t3.^3 +Q3_coeffs(2).*t3.^2 + Q3_coeffs(3).*t3 + Q3_coeffs(4);

% Q4, t4
t4 = [td+dt:dt:1.0];
Q4 = t4.*0;

%% Construct template curve
t = [t1,t2,t3,t4];
Q = [Q1,Q2,Q3,Q4];

end

function Q = Get_Mynard2015_flow_wave

fs = 1000; % Hz
% Data provided by Mynard:
Q = [0.975173175413631;0.973779853988154;0.972395614614815;0.970995717731558;0.969511227019400;0.967822711189793;0.965786748864051;0.963283288352364;0.960257180172596;0.956733614590217;0.952804019043309;0.948592861750212;0.944220671600988;0.939775612934077;0.935299487236573;0.930788206519929;0.926202977984071;0.921487672252196;0.916588056259775;0.911469209683415;0.906128003680296;0.900598361098229;0.894948479110418;0.889270767746236;0.883666691560787;0.878229797375804;0.873030633851194;0.868106794857485;0.863460040023719;0.859060747248515;0.854858328407404;0.850795093069735;0.846820567387548;0.842903445940403;0.839039020574192;0.835250902526769;0.831586960403603;0.828110482041099;0.824888470942132;0.821979531423770;0.819423828948797;0.817237101401127;0.815409739068083;0.813910764962756;0.812695446475572;0.811714543871594;0.810923016496613;0.810286366332797;0.809783558074059;0.809406383225149;0.809155982868884;0.809049312351120;0.822723666912791;1.02233385165230;1.96658629314415;4.58560140254054;9.85427558739735;18.3723567331533;30.1314774777907;44.5666042279459;60.7946355212591;77.8878867053600;95.0770252268905;111.846390198646;127.934762788699;143.278424134962;157.935102195756;172.016187493348;185.639837609637;198.905789610112;211.886315485302;224.626247569543;237.146494154128;249.448001830952;261.515308275138;273.320033907046;284.824893763383;295.988445536841;306.770268021839;317.135901678616;327.060806548574;336.532774620005;345.552564122539;354.132876926804;362.296079123497;370.071215649353;377.490884042716;384.588435958977;391.395812923550;397.942143214896;404.253068924458;410.350659695562;416.253710063464;421.978206529741;427.537777920319;432.943995528724;438.206455565240;443.332643996137;448.327642700327;453.193777638760;457.930329008757;462.533418553688;466.996162750017;471.309138243523;475.461155370642;479.440285080944;483.235041551993;486.835593235815;490.234863026527;493.429385478842;496.419814845701;499.211018863089;501.811744031372;504.233891257801;506.491488218617;508.599479649849;510.572474215038;512.423584988732;514.163481563606;515.799739601955;517.336534070953;518.774681411257;520.111998663707;521.343917882841;522.464274200170;523.466176452167;524.342869998804;525.088510984828;525.698787922101;526.171347706221;526.506006380627;526.704747502486;526.771530539299;526.711946609556;526.532768110292;526.241442186613;525.845576190021;525.352457228893;524.768639135078;524.099619940398;523.349622557002;522.521481755218;521.616632378229;520.635187375853;519.576089782668;518.437320238385;517.216140141572;515.909351597176;514.513556936185;513.025401511467;511.441787659971;509.760050567731;507.978090122726;506.094456391599;504.108390059734;502.019820234977;499.829327332228;497.538078870921;495.147747545238;492.660421408567;490.078515604694;487.404693896395;484.641806405646;481.792847669950;478.860936563018;475.849317075378;472.761376647828;469.600676913928;466.370990500281;463.076337047776;459.721011868244;456.309601576609;452.846982508583;449.338299567697;445.788925134535;442.204399603282;438.590356794328;434.952438775355;431.296205409574;427.627044124906;423.950085234953;420.270128045542;416.591580412164;412.918414146340;409.254137683546;405.601785901964;401.963925786966;398.342675776211;394.739736069239;391.156426918743;387.593731908661;384.052343427062;380.532707909337;377.035068903213;373.559506534790;370.105972482048;366.674320045660;363.264329314927;359.875727739749;356.508206630724;353.161434222293;349.835065959949;346.528752629132;343.242146850871;339.974908347998;336.726708255579;333.497232626731;330.286185183115;327.093289287159;323.918289074735;320.760949682792;317.621056533014;314.498413683235;311.392841324680;308.304172575322;305.232249788152;302.176920648979;299.138034374586;296.115438333779;293.108975399232;290.118482297792;287.143789164435;284.184720425785;281.241097063450;278.312740302484;275.399476098279;272.501140627695;269.617586923616;266.748691925274;263.894363582487;261.054547078738;258.229230595293;255.418450710163;252.622296983934;249.840915478514;247.074511117000;244.323348894336;241.587754005714;238.868110988727;236.164861986677;233.478504240701;230.809586911346;228.158707317984;225.526506668447;222.913665332752;220.320897695917;217.748946607613;215.198577433230;212.670571703443;210.165720358625;207.684816590472;205.228648295138;202.797990168838;200.393595496125;198.016187701204;195.666451751727;193.345025520643;191.052491223716;188.789367057162;186.556099160881;184.353054028231;182.180511473288;180.038658252078;177.927582416260;175.847268457468;173.797593279242;171.778323012440;169.789110670159;167.829494620601;165.898897841200;163.996627905452;162.121877644825;160.273726422047;158.451141948429;156.652982576222;154.877999996753;153.124842275710;151.392057157989;149.678095575415;147.981315291308;146.299984616011;144.632286127094;142.976320327190;141.330109171368;139.691599394906;138.058665571552;136.429112832201;134.800679174373;133.171037294589;131.537795878214;129.898500285375;128.250629578665;126.591552056964;124.918513546968;123.228701286997;121.519279660594;119.787396713776;118.030173344769;116.244685038415;114.427943132719;112.576879319406;110.688334247474;108.759049349989;106.785660427503;104.764691739498;102.692549897179;100.565517355309;98.3797456114618;96.1312483447445;93.8158947430382;91.4294032575011;88.9673360283976;86.4250942542448;83.7979148138351;81.0808684833106;78.2688601109401;75.3566311223238;72.3387647348645;69.2096942693344;65.9637149614027;62.5949996961145;59.0976191088091;55.4655665094076;51.6927880876767;47.7732188551761;43.7008247436815;39.4696510842022;35.0738773265343;30.5078774208866;25.7662847630790;20.8440598589074;15.7365576638450;10.4395896586277;4.94948395886629;-0.735422487056177;-6.61006110317804;-12.6622828791241;-18.8741551014073;-25.2223446279507;-31.6774795652257;-38.2025309338117;-44.7503098863757;-51.2601554986592;-57.6538083796666;-63.8303716876630;-69.6602108183708;-74.9777007689949;-79.5729897866242;-83.1835808216337;-85.4878655563827;-86.1053603202930;-84.6131464119407;-80.5954995994022;-73.7514724136509;-64.0798633763130;-52.1091746662620;-39.0192868202663;-26.4236251384105;-15.7941270663738;-7.91285846891109;-2.74469604605859;0.293602435264025;1.98007377814065;2.92205567366193;3.49397045416749;3.89008117479282;4.14089208439632;4.20614627486643;4.04493717059770;3.65230704785446;3.06654618725146;2.35815282267805;1.61164891584995;0.908193676387704;0.312993750191313;-0.131483838549379;-0.407114525746058;-0.517743427869374;-0.484134357520946;-0.337500761404327;-0.112833084082376;0.156872389987909;0.444636538227652;0.731345027437622;1.00549206445160;1.26133764803882;1.49627099108680;1.70822087230301;1.89379163630702;2.04750482337735;2.16218836931823;2.23026692820176;2.24552638787934;2.20487245926461;2.10966251151507;1.96632527713331;1.78615156071780;1.58430251438172;1.37821353006803;1.18565666114217;1.02275857288413;0.902257504376651;0.832227903567992;0.815419325838948;0.849258943000078;0.926471019680050;1.03618649105827;1.16536244555750;1.30031059679785;1.42814387876862;1.53798635774857;1.62184185248123;1.67507168268247;1.69648274308898;1.68806788064185;1.65446889571512;1.60224849115926;1.53906274445581;1.47282227801065;1.41092024226895;1.35959019164763;1.32343840995741;1.30517473525675;1.30554512513065;1.32344993220714;1.35621592404900;1.39997894562081;1.45012864222408;1.50176691865472;1.55013711262050;1.59098990128905;1.62086310064297;1.63726410106609;1.63875430171567;1.62494357247146;1.59640897161454;1.55455559570685;1.50143877623741;1.43956630142906;1.37169744439211;1.30065285302460;1.22914600799649;1.15964393278636;1.09426168755604;1.03469201194843;0.982169473485648;0.937466421126849;0.900916576754242;0.872461092336484;0.851711324911019;0.838022391974301;0.830571715975153;0.828437194426511;0.830670286624854;0.836360134370805;0.844685779213896;0.854954550765467;0.866625729579164;0.879319585335726;0.892812809676451;0.907022159158339;0.921978757821491;0.937795948512687;0.954633805326337;0.972663416953261;0.992033828224525;1.01284410678493;1.03512242076436;1.05881332183346;1.08377368490528;1.10977702323459;1.13652523530348;1.16366629882222;1.19081604552159;1.21758195002861;1.24358685111199;1.26849068155330;1.29200857560822;1.31392414863640;1.33409737863406;1.35246646033510;1.36904423843737;1.38390978338493;1.39719607347560;1.40907499283471;1.41974096642944;1.42939455153875;1.43822720054384;1.44640822057111;1.45407470397354;1.46132491496614;1.46821531763466;1.47476114298754;1.48094013869504;1.48669894093813;1.49196136425231;1.49663782766939;1.50063512319714;1.50386578038738;1.50625637879319;1.50775429615074;1.50833254006601;1.50799248071200;1.50676446792542;1.50470646655165;1.50190096972003;1.49845054448812;1.49447242448339;1.49009258936109;1.48543976285306;1.48063972419857;1.47581026734594;1.47105706521158;1.46647060968185;1.46212430881871;1.45807373738226;1.45435696080351;1.45099579024122;1.44799778032830;1.44535875286724;1.44306561926812;1.44109928056824;1.43943740433552;1.43805690958957;1.43693603064433;1.43605587501926;1.43540143560977;1.43496206010104;1.43473141839304;1.43470703956439;1.43488951227804;1.43528145584444;1.43588637346798;1.43670749517974;1.43774670667983;1.43900364339042;1.44047500816853;1.44215414823655;1.44403090376470;1.44609171884197;1.44831998672798;1.45069658631969;1.45320055643842;1.45580984905350;1.45850210192629;1.46125537497052;1.46404880206793;1.46686312043640;1.46968105181609;1.47248752262826;1.47526972300285;1.47801701615651;1.48072071934780;1.48337378498477;1.48597041509377;1.48850564418889;1.49097492473993;1.49337374632552;1.49569731442067;1.49794030833073;1.50009673043004;1.50215985127527;1.50412224791871;1.50597592630257;1.50771251338018;1.50932350081571;1.51080051987820;1.51213562656644;1.51332157683547;1.51435207395019;1.51522197314722;1.51592743274461;1.51646600506683;1.51683666492955;1.51703977759086;1.51707701165433;1.51695120544616;1.51666619753251;1.51622663328266;1.51563775984886;1.51490522148164;1.51403486597555;1.51303257136306;1.51190409986457;1.51065498366048;1.50929044472432;1.50781534874546;1.50623419114788;1.50455111150159;1.50276993125144;1.50089420884152;1.49892730597802;1.49687245884812;1.49473284859790;1.49251166618081;1.49021216766126;1.48783771724494;1.48539181653539;1.48287811963546;1.48030043462400;1.47766271257347;1.47496902595788;1.47222353890944;1.46943047221929;1.46659406611873;1.46371854366174;1.46080807705649;1.45786675868910;1.45489857791480;1.45190740409148;1.44889697580608;1.44587089572787;1.44283263024876;1.43978551276761;1.43673274933238;1.43367742530605;1.43062251176056;1.42757087045167;1.42452525641417;1.42148831753469;1.41846259070050;1.41545049448519;1.41245431861558;1.40947621072178;1.40651816109149;1.40358198630023;1.40066931265464;1.39778156042536;1.39491992977262;1.39208538917569;1.38927866693981;1.38650024618245;1.38375036357629;1.38102901198697;1.37833594693885;1.37567069659738;1.37303257478828;1.37042069645756;1.36783399494897;1.36527124044948;1.36273105902113;1.36021195164343;1.35771231280773;1.35523044825332;1.35276459156949;1.35031291945742;1.34787356556832;1.34544463293838;1.34302420511081;1.34061036425278;1.33820123324222;1.33579484431388;1.33338926016340;1.33098259456667;1.32857300400822;1.32615867997684;1.32373784447988;1.32130875074699;1.31886969005915;1.31641900407923;1.31395510053223;1.31147646950788;1.30898169807164;1.30646948191180;1.30393863372965;1.30138808866182;1.29881690734648;1.29622427726826;1.29360951296356;1.29097205556370;1.28831147197233;1.28562745386192;1.28291981639856;1.28018849649590;1.27743355035738;1.27465515011932;1.27185357956203;1.26902922899795;1.26618258951974;1.26331424686165;1.26042487506438;1.25751523010908;1.25458614360014;1.25163851651538;1.24867331300571;1.24569155419815;1.24269431200759;1.23968270293181;1.23665788186434;1.23362103596021;1.23057337859285;1.22751614349328;1.22445057906192;1.22137794288460;1.21829949650900;1.21521650054114;1.21213021004095;1.20904187013147;1.20595271174608;1.20286394751884;1.19977676786138;1.19669233729080;1.19361179106446;1.19053623215835;1.18746672858958;1.18440431111538;1.18134997126825;1.17830465973955;1.17526928511395;1.17224471292174;1.16923176504156;1.16623121941054;1.16324381003570;1.16027022726762;1.15731111829750;1.15436708785409;1.15143869904904;1.14852647435796;1.14563089670433;1.14275241057533;1.13989142319184;1.13704830578597;1.13422339494885;1.13141699401049;1.12862937441509;1.12586077710418;1.12311141392642;1.12038146909308;1.11767110067176;1.11498044212579;1.11230960387566;1.10965867487975;1.10702772420934;1.10441680261584;1.10182594408854;1.09925516737407;1.09670447746561;1.09417386703982;1.09166331783351;1.08917280196709;1.08670228318944;1.08425171804785;1.08182105698174;1.07941024532497;1.07701922423083;1.07464793150520;1.07229630237125;1.06996427014828;1.06765176686126;1.06535872378557;1.06308507192919;1.06083074245840;1.05859566706718;1.05637977831226;1.05418300991596;1.05200529703488;1.04984657649583;1.04770678698755;1.04558586922764;1.04348376609723;1.04140042277914;1.03933578708579;1.03728981028786;1.03526244848414;1.03325366392776;1.03126342537163;1.02929170672478;1.02733848406693;1.02540373170865;1.02348741823828;1.02158950339000;1.01970993615314;1.01784865411397;1.01600558373497;1.01418064114576;1.01237373309989;1.01058475778510;1.00881360530218;1.00706015770754;1.00532428865248;1.00360586280348;1.00190473526114;1.00022075127348;0.998553746415667;0.996903547212885;0.995269972068214;0.993652832340103;0.992051933470198;0.990467076131119;0.988898057416265;0.987344672099656;0.985806714003158;0.984283977483060;0.982776259028291;0.981283358940195;0.979805083034552;0.978341244292943;0.976891664388993;0.975456175004534;0.974034618864419;0.972626849924471;0.971232531803553;0.969847314708597;0.968446559622980;0.966961576526169;0.965273284904924;0.963238556143108;0.960737398614335;0.957714424958803;0.954194370546731;0.950268164283976;0.946059908881930;0.941690011856661;0.937246781822003;0.932772374536440;0.928263160351216;0.923680796836290;0.918969481829967;0.914075101004032;0.908962603132715;0.903628497931898;0.898106190788554;0.892463314366207;0.886791786457542;0.881192755742685;0.875759694103406;0.870563324949767;0.865641617188617;0.860996810907333;0.856599753215754;0.852400202939444;0.848340613182712;0.844370414994066;0.840457985701078;0.836598141523446;0.832813957473801;0.829152817129408;0.825677676834808;0.822455432486473;0.819544826720973;0.816986377513701;0.814796311551718;0.812965539338525;0.811463529533424;0.810245838926853;0.809263321669889;0.808470844851686;0.807833683488993;0.807330517386983;0.806952875718220;0.806701723506212;0.806589307442645;0.818727985867816;1.00588741144649;1.91128708796076;4.45410595443429;9.61193334647688;18.0006459823640;29.6332314140555;43.9622690940319;60.1137754534786;77.1605569936956;94.3277825133242;111.091698659148;127.183197017379;142.532489199134;157.193532957161;171.275966787995;184.897608484247;198.158660509707;211.132219432564;223.864044632443;236.375957796940;248.669758908944;260.730746615467;272.531165797111;284.034165391563;295.198499111895;305.983681531697;316.354941940750];
Q(1:52) = 0;
Q(362:end) = 0;
Q = [Q(53:end); Q(1:52)];

Q = Q(:)';
end

function Qtot_old = Get_parametric_old(dt, Qi, Qmax, Qmin,Qi2, ti,tmax,tmin,ti2, ts, td)


%% Q2 = fourth order (for old wave); no horizontal slope at t=0

% 15 unknowns
old_LHS_O = [
    % First polynomial curve, Q1 (systolic upslope)
    tmax^4      tmax^3      tmax^2      tmax^1  1       0           0           0       0       0       0           0           0       0       0;  %Q1(tmax)   = Qmax
    0           0           0           0       1       0           0           0       0       0       0           0           0       0       0;  %Q1(0)      = 0
    4*tmax^3    3*tmax^2    2*tmax^1    1       0       0           0           0       0       0       0           0           0       0       0;  %Q1'(tmax)  = 0
    4*3*ti^2    3*2*ti      2           0       0       0           0           0       0       0       0           0           0       0       0;  %Q1''(ti)   = 0
    ti^4        ti^3        ti^2        ti      1       0           0           0       0       0       0           0           0       0       0;  %Q1(ti)     = Qi
    % Second polynomial curve, Q2 (systolic downslope)
    0           0           0           0       0       tmax^4      tmax^3      tmax^2  tmax    1       0           0           0       0       0; %Q2(tmax)    = Qmax
    0           0           0           0       0       4*tmax^3    3*tmax^2    2*tmax  1       0       0           0           0       0       0; %Q2'(tmax)   = 0
    0           0           0           0       0       ti2^4       ti2^3       ti2^2   ti2     1       0           0           0       0       0; %Q2(ti2)     = Qi2
    0           0           0           0       0       4*3*ti2^2   3*2*ti2     2       0       0       0           0           0       0       0; %Q2'(tmax)   = 0
    
    0           0           0           0       0       ts^4        ts^3        ts^2    ts      1       -ts^4       -ts^3       -ts^2   -ts     -1; %Q2(ts)     = Q3(ts)
    0           0           0           0       0       4*ts^3      3*ts^2      2*ts    1       0       -4*ts^3     -3*ts^2     -2*ts   -1      0;  %Q2'(ts)    = Q3'(ts)  
    0           0           0           0       0       0           0           0       0       0       td^4        td^3        td^2    td      1; %Q3(td)      = 0
    0           0           0           0       0       0           0           0       0       0       4*td^3      3*td^2      2*td    1       0; %Q3'(td)     = 0
    0           0           0           0       0       0           0           0       0       0       tmin^4      tmin^3      tmin^2  tmin    1; %Q3(tmin)    = Qmin
    0           0           0           0       0       0           0           0       0       0       4*tmin^3    3*tmin^2    2*tmin  1       0  %Q3'(tmin)   = 0
    ];

%RHS
old_RHS_O = [
    Qmax;
    0;
    0;
    0;
    Qi;
    Qmax;
    0;
    Qi2;
    0;
    0;
    0;
    0;
    0;
    Qmin;
    0
];

% Pete's editing

% 15 unknowns
LHS_O = [
    % First polynomial curve, Q1 (systolic upslope) - fourth order polynomial
    tmax^4      tmax^3      tmax^2      tmax^1  1       0           0           0           0       0           0           0           0           0           0;  %Q1(tmax)   = Qmax
    0           0           0           0       1       0           0           0           0       0           0           0           0           0           0;  %Q1(0)      = 0
    4*tmax^3    3*tmax^2    2*tmax^1    1       0       0           0           0           0       0           0           0           0           0           0;  %Q1'(tmax)  = 0
    4*3*ti^2    3*2*ti      2           0       0       0           0           0           0       0           0           0           0           0           0;  %Q1''(ti)   = 0
    ti^4        ti^3        ti^2        ti      1       0           0           0           0       0           0           0           0           0           0;  %Q1(ti)     = Qi
    % Second polynomial curve, Q2 (systolic downslope) - fourth order polynomial
    0           0           0           0       0       tmax^4      tmax^3      tmax^2      tmax    1           0           0           0           0           0; %Q2(tmax)    = Qmax
    0           0           0           0       0       4*tmax^3    3*tmax^2    2*tmax      1       0           0           0           0           0           0; %Q2'(tmax)   = 0
    0           0           0           0       0       ti2^4       ti2^3       ti2^2       ti2     1           0           0           0           0           0; %Q2(ti2)     = Qi2
    0           0           0           0       0       4*3*ti2^2   3*2*ti2     2           0       0           0           0           0           0           0; %Q2''(ti2)   = 0
    0           0           0           0       0       ts^4        ts^3        ts^2        ts      1           -ts^4        -ts^3        -ts^2        -ts          -1; %Q2(ts)      = 0   % Changed: used to be Q2(ts) = Q3(ts)
    % Second and third polynomial curves
    0           0           0           0       0       5*ts^4      4*ts^3      3*ts^2      2*ts    1           -4*ts^3     -3*ts^2     -2*ts       -1          0; %Q2'(ts) - Q3'(ts) = 0
    % Third polynomial curve, Q3 (reverse flow) - fifth order polynomial
    0           0           0           0       0       0           0           0           0       0           ts^4        ts^3        ts^2        ts          1; %Q3(ts)      = 0   % Added
    0           0           0           0       0       0           0           0           0       0           td^4        td^3        td^2        td          1; %Q3(td)      = 0
    % 0           0           0           0       0       0           0           0           0       0       4*td^3      3*td^2      2*td        1           0; %Q3'(td)     = 0
    0           0           0           0       0       0           0           0           0       0           tmin^4      tmin^3      tmin^2      tmin        1; %Q3(tmin)    = Qmin
    0           0           0           0       0       0           0           0           0       0           4*tmin^3    3*tmin^2    2*tmin      1           0  %Q3'(tmin)   = 0
    ];

%RHS
RHS_O = [
    Qmax;
    0;
    0;
    0;
    Qi;
    Qmax;
    0;
    Qi2;
    0;
    0;
    0;
    0;
    0;
    Qmin;
    0
];

%Solve system A*X=B: x = A\B;
clear X;
X = old_LHS_O\old_RHS_O;


%% Equation of flow curve solution - OLD
dt=1e-3;
%Q1, t1
t1=[0:dt:tmax];
Q1 = X(1).*t1.^4 + X(2).*t1.^3 + X(3).*t1.^2 +X(4).*t1 + X(5);

%Q2,t2
t2=[tmax+dt:dt:ts];
% old
% Q2 = X(6).*t2.^4 + X(7).*t2.^3 + X(8).*t2.^2 + X(9).*t2 + X(10);
% new
Q2 = X(6).*t2.^4 + X(7).*t2.^3 + X(8).*t2.^2 + X(9).*t2 + X(10);

%Q3,t3
t3 = [ts+dt:dt:td];
% old
% Q3 = X(11).*t3.^4 + X(12).*t3.^3 + X(13).*t3.^2 + X(14).*t3 + X(15);
Q3 = X(11).*t3.^4 + X(12).*t3.^3 + X(13).*t3.^2 + X(14).*t3 + X(15);

%Q4, t4
t4 = [td+dt:dt:1.0];
Q4 = t4.*0;

Qtot_old = [Q1, Q2, Q3, Q4]; 

end

function Qtot_young = Get_parametric_young(dt, Qi, Qmax, Qmin, ti,tmax,tmin, ts, td)


%% Q2:second order, no horizontal slope = 0 at t=0 - YOUNG wave

% [P14 P13 P12 P11 P10     P22 P21 P20    P34 P33 P32 P31 P30]
LHS_Y = [
    tmax^4      tmax^3      tmax^2      tmax^1  1        0 0 0    0 0 0 0 0;  %Q1(tmax) = Qmax
    0           0           0           0       1        0 0 0    0 0 0 0 0;  %Q1(0) = 0
%    0           0           0           1       0        0 0 0    0 0 0 0 0;  %Q1'(0) = 0
    4*tmax^3    3*tmax^2    2*tmax^1    1       0        0 0 0    0 0 0 0 0;  %Q1'(tmax) = 0
    4*3*ti^2    3*2*ti      2           0       0        0 0 0    0 0 0 0 0;  %Q1''(ti) =0
    ti^4        ti^3        ti^2        ti      1        0 0 0    0 0 0 0 0;  %Q1(ti) = Qi  
    0           0           0           0       0          tmax^2   tmax 1       0 0 0 0 0; %Q2(tmax)=Qmax
    0           0           0           0       0          2*tmax   1    0       0 0 0 0 0; %Q2'(tmax)=0
    0           0           0           0       0          ts^2     ts   1       -ts^4      -ts^3     -ts^2    -ts -1; %Q2(ts)=Q3(ts)
    0           0           0           0       0           2*ts     1    0       -4*ts^3    -3*ts^2   -2*ts    -1  0;  %Q2'(ts)=Q3'(ts)  
    0 0 0 0 0    0 0 0     td^4        td^3        td^2    td      1; %Q3(td) = 0
    0 0 0 0 0    0 0 0     4*td^3      3*td^2      2*td    1       0; %Q3'(td)=0
    0 0 0 0 0    0 0 0     tmin^4      tmin^3      tmin^2  tmin    1; %Q3(tmin) = Qmin
    0 0 0 0 0    0 0 0     4*tmin^3    3*tmin^2    2*tmin  1       0  %Q3'(tmin) = 0
    ];

%RHS
RHS_Y = [
    Qmax;
    0;
%     0; %Q1'(0) = 0
    0;
    0;
    Qi;
    Qmax;
    0;
    0;
    0;
    0;
    0;
    Qmin;
    0
];

%Solve system A*X=B: x = A\B;
X = LHS_Y\RHS_Y;

%% Equation of flow curve solution - YOUNG
%Q1, t1
t1=[0:dt:tmax];
Q1 = X(1).*t1.^4 + X(2).*t1.^3 + X(3).*t1.^2 +X(4).*t1 + X(5);

%Q2,t2
t2=[tmax+dt:dt:ts];
Q2 = X(6).*t2.^2 + X(7).*t2 + X(8);

%Q3,t3
t3 = [ts+dt:dt:td];
Q3 = X(9).*t3.^4 + X(10).*t3.^3 + X(11).*t3.^2 + X(12).*t3 + X(13);

%Q4, t4
t4 = [td+dt:dt:1.0];
Q4 = t4.*0;

Qtot_young = [Q1, Q2, Q3, Q4]; 

end

function inflow = adjust_sampling_frequency(inflow)

old_t = inflow.t;
inflow.t = 0:(1/inflow.input_params.fs):inflow.t(end);
inflow.v = interp1(old_t, inflow.v, inflow.t);

end

function inflow = calc_characteristics_inflow_wave(inflow)

% NB: inflow.v is in m3 per sec

% Sampling freq (Hz)
chars.fs = inflow.input_params.fs;

% Duration of cardiac cycle (secs)
chars.T = (length(inflow.t)-1)/inflow.input_params.fs;

% Heart rate (bpm)
chars.HR = 60/chars.T;

% Stoke volume (in ml)
chars.SV = sum(inflow.v)*(1/chars.fs)*1000*1000;

% Cardiac output (in l/min)
chars.CO = chars.HR*chars.SV/1000;

% LVET (in ms)
end_systole_el = find(inflow.v(1:end-1)>0 & inflow.v(2:end) <=0, 1);
chars.LVET = inflow.t(end_systole_el);

% Forward volume (in ml)
chars.vol_forward = sum(inflow.v(inflow.v>0))*(1/chars.fs)*1000*1000;

% Regurgitation Volume
chars.Reg_Vol = abs(sum(inflow.v(inflow.v<0)))*(1/chars.fs)*1000*1000;

% Proportion regurgitation (%)
chars.perc_reg = 100*chars.Reg_Vol/chars.vol_forward;

% Measure of contractility
[chars.max_dq_dt, max_contr_el] = max(diff(inflow.v(1:end_systole_el)));
chars.max_dq_dt = chars.max_dq_dt/(1/inflow.input_params.fs);

% Maximum flow rate
chars.Qmax = max(inflow.v);

% Additional timings:
%  - Time of systolic peak
[~, sys_peak_el] = max(inflow.v);
chars.T_Peak_Flow = inflow.t(sys_peak_el);
%  - Time of max regurgitation
[max_reg, max_reg_el] = min(inflow.v);
if max_reg < 0
    chars.T_Min_Flow = inflow.t(max_reg_el);
else
    chars.T_Min_Flow = nan;
end
%  - Time of end of regurgitation
if max_reg < 0
    end_reg_el = find(inflow.v(1:end-1)<0 & inflow.v(2:end) ==0, 1);
    chars.T_dia = inflow.t(end_reg_el);
else
    chars.T_dia = nan;
end
%  - Time of max gradient of systolic upslope
chars.T_max_contr = inflow.t(max_contr_el);

inflow.chars = orderfields(chars);

end

function plot_figures(inflow)

% plot figure if requested
if inflow.input_settings.do_plot
    
    paper_size = [900 600];
        
    % Final inflow waveform in ml/s
    fig_handles = findobj('Type', 'figure');
    if ~inflow.input_settings.plot_multiple || isempty(fig_handles)
        figure('Position',[600 400 paper_size])
    end
    clear figHandles
    set(gca,'fontsize',28)
    hold on
    plot([inflow.t(1), inflow.t(end)], [0,0], '--', 'Color', 0.2*ones(1,3));
    plot(inflow.t,inflow.v*1e6,'k','linewidth',2);
    plot(inflow.t(end),inflow.v(end)*1e6,'ok','linewidth',2);
    ylabel('Q (ml/s)')
    xlabel('time (s)')
    
    %- adjust color if multiple plots
    if inflow.input_settings.plot_multiple
        h_lines = findall(gcf, 'Type', 'line');
        color_range = [0,0.75];
        if length(h_lines) > 1
            rel_colors = linspace(color_range(1), color_range(2), length(h_lines));
            for s = 1 : length(h_lines)
                h_lines(s).Color = rel_colors(s)*[1,1,1];
            end
        end      
    end
    
    if inflow.input_settings.save_plot
        save_plot(gcf, paper_size, inflow.input_settings.plot_name, inflow.input_settings.file_path)
    end
    
end

end

function save_plot(h_fig, paper_size, filename, filepath)

savepath = [filepath, filename];
set(gcf,'color','w');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize', [paper_size(1), paper_size(2)]./40);
set(gcf,'PaperPosition',[0 0 paper_size(1) paper_size(2)]./40);
print(gcf,'-dpdf',savepath)
print(gcf,'-dpng',savepath)

end
