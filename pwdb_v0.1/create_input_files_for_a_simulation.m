function FileName = create_input_files_for_a_simulation(sim_settings, up)
% CREATE_INPUT_FILES_FOR_A_SIMULATION creates a set of Nektar1D input files
% for this particular virtual subject.
%
%               create_input_files_for_a_simulation
%
%   Inputs:          (all passed from create_pwdb_input_files.m)
%               - sim_settings: a structure of settings for this particular
%                               virtual subject.
%               - up: a structure of parameters (mainly folder paths) used
%                               when creating the input files.
%
%   Outputs:    - input files for Nektar simulations:
%                           - filename.in: the main input file (containing details of vascular properties)
%                           - filename_IN_1.bcs: an additional input file (containing the input flow waveform)
%                           - filename.mat: a Matlab file containing details of the settings used to create these input files
%               - shell script files to make it easier to run a batch of simulations
%
%   Usage:      This script is called by create_pwdb_input_files.m
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
%  Copyright (C) 2019  King's College London
% 
% v.1.0  Contributed to by: Marie Willemet, Jordi Alastruey, and Peter H. Charlton
 
%% Setup parameters

% setup paths
up = setup_up(sim_settings, up);

% Setup Geometry, including
% - lengths
% - diameters
params.geo = setup_geo_params(sim_settings, up);

% Setup Physical params, including
% - Eh
% - peripheral resistance
params.phys = setup_phys_params(sim_settings, up);

% setup Nektar mesh
% (influenced by geometry)
params.mesh = setup_mesh(params, sim_settings, up);

% Store Inflow Waveform
params.inflow = sim_settings.inflow;

% Calculate Outflow
% (influenced by all sorts: PVR, lengths ...)
params.outlet = calc_outflow(params, sim_settings, up);

%% Write Nektar input files
up = write_input_files(params, up, sim_settings);
create_inflow_file(params, up);

%% Save workspace
sim_up = up; clear up
save([sim_up.paths.Matlab_computer_PathRoot, sim_up.paths.OutputFile_workspace]);

%% output filename
FileName = sim_up.paths.OutputFile(1:end-3);

end

function up = setup_up(sim_settings, up)

% prepare filename of output files

% output files (.in and .bcs)
up.paths.OutputFile = [sim_settings.filename, '.in'];                     % Generated input file for Nektar
up.paths.OutputFile_inflow = [sim_settings.filename, '_IN_1.bcs'];        % Generated inflow BC file for Nektar

% other output files
up.paths.OutputFile_workspace = [sim_settings.filename, '.mat'];          % Saved workspace variables

end

function geo_params = setup_geo_params(sim_settings, up)

%% Setup geometry

% - Extracting parameters from raw data
geo_params.node1 = sim_settings.network_spec.inlet_node;        %Nodes at proximal ends - start at 1 to 56
geo_params.node2 = sim_settings.network_spec.outlet_node;       %Nodes at distal ends
geo_params.L = sim_settings.network_spec.length;                %Segment length [m]
geo_params.Rad_prox = sim_settings.network_spec.inlet_radius;   %Proximal radius [m]
geo_params.Rad_dist = sim_settings.network_spec.outlet_radius;  %Distal radius [m]
geo_params.nbSeg = length(sim_settings.network_spec.length);    %Number of arterial segments

% - Specifying network connections
geo_params.GEO.E = geo_params.nbSeg;                            %Number of edges
geo_params.GEO.N = max([max(geo_params.node1) max(geo_params.node2)]);    %Number of nodes
geo_params.GEO.NodesAtEdges = [geo_params.node1 geo_params.node2];
geo_params.GEO = DefineGeometry(geo_params.GEO);

end

function phys_params = setup_phys_params(sim_settings, up)

%% ==== Arterial wall properties

% Empirical parameters for PWV from Mynard et al. 2015
phys_params.wall.Eh_k1 = sim_settings.network_spec.k(1); %g/s^2/cm
phys_params.wall.Eh_k2 = sim_settings.network_spec.k(2); %cm^(-1)
phys_params.wall.Eh_k3 = sim_settings.network_spec.k(3); %g/s^2/cm

% Empirical parameters for Gamma from Mynard et al. 2015
phys_params.wall.Gamma_b0 = sim_settings.gamma_b0; % baseline: 400 g/s
phys_params.wall.Gamma_b1 = sim_settings.gamma_b1; % baseline: 100 g cm/s=

%% ==== Fluid properties

phys_params.fluid.mu = sim_settings.mu; % 3.5e-3;        %Fluid Viscosity [Pa*s]
phys_params.fluid.density = sim_settings.rho; % usually 1060;     %Fluid Density [kg/m3]
phys_params.fluid.alpha = sim_settings.alpha; % 1.1;        %Velocity profile coef

end

function mesh_params = setup_mesh(params, sim_settings, up)

mesh_params.p_orderV = sim_settings.fixed.p_order.*ones(params.geo.nbSeg,1); %Vector with polynomial orders
mesh_params.q_orderV = sim_settings.fixed.q_order.*ones(params.geo.nbSeg,1); %Vector with quadrature orders
mesh_params.nbElemBySeg = ceil(params.geo.L./sim_settings.fixed.elmt); %number of elements in each segment

% If the length of the element is < 1cm, 
% - group the two last elements of the domain (if there are more than 1 element in the domain)
% - or reduce the polynomial and quadrature order to 2 (if there is only one element in the domain)
displayed_red_order = 0;
% cycle through each element
for i=1:params.geo.nbSeg
    
    % identify distances of each element. NB: this elaborate scheme avoids the
    % rounding errors sometimes produced by:      temp = [[0:sim_settings.fixed.elmt:params.geo.L(i)],params.geo.L(i)];
    no_els = ceil(params.geo.L(i)/sim_settings.fixed.elmt);
    no_complete_els = floor(params.geo.L(i)/sim_settings.fixed.elmt);
    temp = 0;
    for el_no = 1 : no_complete_els
        temp = [temp, el_no*sim_settings.fixed.elmt];
    end
    if no_els>no_complete_els
        temp(end+1) = params.geo.L(i);
    end
    
    %ensure you don't repeat the last point twice (if L=multiple of sim_settings.fixed.elmt)
    if (temp(end) == temp(end-1))
        temp=temp(1:end-1);
    end
    
    % calculate shortened length of the last element of this segment
    shortL = temp(end)-temp(end-1);
    
    % if the length of this last element is less than a threshold value, then
    if (shortL <sim_settings.fixed.elmt_shortL_reduced)
        % if there is more than one element in this segment
        if(not(mesh_params.nbElemBySeg(i) == 1))
            % then join the second last and last segment into one
            temp = [temp(1:end-2), temp(end)];
            mesh_params.nbElemBySeg(i) = mesh_params.nbElemBySeg(i)-1;
        else
            % use reduced order for this segment
            mesh_params.p_orderV(i) = sim_settings.fixed.p_reduced;
            mesh_params.q_orderV(i) = sim_settings.fixed.q_reduced;
        end
    end
    clear shortL
    
    % store the distances the elements in this segment
    mesh_params.x_elem{i} = temp; clear temp
end
clear i

% Radius(x) = (Rad_dist - Rad_prox)*x/L + Rad_prox;
%           = params.mesh.Rad_coef(:,1)*x + params.mesh.Rad_coef(:,2);
% Area(x) = pi*Radius(x)*Radius(x) 
mesh_params.Rad_coef(:,1) = (params.geo.Rad_dist - params.geo.Rad_prox)./params.geo.L;
mesh_params.Rad_coef(:,2) = params.geo.Rad_prox;

end

function create_inflow_file(params, up)

%% Save into frequency file
OutputFile = [up.paths.Matlab_computer_PathRoot, up.paths.OutputFile_inflow];
inflow.t = [0:length(params.inflow.v)-1]/params.inflow.chars.fs;
Harmonics_filter(inflow.t, params.inflow.v, length(params.inflow.v)-1, 0, OutputFile);

end

function outlet = calc_outflow(params, sim_settings, up)

%% Outlet Baseline Constants
outlet.target_PVR = sim_settings.pvr;        %Target peripheral vascular resistance [Pa s/m^3]
outlet.p_WK_outlet = sim_settings.p_out;                  %P outlet windkessel [Pa]
outlet.Pdiastolic = sim_settings.dbp;                  %Diastolic pressure [Pa]
%outlet.Pdrop = sim_settings.fixed.NEKTAR.drop_p;                     %Pressure drop from root to periphery [Pa]
outlet.R1_FIXED = sim_settings.fixed.R1_FIXED;                            %=1 if R1 is fixed via input dataFile
                                                %=0 if R1 = characteristic impedance Z0

%% Calculate properties of each arterial segment
for d=1:params.geo.nbSeg
    outlet.c_prox(d) = sqrt(2/3/params.phys.fluid.density*(0.1*(params.phys.wall.Eh_k1*exp(params.phys.wall.Eh_k2*100*params.geo.Rad_prox(d))+params.phys.wall.Eh_k3)));
    outlet.c_dist(d) = sqrt(2/3/params.phys.fluid.density*(0.1*(params.phys.wall.Eh_k1*exp(params.phys.wall.Eh_k2*100*params.geo.Rad_dist(d))+params.phys.wall.Eh_k3)));
    outlet.c_avg(d) = 0.5*(outlet.c_prox(d)+outlet.c_dist(d));
    outlet.A_avg(d) = 0.5*pi*(params.geo.Rad_prox(d)^2+params.geo.Rad_dist(d)^2);
    outlet.C_Art(d) = outlet.A_avg(d)*params.geo.L(d)/params.phys.fluid.density/outlet.c_avg(d)^2;
    outlet.Z_WK(d) = params.phys.fluid.density*outlet.c_dist(d)/pi/params.geo.Rad_dist(d)^2;
end
clear d                                            

%% Define terminal vessels belonging to each vascular bed
% This is helpful because the resistances and compliances of each vascular
% bed are known, but not those of each individual artery in leading into
% that vascular bed
outlet.VB.Domains{1} = [62 63 69 70 72 73 74 76 93 95];     outlet.VB.Name{1} = 'Brain';
outlet.VB.Domains{2} = [78 80 82 84 86 88 89 90 91 92];     outlet.VB.Name{2} = 'FaceNeck';
outlet.VB.Domains{3} = [10 102 103 104 105];                outlet.VB.Name{3} = 'RShoulderRArm';
outlet.VB.Domains{4} = [24 112 113 114 115];                outlet.VB.Name{4} = 'LShoulderLArm';
outlet.VB.Domains{5} = [26];                                outlet.VB.Name{5} = 'Chest';
outlet.VB.Domains{6} = [36];                                outlet.VB.Name{6} = 'RKidney';
outlet.VB.Domains{7} = [38];                                outlet.VB.Name{7} = 'LKidney';
outlet.VB.Domains{8} = [32];                                outlet.VB.Name{8} = 'Liver';
outlet.VB.Domains{9} = [33];                                outlet.VB.Name{9} = 'Spleen';
outlet.VB.Domains{10} = [31];                               outlet.VB.Name{10} = 'Stomach';
outlet.VB.Domains{11} = [34 40];                            outlet.VB.Name{11} = 'Intestines';
outlet.VB.Domains{12} = [45 51];                            outlet.VB.Name{12} = 'Pelvis';
outlet.VB.Domains{13} = [53];                               outlet.VB.Name{13} = 'RUppLeg';
outlet.VB.Domains{14} = [47];                               outlet.VB.Name{14} = 'LUppLeg';
outlet.VB.Domains{15} = [54 55];                            outlet.VB.Name{15} = 'RLowLeg';
outlet.VB.Domains{16} = [48 49];                            outlet.VB.Name{16} = 'LLowLeg';

% load R and C data for each vascular bed from the literature
WK_data = sim_settings.wk_params;
outlet.VB.R = WK_data(:,1)*133.33*1e6; %Vascular bed resistance [Pa s/m^3]
outlet.VB.C = sim_settings.pvc*WK_data(:,2)/133.33/1e6; %Vascular bed compliance [m^3/Pa]
clear WK_data

%% OUTLET: RCR Windkessel parameters
% distribute Z, R and C between individual arteries leading into each vascular bed

outlet.Rp_tot = 0;  % initialise total peripheral resistance across entire model

% cycle through vascular beds
for j=1:length(outlet.VB.Domains)
    
    %- Find total area across this vascular bed
    A_tot = 0;  % initialise total area
    for d=1:length(outlet.VB.Domains{j})
        A_tot = A_tot + pi*params.geo.Rad_dist(outlet.VB.Domains{j}(d))^2;
    end
    clear d
    
    %- Calculate Z, R and C for each artery leading into this vascular bed
    for d=1:length(outlet.VB.Domains{j})
        
        % R outflow = share of vascular bed's R
        outlet.R_WK_outflow(outlet.VB.Domains{j}(d)) = outlet.VB.R(j)*A_tot/pi/params.geo.Rad_dist(outlet.VB.Domains{j}(d))^2;
        % R total = R outflow + Z
        outlet.R_WK_tot(outlet.VB.Domains{j}(d)) = outlet.Z_WK(outlet.VB.Domains{j}(d)) + outlet.R_WK_outflow(outlet.VB.Domains{j}(d));
        % C = share of vascular bed's C
        outlet.C_WK(outlet.VB.Domains{j}(d)) = outlet.VB.C(j)*pi*params.geo.Rad_dist(outlet.VB.Domains{j}(d))^2/A_tot;
        % C_WK_IW = C * R outflow / R total
        outlet.C_WK_IW(outlet.VB.Domains{j}(d)) = outlet.C_WK(outlet.VB.Domains{j}(d))*outlet.R_WK_outflow(outlet.VB.Domains{j}(d))/ outlet.R_WK_tot(outlet.VB.Domains{j}(d));
        
        % update total peripheral resistance across entire model
        outlet.Rp_tot = outlet.Rp_tot + 1/outlet.R_WK_tot(outlet.VB.Domains{j}(d));
        
    end
    clear d
    
end
clear j

% Adjust terminal resistances to obtain the target PVR (accounting for
% pressure drop within the 1-D model network)
outlet.Rp_tot = 1/outlet.Rp_tot; %Peripheral vascular resistance (Pa s/m^3)
%outlet.Rp_tot_drop = outlet.Rp_tot + outlet.Pdrop/params.inflow.CO_m3_per_sec; %Net peripheral vascular resistance (Pa s/m^3)
CO_m3_per_sec = params.inflow.chars.CO/(60*1000);  % converting from l/min to m3 / sec
%scale_WK = (outlet.target_PVR - outlet.Pdrop/CO_m3_per_sec)/outlet.Rp_tot; %Calculate the scaling factor for R_WK_tot
scale_WK = outlet.target_PVR/outlet.Rp_tot; %Calculate the scaling factor for R_WK_tot
outlet.R_WK_tot = outlet.R_WK_tot * scale_WK; %Scale WK resitances to obtain the target PVR
outlet.Rp_tot = outlet.Rp_tot * scale_WK; %Scaled peripheral vascular resistance (Pa s/m^3)
%Rp_tot/133.33/1e6; %(mmHg s/ml)

% Calculate total compliances
outlet.Cp_tot = sum(outlet.VB.C); %Total peripheral compliance (m^3/Pa)
outlet.Cp_IW_tot = sum(outlet.C_WK_IW); %Total impedance-weighted peripheral compliance (m^3/Pa)

outlet.C_Art_tot = sum(outlet.C_Art); %Total arterial compliance (m^3/Pa)
outlet.C_tot = outlet.C_Art_tot+outlet.Cp_IW_tot; %Total compliance (m^3/Pa)
outlet.C_tot*133.33*1e6; %(ml/mmHg)

end

function up = write_input_files(params, up, sim_settings)


%% Open files
fileIn = fopen([up.paths.Matlab_computer_PathRoot, up.paths.OutputFile],'w');

%% Determine time step for simulation
sim_settings.dt = sim_settings.time_step;

%% Constant parameters and variables

% === Arrange

nbP = 0;   %number of parameters
nbP=nbP+1; varType{nbP} = '%d';    varName{nbP} = 'EQTYPE';     varVal(nbP) = sim_settings.fixed.LinearForm;
nbP=nbP+1; varType{nbP} = '%1.3e'; varName{nbP} = 'DT';         varVal(nbP) = sim_settings.dt;
nbP=nbP+1; varType{nbP} = '%1.3e'; varName{nbP} = 'NSTEPS';     varVal(nbP) = sim_settings.fixed.duration/sim_settings.dt;
nbP=nbP+1; varType{nbP} = '%d';    varName{nbP} = 'HISSTEP';    varVal(nbP) = round(sim_settings.fixed.output_dt/sim_settings.dt);
nbP=nbP+1; varType{nbP} = '%d';    varName{nbP} = 'INTTYPE';    varVal(nbP) = sim_settings.fixed.int_time;
nbP=nbP+1; varType{nbP} = '%3.1f'; varName{nbP} = 'Rho';        varVal(nbP) = params.phys.fluid.density;
nbP=nbP+1; varType{nbP} = '%1.4f'; varName{nbP} = 'Viscosity';  varVal(nbP) = params.phys.fluid.mu;
nbP=nbP+1; varType{nbP} = '%1.3f'; varName{nbP} = 'Alpha';      varVal(nbP) = params.phys.fluid.alpha;
nbP=nbP+1; varType{nbP} = '%2.3f'; varName{nbP} = 'T_initial';  varVal(nbP) = sim_settings.fixed.duration - sim_settings.inflow.chars.T;
nbP=nbP+1; varType{nbP} = '%2.3f'; varName{nbP} = 'T_final';    varVal(nbP) = sim_settings.fixed.duration;
nbP=nbP+1; varType{nbP} = '%3.3f'; varName{nbP} = 'pinf';       varVal(nbP) = params.outlet.p_WK_outlet;
nbP=nbP+1; varType{nbP} = '%2.3e'; varName{nbP} = 'Pext';       varVal(nbP) = params.outlet.Pdiastolic;

% === Print to files

% number of constant parameters and variables (no 'pinf' and 'Pext' for conduit models)
if sim_settings.reflect_log % a value of 1 indicates 'reflect'
    no_vars = nbP;
    outlet_wk_type = 'W';   % used when specifying the outlet boundary conditions
elseif ~sim_settings.reflect_log  % a value of 0 indicates 'absorb'
    no_vars = nbP-2;
    no_vars = nbP;  % trial by PC, wasn't here before
    outlet_wk_type = 'R';
else
    error('Unrecognised terminal type')
end

% header line (stating how many constant parameters and variables there are)
fprintf(fileIn,'%d \t %s \n',no_vars,' parameter list');
% print each constant parameters and variable
for k=1:no_vars
    fprintf(fileIn,[varType{k},' \t %s\n'],varVal(k),varName{k});
end
clear k

%% Mesh
% Provide values of area, beta and gamma for each arterial segment (each domain)

% header line
fprintf(fileIn,'%s %d \n','Mesh -- expansion order --  quadrature order Ndomains =',params.geo.nbSeg);

% cycle through each segment (domain)
for e=1:params.geo.nbSeg
    
    % Make string for radius
    Radius = ['(',num2str(params.mesh.Rad_coef(e,1)),'*x+',num2str(params.mesh.Rad_coef(e,2)),')'];
        
    % Make string for Eh
    Eh = [num2str(0.1),'*(',num2str(params.phys.wall.Eh_k1),'*exp(',num2str(params.phys.wall.Eh_k2),'*100*',Radius,') +',num2str(params.phys.wall.Eh_k3),')*',Radius];
    
    % Make string for gamma
    Gamma = ['4/PI*(',num2str(params.phys.wall.Gamma_b1),'/(2*100*',Radius,')+',num2str(params.phys.wall.Gamma_b0),')/1000/(2*',Radius,')^2'];
    
    % Trial version: ghigo_gamma = (2/3)*sqrt(pi)*0.017.*Eh./(pi*(r.^2));
    
    %Gamma = ['(2/3)*sqrt(PI)*0.017*' Eh '/(PI*(' Radius '*' Radius '))'];
    
    %    Beta = [num2str(0.1*4/3/sqrt(pi)),'*(',num2str(params.phys.wall.Eh_k1),'*exp(',num2str(params.phys.wall.Eh_k2),'*0.01*',Radius,') +',num2str(params.phys.wall.Eh_k3),')/',Radius];
    
    % Print header line stating: Eh, Area (and gamma if visco-elastic)
    if (params.geo.L(e)>sim_settings.fixed.elmt_shortL_vw) && sim_settings.visco_elastic_log
        fprintf(fileIn,'%d %s %d %s\n',params.mesh.nbElemBySeg(e),'    nel domain ',e, 'Eh Area Gamma');
    elseif (params.geo.L(e)<=sim_settings.fixed.elmt_shortL_vw) || ~sim_settings.visco_elastic_log
        fprintf(fileIn,'%d %s %d %s\n',params.mesh.nbElemBySeg(e),'    nel domain ',e, 'Eh Area');
    else
        error('\n Unrecognised wall type')
    end
    
    % Print individual values for each element in this segment
    x_temp = params.mesh.x_elem{e};
    % cycle through elements
    for j=1:params.mesh.nbElemBySeg(e)
        
        % identify distances (x values)
        l_start = x_temp(j);
        l_end   = x_temp(j+1);
        
        % print mesh settings
        fprintf(fileIn,'%3.4f %3.4f %d %d %s\n', l_start, l_end, params.mesh.p_orderV(e), params.mesh.q_orderV(e), '# x_prox x_dist p q');
        
        % print area
        fprintf(fileIn,'%s %s %s %s\n','Area = PI*', Radius,'*',Radius);
        
        % print Eh
        fprintf(fileIn,'%s %s\n','Eh = ', Eh); %If prescribed with x-varying Area and Eh
        
        % print gamma if this segment is visco-elastic
        if (params.geo.L(e)>sim_settings.fixed.elmt_shortL_vw) && sim_settings.visco_elastic_log
            fprintf(fileIn,'%s %s\n','Gamma = ', Gamma); %Visco-elastic coefficient
        end
    end
end

%% Boundary Conditions

% header
fprintf(fileIn,'%s\n','Boundary conditions');

% cycle through each segment (domain)
for i=1:params.geo.nbSeg

    %--- Inlet Boundary Conditions

    if i==1
        %- Domain 1 : Specify inlet as reflective flow
        fprintf(fileIn,'%s \t%s\n','F  0',' # Domain 1'); %Frequency input (sum of harmonics in bcs file (F 0 or F 3))
        fprintf(fileIn,'%s\n','F  0');
        
    else
        %- all other domains: specify which segments they are connected to
        
        % properties of inlet
        inNode = params.geo.node1(i);
        deg = params.geo.GEO.degree(inNode);
        inConnect = params.geo.GEO.EdgesByNode(inNode,:);
        
        % switch according to the type of junction at the inlet
        switch deg
            
            case 1 %no connectivity at inlet
                disp(['there is a problem at segment ',num2str(i)]);
                return;
                
            case 2 %junction
                ij = find(not(inConnect([1 2]) == i));
                fprintf(fileIn,'%s %d %d  \t %s%d\n','J',inConnect(ij),inConnect(ij), ' # Domain ',i);
                fprintf(fileIn,'%s %d %d\n','J',inConnect(ij),inConnect(ij));
                
            case 3 %bifurcation
                ib = find(not(inConnect == i));  % the segments which join to the current segment
                % determine type of bifurcation
                if (params.geo.GEO.BifType(inNode)<0)
                    BifType = 'B';
                else
                    BifType = 'C';
                end
                fprintf(fileIn,'%s %d %d  \t %s%d\n',BifType,inConnect(ib(1)),inConnect(ib(2)), ' # Domain ',i);
                fprintf(fileIn,'%s %d %d\n',BifType,inConnect(ib(1)),inConnect(ib(2)));
        end
    end
    
    %--- Outlet Boundary Conditions
    
    % properties of outlet
    outNode = params.geo.node2(i);
    deg = params.geo.GEO.degree(outNode);
    outConnect = params.geo.GEO.EdgesByNode(outNode,:);
    
    % switch according to the type of junction at the outlet
    switch deg
        
        case 1 %outlet WK
            % varies according to whether the terminal boundary conditions absorb reflections (purely reistance), or reflect (windkessel)
            
            % imaginary part
            if strcmp(outlet_wk_type, 'W')
                fprintf(fileIn,'%s %1.4e\n', outlet_wk_type, params.outlet.C_WK(i)); %WK Compliance
            else
                fprintf(fileIn,'%s %1.4e\n',outlet_wk_type, params.outlet.Z_WK(i)); % R = Z0
            end
            
            % real part
            if (params.outlet.R1_FIXED == 1)
                fprintf(fileIn,'%s %1.4e %1.4f\n','W',params.outlet.R_WK_tot(i), params.outlet.r_WK(i)); %WK Resistance - R1 value fixed
            else
                fprintf(fileIn,'%s %1.4e\n','W',params.outlet.R_WK_tot(i)); %WK Resistance - R1 = Z0
            end
            
            % resistance for absorbing terminal
            if strcmp(outlet_wk_type, 'R')
                fprintf(fileIn,'%s %1.4e\n','R',params.outlet.Z_WK(i)); % R = Z0
            end
            
        case 2 %junction
            ij = find(not(outConnect([1 2]) == i));
            fprintf(fileIn,'%s %d %d\n','J',outConnect(ij),outConnect(ij));
            fprintf(fileIn,'%s %d %d\n','J',outConnect(ij),outConnect(ij));
            
        case 3 %bifurcation
            ib = find(not(outConnect == i));  % the segments which join to the current segment
            % determine type of bifurcation
            if (params.geo.GEO.BifType(outNode)<0)
                BifType = 'B';
            else
                BifType = 'C';
            end
            fprintf(fileIn,'%s %d %d\n',BifType,outConnect(ib(1)),outConnect(ib(2)));
            fprintf(fileIn,'%s %d %d\n',BifType,outConnect(ib(1)),outConnect(ib(2)));
    end
end

%% Initial conditions

% header
fprintf(fileIn,'%s\n','Initial condition');

% cycle through each segment (domain)
for e=1:params.geo.nbSeg
    
    % Make string for radius
    Radius = ['(',num2str(params.mesh.Rad_coef(e,1)),'*x+',num2str(params.mesh.Rad_coef(e,2)),')'];
    
    % Make string for Eh
    Eh = [num2str(0.1),'*(',num2str(params.phys.wall.Eh_k1),'*exp(',num2str(params.phys.wall.Eh_k2),'*100*',Radius,') +',num2str(params.phys.wall.Eh_k3),')*',Radius];
    
    % print initial area
    %if strcmp(outlet_wk_type, 'W')
        fprintf(fileIn,'%s %s %s %s %s %s %s %s %s %s %s\n','a = PI*',Radius,'*',Radius,'*(1-3/4*',Radius,'*',num2str(params.outlet.Pdiastolic),'/(',Eh,'))^2');
    %else
    %    fprintf(fileIn_cond,'%s\n','a = Ao');
    %end
    
    % print initial velocity
    fprintf(fileIn,'%s\n','u = 0');
    
end

%% History points

% header
fprintf(fileIn,'%s\n','History Pts');
fprintf(fileIn,'%d %s\n',length(sim_settings.fixed.artery_output),' #Number of Domains with history Points');

%- Outputs at inlet, middle and outlet of some segment
locs = cell(length(sim_settings.fixed.artery_output),1);
for seg_no_el = 1 : length(sim_settings.fixed.artery_output)
    curr_seg_no = sim_settings.fixed.artery_output(seg_no_el);
    temp = [0, params.geo.L(curr_seg_no)./2, params.geo.L(curr_seg_no)];
    % Add additional points every 2 cm along arterial segments in path from asc aorta to digital
    if sum(sim_settings.fixed.Path_dig == curr_seg_no) ~= 0 ...
            || sum(sim_settings.fixed.Path_ankle == curr_seg_no) ~= 0 ...
        || sum(sim_settings.fixed.Path_brain == curr_seg_no) ~= 0
        temp = [temp, 0.04:0.04:params.geo.L(curr_seg_no)];
    end
    % 3/4 of way along brachial arteries
    if curr_seg_no == 7 || curr_seg_no == 21
        temp = [temp, 0.75*params.geo.L(curr_seg_no)];
    end
    locs{seg_no_el} = unique(temp); clear temp
end

% cycle through each segment at which to write output measurements
for seg_no_el=1:length(sim_settings.fixed.artery_output)
    
    % print to file
    fprintf(fileIn,'%d %d\n',length(locs{seg_no_el}), sim_settings.fixed.artery_output(seg_no_el));
    format_spec = [repmat('%1.4f ', [1,length(locs{seg_no_el})-1]), '%1.4f\n'];
    fprintf(fileIn,format_spec,locs{seg_no_el});
    
end

fprintf(fileIn,'\n \n \n %s %s\n','oneDbio -a -R', up.paths.OutputFile);

%% Close file
fclose(fileIn);

end

function [GEO] = DefineGeometry(GEO)
% GEO contains the variables
%   GEO.N : number of nodes in network
%   GEO.E : number of edges in network
%   GEO.NodesAtEdges (E x 2): Start and end nodes for each edge

% The function constructs the following geometry matrices:
%   GEO.B            (NxE)  : Incidence matrix
%   GEO.EdgesByNode  (N x 3): All edges connected to each node
%   GEO.degree       (1xN)  : Degree of each node
%   GEO.BifType      (1xN)  : Type of bifurcation for each node

%--- Incidence Matrix B
GEO.B = zeros(GEO.N, GEO.E);
GEO.BifType = zeros(GEO.N, 1);

for (i = 1:GEO.E)
    %start of node
    GEO.B(GEO.NodesAtEdges(i,1),i) = -1;
    %end of node
    GEO.B(GEO.NodesAtEdges(i,2),i) = 1;
end
    

%-- EdgesByNodes & degree of nodes
GEO.EdgesByNode(1:GEO.N,1:3) = 0;

for i=1:GEO.N
    % Get the degree of node i
    GEO.degree(i) = sum(abs(GEO.B(i,:)));
    if (GEO.degree(i)==3) 
        GEO.BifType(i) = sum(GEO.B(i,:));
    end
    % Get the edges connected to the internal node i
    clear e Y;
    s=1; 
    for j=1:GEO.E  
        if abs(GEO.B(i,j))>0
            e(s)=j;
            s=s+1;
        end
    end
    GEO.EdgesByNode(i,:) = [e zeros(1,4-s)];
end

end