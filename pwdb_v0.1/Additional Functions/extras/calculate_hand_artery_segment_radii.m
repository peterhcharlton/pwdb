function calculate_hand_artery_segment_radii

up.matched_ratio = 1.15; % as stated in the introduction of: Greenwald SE, Newman DL. Impulse propagation through junctions. Med Biol Eng Comput 20: 343?350, 1982.

% %% Load data
% 
% % Give details of the Excel file which contains the baseline model geometry:
% up.network_spec_sheet = 'Haemod 116 segments';
% curr_filepath = fileparts(mfilename('fullpath'));
% up.network_spec_file = [curr_filepath, '/input_data/116_artery_model_v10_before_hand_adj.xlsx'];
% [num,txt,raw] = xlsread(up.network_spec_file, up.network_spec_sheet);
% rel_col = strcmp(raw(1,:), 'Segment No');
% network_spec.seg_no = num(:, rel_col);
% rel_col = strcmp(raw(1,:), 'Inlet node');
% network_spec.inlet_node = num(:, rel_col);
% rel_col = strcmp(raw(1,:), 'Outlet node');
% network_spec.outlet_node = num(:, rel_col);
% rel_col = strcmp(raw(1,:), 'Length [m]');
% network_spec.length = num(:, rel_col);
% rel_col = strcmp(raw(1,:), 'Inlet Radius [m]');
% network_spec.inlet_radius = num(:, rel_col);
% rel_col = strcmp(raw(1,:), 'Outlet Radius [m]');
% network_spec.outlet_radius = num(:, rel_col);
% rel_col = strcmp(raw(1,:), 'Mynard Name');
% network_spec.segment_name = raw(2:size(num,1)+1, rel_col);

%% Data from Alastruey 2006 (model 2)
old_outlet_radius = nan(116,1);
% Left Hand Side
old_outlet_radius(22) = 0.00174;  % matches article (radial)
old_outlet_radius(25) = 0.00203;  % matches article (ulnar)
old_outlet_radius(107) = 0.00093;  % matches DPA 3 (which is the correct segment when comparing the Nektar sketch to the diagram of the model in the paper)
old_outlet_radius(108) = 0.00161;  % matches SPA 3 (which is the correct segment when comparing the Nektar sketch to the diagram of the model in the paper)
old_outlet_radius(109) = 0.00130;  % matches DPA 1 (which is the correct segment when comparing the Nektar sketch to the diagram of the model in the paper)
old_outlet_radius(110) = 0.00169;  % matches SPA 1 (which is the correct segment when comparing the Nektar sketch to the diagram of the model in the paper)
old_outlet_radius(111) = 0.00114;  % matches SPA II (which is the correct segment when comparing the Nektar sketch to the diagram of the model in the paper)
old_outlet_radius(112) = 0.00132;  % matches DIG 3 (which is the correct segment when comparing the Nektar sketch to the diagram of the model in the paper)
old_outlet_radius(113) = 0.00136;  % matches DIG 2 (which is the correct segment when comparing the Nektar sketch to the diagram of the model in the paper)
old_outlet_radius(114) = 0.00093;  % matches DIG 4 (which is the correct segment when comparing the Nektar sketch to the diagram of the model in the paper)
old_outlet_radius(115) = 0.00130;  % matches DIG 1 (which is the correct segment when comparing the Nektar sketch to the diagram of the model in the paper)
old_outlet_radius(116) = 0.00075;  % matches DPA 2 (which is the correct segment when comparing the Nektar sketch to the diagram of the model in the paper)
% Right Hand Side
old_outlet_radius(8) = 0.00174;
old_outlet_radius(11) = 0.00203;
old_outlet_radius(97) = 0.00093;
old_outlet_radius(98) = 0.00161;
old_outlet_radius(99) = 0.00130;
old_outlet_radius(100) = 0.00169;
old_outlet_radius(101) = 0.00114;
old_outlet_radius(102) = 0.00132;
old_outlet_radius(103) = 0.00136;
old_outlet_radius(104) = 0.00093;
old_outlet_radius(105) = 0.00130;
old_outlet_radius(106) = 0.00075;
old_areas = pi*(old_outlet_radius.^2);

%% New data
% These are the data for the network before scaling by sqrt(1.5)
new_areas = old_areas;
new_areas(22) = pi*(0.00107^2);
new_areas(25) = pi*(0.00183^2);
new_areas(8) = pi*(0.00107^2);
new_areas(11) = pi*(0.00183^2);

%% Find areas

% modify areas for hand segments
new_areas = modify_bifurcations(new_areas, 22, 107, 108, up);   % junction at end of left radial
new_areas = modify_bifurcations(new_areas, 25, 109, 110, up);   % junction at end of left ulnar
new_areas = modify_bifurcations(new_areas, 8, 97, 98, up);      % junction at end of right radial
new_areas = modify_bifurcations(new_areas, 11, 99, 100, up);    % junction at end of right ulnar

% modify areas which are common (anastomoses)
new_areas = modify_common_areas(old_areas, new_areas, 111, 108, 110, up);   % central section of left superficial palmar arch
new_areas = modify_common_areas(old_areas, new_areas, 116, 109, 107, up);   % left deep palmar arch
new_areas = modify_common_areas(old_areas, new_areas, 101, 98, 100, up);    % central section of right superficial palmar arch
new_areas = modify_common_areas(old_areas, new_areas, 106, 99, 97, up);     % right deep palmar arch

% modify terminal branches according to proportion of area
new_areas = modify_terminal_branches(old_areas, new_areas, 112, 108);
new_areas = modify_terminal_branches(old_areas, new_areas, 113, 110);
new_areas = modify_terminal_branches(old_areas, new_areas, 102, 98);
new_areas = modify_terminal_branches(old_areas, new_areas, 103, 100);
new_areas = modify_terminal_branches(old_areas, new_areas, 115, 109);
new_areas = modify_terminal_branches(old_areas, new_areas, 114, 107);
new_areas = modify_terminal_branches(old_areas, new_areas, 105, 99);
new_areas = modify_terminal_branches(old_areas, new_areas, 104, 97);

new_radii = sqrt(new_areas/pi);

% Output new radii
fprintf('\n\n~~~~~~~~~~ New Radii ~~~~~~~~~')
for seg_no = [8,11,22,25,97:116]
    fprintf(['\n new_radius(' num2str(seg_no) ') = ' num2str(new_radii(seg_no),3) ]);
end

fprintf('\n\n~~~~~~~~~~ New Radii, after scaling by sqrt(1.5) ~~~~~~~~~')
for seg_no = [8,11,22,25,97:116]
    fprintf(['\n new_radius(' num2str(seg_no) ') = ' num2str(new_radii(seg_no),3) ]);
end

end

function new_areas = modify_bifurcations(new_areas, inlet_seg, outlet_seg1, outlet_seg2, up)

prop1 = new_areas(outlet_seg1)/(new_areas(outlet_seg2) + new_areas(outlet_seg1));

new_total_area = new_areas(inlet_seg)*up.matched_ratio;

new_areas(outlet_seg1) = prop1*new_total_area;
new_areas(outlet_seg2) = (1-prop1)*new_total_area;

end

function new_areas = modify_common_areas(old_areas, new_areas, common_seg, inlet1, inlet2, up)

old_prop = mean([old_areas(common_seg)/old_areas(inlet1), old_areas(common_seg)/old_areas(inlet2)]);

new_areas(common_seg) = old_prop*mean([new_areas(inlet1), new_areas(inlet2)]);

end

function new_areas = modify_terminal_branches(old_areas, new_areas, terminal_seg, inlet_seg)

old_prop = old_areas(terminal_seg)/old_areas(inlet_seg);
new_areas(terminal_seg) = old_prop*new_areas(inlet_seg);

end