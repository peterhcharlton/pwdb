function calculate_aorta_leg_segment_radii

close all
up.matched_ratio = 1.15;
up.Dsys_Ddia_ratio = 1.08; % looks approximately right from an initial simulation

%% Load data

% Give details of the Excel file which contains the baseline model geometry:
up.network_spec_sheet = 'Haemod 116 segments';
curr_filepath = fileparts(mfilename('fullpath'));
up.network_spec_file = [curr_filepath, '/input_data/116_artery_model_v10_after_hand_adj_bef_aorta_adj.xlsx'];
[num,txt,raw] = xlsread(up.network_spec_file, up.network_spec_sheet);
rel_col = strcmp(raw(1,:), 'Segment No');
network_spec.seg_no = num(:, rel_col);
rel_col = strcmp(raw(1,:), 'Inlet node');
network_spec.inlet_node = num(:, rel_col);
rel_col = strcmp(raw(1,:), 'Outlet node');
network_spec.outlet_node = num(:, rel_col);
rel_col = strcmp(raw(1,:), 'Length [m]');
network_spec.length = num(:, rel_col);
rel_col = strcmp(raw(1,:), 'Inlet Radius [m]');
network_spec.inlet_radius = num(:, rel_col);
rel_col = strcmp(raw(1,:), 'Outlet Radius [m]');
network_spec.outlet_radius = num(:, rel_col);
rel_col = strcmp(raw(1,:), 'Mynard Name');
network_spec.segment_name = raw(2:size(num,1)+1, rel_col);
a = network_spec;

%% Identify segments to be changed

aorta_segs = [1 2 14 18 27 28 35 37 39 41];
leg_segs = 42:55;

%% Initial data

lengths = a.length;
old_inlet_areas = pi*(a.inlet_radius.^2);
old_outlet_areas = pi*(a.outlet_radius.^2);

%% Desired data (diameters)
desired_rads = [];

% Asc Aorta (Hickson L1)
rel_el = 1;
rel_dist_prop = 0.5;
article_data.val.age = [24,34,45,57,63,73];
article_data.val.v = 10 + (25*[35.5,38.5,41.75,42.0,44.0,44.5]./47);
article_dia_mm = interp1(article_data.val.age, article_data.val.v, 25);
desired_rad = article_dia_mm/(2000*up.Dsys_Ddia_ratio);
desired_rads = [desired_rads; rel_el, rel_dist_prop, desired_rad];

% Desc Thor Aorta (Hickson L2)
rel_el = 18;
rel_dist_prop = 0.5;
article_data.val.age = [24,34,45,57,63,73];
article_data.val.v = 10 + (25*[48,53.5,58.5,58,65,64.5]./108);
article_dia_mm = interp1(article_data.val.age, article_data.val.v, 25);
desired_rad = article_dia_mm/(2000*up.Dsys_Ddia_ratio);
desired_rads = [desired_rads; rel_el, rel_dist_prop, desired_rad];

% Diaphragm (Hickson L3)   %% CHECK
rel_el = 27;
rel_dist_prop = 0.9;
article_data.val.age = [24,34,45,57,63,73];
article_data.val.v = 10*[1.85, 1.9, 2.05, 2.05, 2.2, 2.25];
article_dia_mm = interp1(article_data.val.age, article_data.val.v, 25);
desired_rad = article_dia_mm/(2000*up.Dsys_Ddia_ratio);
desired_rads = [desired_rads; rel_el, rel_dist_prop, desired_rad];

% Abd Aorta (Hickson L4)
rel_el = 39;
rel_dist_prop = 0.1;
article_data.val.age = [24,34,45,57,63,73];
article_data.val.v = 10 + (25*[29.5,35,38,37,42.5,44]./108);
article_dia_mm = interp1(article_data.val.age, article_data.val.v, 25);
desired_rad = article_dia_mm/(2000*up.Dsys_Ddia_ratio);
desired_rads = [desired_rads; rel_el, rel_dist_prop, desired_rad];

% 3cm above Aortic Bifurcation (Hickson L5)   %% CHECK
rel_el = 39;
rel_dist_prop = 0.9;
article_data.val.age = [24,34,45,57,63,73];
article_data.val.v = 10*[1.6, 1.65, 1.65, 1.7, 1.75, 1.75];
article_dia_mm = interp1(article_data.val.age, article_data.val.v, 25);
desired_rad = article_dia_mm/(2000*up.Dsys_Ddia_ratio);
desired_rads = [desired_rads; rel_el, rel_dist_prop, desired_rad];

% Calculate lengths along aorta
for s = 1 : size(desired_rads,1)
    rel_seg_no = desired_rads(s,1);
    rel_prop = desired_rads(s,2);
    temp = find(aorta_segs == rel_seg_no);
    rel_full_els = aorta_segs(1:temp-1);
    lens(s,1) = sum(lengths(rel_full_els)) + rel_prop*lengths(rel_seg_no);
end

%% Fit curve to measured radii
rads = desired_rads(:,3);
tbl = table(lens,rads);
%mdl = fitlm(tbl,'quadratic');
f = fit(lens,rads,'cubicinterp');
%plot(f,lens,rads)

%% Plot fitted and actual data
% Calculate fitted
rel_segs = aorta_segs; 
dists = [0; cumsum(a.length(rel_segs))];
fitted = feval(f,0:0.01:dists(end));
% Plot fitted
subplot(1,2,1)
plot(0:0.01:dists(end), fitted*2000); hold on
plot(lens, rads*2000, '*k'); hold on
% Calculate actual
d1 = dists(1:end-1); r1 = a.inlet_radius(rel_segs);
d2 = dists(2:end); r2 = a.outlet_radius(rel_segs);
[~,order] = sort([d1;d2]);
for s = 1 : length(d1)
    d(2*s-1) = d1(s); d(2*s) = d2(s); r(2*s-1) = r1(s); r(2*s) = r2(s);
end
% Plot actual
subplot(1,2,1)
plot(d,2*1000*r)

%% Find old reflection ratios
for seg_no = 1 : length(a.outlet_node)
    curr_outlet_node = a.outlet_node(seg_no);
    inlet_segs = find(a.inlet_node == curr_outlet_node);
    inlet_seg_radii = a.inlet_radius(inlet_segs);
    outlet_area = pi*a.outlet_radius(seg_no)^2;
    inlet_area = sum(pi*inlet_seg_radii.^2);
    refl_ratios(seg_no,1) = outlet_area / inlet_area;
end

%% New data
new = a;
for seg_no = 1 : length(aorta_segs)
    curr_seg = aorta_segs(seg_no);
    % Remove any tapering from this segment
    new.outlet_radius(curr_seg) = new.inlet_radius(curr_seg);
    % Adjust the next segment's inlet radius to ensure that the same ratio is kept
    outlet_area = pi*new.outlet_radius(curr_seg)^2;
    curr_outlet_node = a.outlet_node(curr_seg);
    inlet_segs = find(a.inlet_node == curr_outlet_node);
    inlet_seg_radii = a.inlet_radius(inlet_segs);
    inlet_area = sum(pi*inlet_seg_radii.^2);
    curr_refl_ratio = outlet_area / inlet_area;
    desired_refl_ratio = refl_ratios(curr_seg);
    area_scaling = curr_refl_ratio/desired_refl_ratio;
    radius_scaling = sqrt(area_scaling);
    new.inlet_radius(inlet_segs) = inlet_seg_radii*radius_scaling;
    
end

%% Find new reflection ratios (just to check)
for seg_no = 1 : length(aorta_segs)
    curr_seg = aorta_segs(seg_no);
    curr_outlet_node = new.outlet_node(curr_seg);
    inlet_segs = find(new.inlet_node == curr_outlet_node);
    inlet_seg_radii = new.inlet_radius(inlet_segs);
    outlet_area = pi*new.outlet_radius(curr_seg)^2;
    inlet_area = sum(pi*inlet_seg_radii.^2);
    new_refl_ratios(seg_no,1) = outlet_area / inlet_area;
end

%% Plot new actual (aorta)

% Calculate actual
d1 = dists(1:end-1); r1 = new.inlet_radius(rel_segs);
d2 = dists(2:end); r2 = new.outlet_radius(rel_segs);
for s = 1 : length(d1)
    d(2*s-1) = d1(s); d(2*s) = d2(s); r(2*s-1) = r1(s); r(2*s) = r2(s);
end
% Plot actual
subplot(1,2,1)
plot(d,2*1000*r)

%% Find old leg tapering proportions
for s = 1 : length(leg_segs)
    
    % Inlet radius has already been increased
    
    % Find old tapering proportions
    curr_seg = leg_segs(s);
    old_outlet_area = pi*(a.outlet_radius(curr_seg)^2);
    old_inlet_area = pi*(a.inlet_radius(curr_seg)^2);
    old_prop_tapering = old_outlet_area/old_inlet_area;
    
    % Adjust outlet radius to keep old tapering proportions
    new_inlet_area = pi*(new.inlet_radius(curr_seg)^2);
    new_outlet_area = new_inlet_area*old_prop_tapering;
    new.outlet_radius(curr_seg) = sqrt(new_outlet_area/pi);
    
    % Adjust inlet radii of next segments to maintain reflection ratio
    curr_outlet_node = new.outlet_node(curr_seg);
    inlet_segs = find(new.inlet_node == curr_outlet_node);
    inlet_seg_radii = new.inlet_radius(inlet_segs);
    inlet_area = sum(pi*inlet_seg_radii.^2);
    curr_refl_ratio = new_outlet_area / inlet_area;
    desired_refl_ratio = refl_ratios(curr_seg);
    area_scaling = curr_refl_ratio/desired_refl_ratio;
    radius_scaling = sqrt(area_scaling);
    new.inlet_radius(inlet_segs) = inlet_seg_radii*radius_scaling;
end

%% Plot old actual (legs)

init_dist = dists(end);
rel_segs = [42, 44, 46, 48];
rel_segs = [43, 50, 52, 55];
dists = cumsum([init_dist; a.length(rel_segs)]);

% Calculate actual
d1 = dists(1:end-1); r1 = a.inlet_radius(rel_segs);
d2 = dists(2:end); r2 = a.outlet_radius(rel_segs);
clear d r
for s = 1 : length(d1)
    d(2*s-1) = d1(s); d(2*s) = d2(s); r(2*s-1) = r1(s); r(2*s) = r2(s);
end
% Plot actual
subplot(1,2,2)
plot(d,2*1000*r), hold on

%% Plot new actual (legs)

% Calculate actual
d1 = dists(1:end-1); r1 = new.inlet_radius(rel_segs);
d2 = dists(2:end); r2 = new.outlet_radius(rel_segs);
[~,order] = sort([d1;d2]);
for s = 1 : length(d1)
    d(2*s-1) = d1(s); d(2*s) = d2(s); r(2*s-1) = r1(s); r(2*s) = r2(s);
end
% Plot actual
subplot(1,2,2)
plot(d,2*1000*r)



%% Output new radii

for seg_no = [aorta_segs, leg_segs]
    fprintf(['\n new_inlet_radius(' num2str(seg_no) ') = ' num2str(new.inlet_radius(seg_no),4) ]);
end
for seg_no = [aorta_segs, leg_segs]
    fprintf(['\n new_outlet_radius(' num2str(seg_no) ') = ' num2str(new.outlet_radius(seg_no),4) ]);
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