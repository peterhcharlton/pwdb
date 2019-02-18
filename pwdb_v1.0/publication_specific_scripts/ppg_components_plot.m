function ppg_components_plot
% ppg_components_plot creates a plot of the PPG signal illustrating the
% origins of light attenuation (arterial blood, venous blood and other
% tissues).
%
%               ppg_components_plot
%
%	This file creates an image adapted from:
%           Peter H Charlton. (2018, August). Capitalising on Smart
%           Wearables to Improve Health Monitoring. Zenodo.
%           DOI: http://doi.org/10.5281/zenodo.1406011  
%   Please cite this publication when using this image.
%   
%   Output:
%       an EPS image in the same folder as this script
%     
%   Comments, Questions, Criticisms, Feedback, Contributions:
%       See: http://peterhcharlton.github.io/RRest/contributions.html
%
%   Licence:
%       please see the accompanying file named "LICENSE"
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function contains the path of the folder in which to store the data
% (which requires modification)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
up = setup_up;

create_folders(up);

download_data(up);

extract_data(up);

plot_data(up);

end

function up = setup_up

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         The path of the folder in which to store the data
%                    (which requires modification)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
up.paths.folders.root = ['/Users/petercharlton/Desktop/temp/ppg_components_signals', filesep];

% remaining paths
up.paths.folders.data = [up.paths.folders.root, 'raw_data', filesep];
up.paths.folders.plot = [up.paths.folders.root, 'plots', filesep];
up.paths.converted_data = [up.paths.folders.data, 'converted_data'];

% download paths
up.files.web = {'https://www.physionet.org/physiobank/database/bidmc/bidmc05'};
up.files.times.deb = [3];
up.files.times.fin = [8];
up.files.sigs.names = {{'ppg'}};
up.files.sigs.nos = {[2]};

close all

% Check that the WFDB Matlab Toolbox has been installed
if ~exist('getWfdbClass', 'file')
    error('Couldn''t find the WFDB Matlab Toolbox. Please install as described at the top of the file.')
end

% Check that the folder in which to store downloaded data exists
if ~exist(up.paths.folders.root, 'dir')
%     error('Couldn''t find the folder in which to store the data.')
end

end

function create_folders(up)

folders = fieldnames(up.paths.folders);
for folder_no = 1 : length(folders)
    eval(['curr_folder = up.paths.folders.' folders{folder_no} ';'])
    
    if ~exist(curr_folder)
        mkdir(curr_folder)
    end
   
end

end

function download_data(up)

exts = {'.hea', '.dat'};

for file_no = 1 : length(up.files.web)
    
    url = up.files.web{file_no};
    
    temp = strfind(url, '/');
    filename = url(temp(end)+1:end);
    filepath =  [up.paths.folders.data, filename];
    
    for ext_no = 1 : length(exts)
        curr_ext = exts{ext_no};
        expected_outfilename = [filepath,curr_ext];
        if ~exist(expected_outfilename, 'file')
            outfilename = websave([filepath,curr_ext],[url, curr_ext]);
        end
    end
    
end

end

function extract_data(up)

if exist([up.paths.converted_data, '.mat'])
    return
end

counter_no = 0;
for file_no = 1 : length(up.files.web)
    
    % Extract data from downloaded file
    url = up.files.web{file_no};
    temp = strfind(url, '/');
    filename = url(temp(end)+1:end);
    filepath =  [up.paths.folders.data, filename];
    cd(up.paths.folders.data)
    [signal,Fs,tm]=rdsamp(filename, [],[]);
    
    % Extract required data
    no_sigs = length(up.files.sigs.names{file_no});
    for sig_no = 1 : no_sigs
        curr_sig_no = up.files.sigs.nos{file_no}(sig_no);
        curr_times = [up.files.times.deb(file_no), up.files.times.fin(file_no)];
        rel_els = find(tm>= curr_times(1) & tm <= curr_times(2));
        curr_sig.v = signal(rel_els,curr_sig_no);
        curr_sig.t = [0:length(curr_sig.v)-1]./Fs;
        
        % store data
        counter_no = counter_no+1;
        data(counter_no).sig = curr_sig.v;
        data(counter_no).t = curr_sig.t;
        data(counter_no).db = strrep(url(temp(end-1)+1:temp(end)-1), 'db', '');
        data(counter_no).name = up.files.sigs.names{file_no}{sig_no};
    end
    
end

save(up.paths.converted_data, 'data')

end

function plot_data(up)


%% Make short PPG Plot

% setup
load(up.paths.converted_data)
ftsize = 24;
pos_long = [20,20,500,250];
pos_ppg = [20,20,600,600];
pos_short = [20,20,1000,470];

% Extract data
temp = extractfield(data, 'name');
sig_no = find(strcmp(temp, 'ppg')); clear temp
curr = data(sig_no);

curr.fix = 1*ones(size(curr.t));
curr.ven = curr.fix + 0.7 + 0.1*sin(curr.t*2*pi/4);
curr.sig = (curr.sig-min(curr.sig))/(max(curr.sig)-min(curr.sig));
curr.ppg = curr.ven + 0.5 + curr.sig(:)';

%plot(curr.t, curr.fix), hold on, plot(curr.t, curr.ven), plot(curr.t, curr.ppg)

% - Make figure (just PPG)
figure('Position', pos_short)
plot(curr.t, curr.ppg, 'k', 'LineWidth', 2)
ylim([0 3.31])
set(gca, 'FontSize', ftsize, 'YTick', [])
xlabel('Time (s)', 'FontSize', ftsize)
title('Raw Photoplethysmogram', 'FontSize', ftsize)
box off
savepath = [up.paths.folders.plot, 'short_raw_ppg_init'];
print(savepath, '-depsc')
print(savepath, '-dpng')
close all

% - Make figure (components)

% fill in
figure('Position', pos_short)
h = area(curr.t, curr.ppg, 'LineStyle', 'none'); hold on
h.FaceColor = [1,0,0];
plot(curr.t, curr.ppg, 'k', 'LineWidth', 2)
%plot(curr.t, curr.ven, '--k', 'LineWidth', 2)
h = area(curr.t, curr.ven, 'LineStyle', 'none'); hold on
h.FaceColor = [0,0,1];
plot(curr.t, curr.ven, '--k', 'LineWidth', 2)
h = area(curr.t, curr.fix, 'LineStyle', 'none'); hold on
h.FaceColor = [0.8,0.8,0.8];
plot(curr.t, curr.fix, '--k', 'LineWidth', 2)
ylim([0 3.31])
set(gca, 'FontSize', ftsize, 'YTick', [])
xlabel('Time (s)', 'FontSize', ftsize)
ylab = ylabel({'PPG','(au)'}, 'FontSize', ftsize, 'Rotation', 0);
set(ylab, 'position', get(ylab,'position')-[0.2,0.1,0]);
%title('Raw Photoplethysmogram', 'FontSize', ftsize)
box off

% annotate
dim = [.2 .19 .1 .1];
str = 'Other Tissues';
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'BackgroundColor', [1,1,1], 'FontSize', ftsize, 'HorizontalAlignment', 'center');
dim = [.2 .40 .1 .1];
str = 'Venous Blood';
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'BackgroundColor', [1,1,1], 'FontSize', ftsize, 'HorizontalAlignment', 'center');
dim = [.2 .58 .1 .1];
str = 'Arterial Blood';
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'BackgroundColor', [1,1,1], 'FontSize', ftsize, 'HorizontalAlignment', 'center');

savepath = [up.paths.folders.plot, 'short_raw_ppg'];
print(savepath, '-depsc')
print(savepath, '-dpng')
close all

end

function sig_info = get_sig_info(curr_name)

switch curr_name
    case 'ecg'
        sig_info.name = 'Electrocardiogram';
        sig_info.units = 'mV';
    case 'scg'
        sig_info.name = 'Seismocardiogram';
        sig_info.units = '';
    case 'ppg'
        sig_info.name = 'Photoplethysmogram';
        sig_info.units = '';
    case 'abp'
        sig_info.name = 'Arterial Blood Pressure';
        sig_info.units = '';
    case 'imp'
        sig_info.name = 'Impedance Pneumography';
        sig_info.units = '';
    case 'gyro'
        sig_info.name = 'Gyroscope';
        sig_info.units = '';
    case 'accel'
        sig_info.name = 'Acceleration';
        sig_info.units = '';
        
end

end
