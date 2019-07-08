function [cv_inds, fid_pts, pulses, sigs] = PulseAnalyse10(S, options)
% PULSEANALYSE  Extracts cardiovascular (CV) indices from pulse waves.
%   
%  Inputs:
%
%    S        -  a pulsatile signal, consisting of either a signal pulse, or a
%                 signal containing several pulses. S should be a structure, containing
%                 a vector of amplitudes, S.v, and the sampling frequency (in Hz), S.fs,
%                 and optionally the subject's height in metres, S.ht.
%    options  -  (optional) a structure of options, which may contain any of:
%                    options.exclude_low_quality_data   - a logical (true or false)
%                    options.do_plot                    - a logical (true or false)
%
%  Outputs:
%
%    cv_inds  -  a structure containing the calculate cardiovascular
%                   indices. For instance, cv_inds.AI.v is the value of the
%                   augmentation index (AI). If the input S contains
%                   several pulses, then  cv_inds.AI.v is the median AI
%                   value, and cv_inds.AI.raw is a vector of raw values for
%                   each pulse.
%    fid_pts  -  a structure containing the indices of the fiducial points.
%                   For instance, fid_pts.dic provides the indices of the
%                   dicrotic notch(es).
%    pulses   -  a structure containing information on the pulses:
%                   pulses.peaks    - the indices of systolic peaks
%                   pulses.onsets   - the indices of pulse onsets (i.e. beginning of systole)
%                   pulses.quality  - a logical indicating the signal quality of each pulse (true indicates high quality)
%    S_filt   -  the filtered pulsatile signal to which the "fid_pts" and
%                   "pulses" indices correspond.
%
%  Exemplary usage:
%
%    cv_inds = PulseAnalyse(S)                                  extracts CV indices from the pulsatile signal, S.
%    cv_inds = PulseAnalyse(S, options)                         uses options to specify the analysis.
%    [cv_inds, fid_pts, pulses, S_filt] = PulseAnalyse(___)     also outputs fiducial points, pulse locations, and the filtered signal.
%
%  For further information please see the accompanying manual.
%
%  This script contains items either copied or modified from the RRest
%  toolbox which is covered by the GNU public licence (<a href="http://github.com/peterhcharlton/RRest/">link</a>).
%
% Peter H. Charlton, King's College London, August 2017

%% Setup

% Setup options
if nargin < 2
    options = struct;
end
options = setup_options(options);

% Setup universal parameters
up = setup_up(options);

% Determine whether this is a single pulse, or several pulses
no_of_pulses = determine_no_of_pulses(S, up);

% Make amplitudes into a column vector
S.v = S.v(:);

if sum(strcmp(fieldnames(S), 'ht'))
    ht = S.ht;
end

%% Pre-processing

% Pre-processing is peformed according to whether this is a single pulse or multiple pulses.
switch no_of_pulses
    
    case 'multiple'
        
        % Eliminate very low frequency content
        if options.do_filter
            try
                S_evlf = elim_vlfs(S, up);
            catch
                S_evlf = S;
                fprintf('Signal too short to eliminate very low frequencies\n')
            end
            
            % Eliminate very high frequency content
            S_filt = elim_vhfs(S_evlf, up);
        else
            S_filt = S;
        end
        
        % Identify individual pulse waves
        [pulses.peaks, pulses.onsets, ~] = adaptPulseSegment(S_filt.v,S.fs);
        
        % Assess signal quality
        pulses.quality = assess_signal_quality(S_filt, pulses.peaks);
        
    case 'single'
        
        if options.do_filter
            % Eliminate high frequency content
            S_evhf = elim_vhfs2(S, up);
            
            % Eliminate low frequency content
            S_elim_lf = eliminate_low_freq_from_single_beat(S_evhf, up);
        else
            S_elim_lf = S;
        end
        
        % Ensure that signal commences at start of systole
        S_filt = align_pulse(S_elim_lf, up);
        S = align_pulse(S, up);
        
        % Generate "pulses" information
        pulses = generate_pulses_for_single_beat(S_filt);
        
end

%% Calculate derivatives

sigs = calc_derivs2(S, S_filt, no_of_pulses, up);
% derivs = calc_derivs(S_filt, no_of_pulses, up);

%% Identify fiducial points

if exist('ht', 'var')
    sigs.ht = ht;
end

if ~options.annotate
    try
        % Find fiducial points
        fid_pts = identify_fiducial_pts(sigs, pulses.onsets, up, no_of_pulses, options);
    catch
        fprintf('\n Couldn''t identify fiducial points')
        
        [fid_pts.a, fid_pts.b, fid_pts.c, fid_pts.d, fid_pts.e, fid_pts.f, fid_pts.p1pk, fid_pts.p2pk, fid_pts.p1in, fid_pts.t2t, fid_pts.s, fid_pts.dic, fid_pts.dia, fid_pts.ms] = deal(nan(length(pulses.onsets)-1,1));
        fid_pts.f1 = pulses.onsets(1:end-1);
        fid_pts.f2 = pulses.onsets(2:end);
    end
else
    % Annotate fiducial points
    fid_pts = annotate_fiducial_pts(S_filt, pulses.onsets, derivs, up, no_of_pulses);
end

%% Calculate CV indices

cv_inds = calc_stiffness_inds(sigs, pulses, fid_pts, options);

%% Plot the fiducial points
if options.do_plot
    plot_fiducial_points(sigs, pulses, fid_pts, no_of_pulses, options)
    plot_cv_inds(sigs, pulses, fid_pts, cv_inds, no_of_pulses, options)
end

%% Make a demo plot
if options.do_demo_plot
    make_demo_plot(S_filt, derivs, pulses, fid_pts, cv_inds, no_of_pulses, options)
end

end

function options = setup_options(options)

if isempty(fieldnames(options)) | ~strcmp(fieldnames(options), 'do_plot')
    options.do_plot = 1;
end

if ~strcmp(fieldnames(options), 'exclude_low_quality_data')
    options.exclude_low_quality_data = 1;
end
if ~strcmp(fieldnames(options), 'do_plot')
    options.do_plot = 1;
end
if ~strcmp(fieldnames(options), 'plot_third_deriv')
    options.plot_third_deriv = 1;
end
if ~strcmp(fieldnames(options), 'annotate')
    options.annotate = 0;
end
if ~strcmp(fieldnames(options), 'manual_adjustment')
    options.manual_adjustment = 0;
end
if ~strcmp(fieldnames(options), 'close_figures')
    options.close_figures = 1;
end
if ~strcmp(fieldnames(options), 'do_filter')
    options.do_filter = 1;
end
if ~strcmp(fieldnames(options), 'save_folder')
    options.save_folder = '';
end
if ~strcmp(fieldnames(options), 'save_file')
    options.save_file = '';
end
if ~strcmp(fieldnames(options), 'do_demo_plot')
    options.do_demo_plot = 0;
end
if ~strcmp(fieldnames(options), 'demo_plot_wave')
    options.demo_plot_wave = 'orig';
end

if options.do_demo_plot
    options.do_plot = 0;
end

end

function up = setup_up(options)

if options.close_figures
    close all
end

%% Analysis settings

% Threshold signal duration to distinguish between a single pulse or multiple pulses: 
up.analysis.max_duration_of_single_pulse = 2.5;   % in secs

% Filter characteristics: Eliminate VLFs (below resp freqs): For 4bpm cutoff
up.paramSet.elim_vlf.Fpass = 0.157;  % in Hz
up.paramSet.elim_vlf.Fstop = 0.02;   % in Hz     (0.157 and 0.02 provide a - 3dB cutoff of 0.0665 Hz)
up.paramSet.elim_vlf.Dpass = 0.05;
up.paramSet.elim_vlf.Dstop = 0.01;

% Filter characteristics: Eliminate VHFs (above frequency content of signals)
% up.paramSet.elim_vhf.Fpass = 38.5;  % in HZ
% up.paramSet.elim_vhf.Fstop = 33.12;  % in HZ   (33.12 and 38.5 provide a -3 dB cutoff of 35 Hz)
up.paramSet.elim_vhf.Fpass = 20;  % in HZ
up.paramSet.elim_vhf.Fstop = 15;  % in HZ
up.paramSet.elim_vhf.Dpass = 0.05;
up.paramSet.elim_vhf.Dstop = 0.01;

% No of times to repeat a single pulse to perform VHF filtering
up.paramSet.no_pulse_repeats = 5;

end

function no_of_pulses = determine_no_of_pulses(S, up)
% DETERMINE_NO_OF_PULSES  Determines whether this is a single pulse wave,
% of a pulsatile signal containing multiple pulse waves.

signal_duration = (length(S.v)-1)/S.fs;

% If duration of signal is greater than a threshold then assume this is
% multiple pulse waves:
if signal_duration > up.analysis.max_duration_of_single_pulse
    no_of_pulses = 'multiple';
else
    no_of_pulses = 'single';
end

end

function s_filt = elim_vlfs(s, up)
%% Filter pre-processed signal to remove frequencies below resp
% Adapted from RRest

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Make filter
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord([up.paramSet.elim_vlf.Fstop up.paramSet.elim_vlf.Fpass]/(s.fs/2), [1 0], [up.paramSet.elim_vlf.Dstop up.paramSet.elim_vlf.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% % Gives a -3 dB cutoff at ? Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.0266;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(fs/2);

s_filt.v = filtfilt(AMfilter.numerator, 1, s.v);
s_filt.v = s.v-s_filt.v;
s_filt.fs = s.fs;
end

function s_filt = elim_vhfs(s, up)
%% Filter signal to remove VHFs
% Adapted from RRest

s_filt.fs = s.fs;

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Check to see if sampling freq is at least twice the freq of interest
if (up.paramSet.elim_vhf.Fpass/(s.fs/2)) >= 1
    % then the fs is too low to perform this filtering
    s_filt.v = s.v;
    return
end

%% Create filter
% parameters for the low-pass filter to be used
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord([up.paramSet.elim_vhf.Fstop up.paramSet.elim_vhf.Fpass]/(s.fs/2), [1 0], [up.paramSet.elim_vhf.Dstop up.paramSet.elim_vhf.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% Gives a -3 dB cutoff at cutoff_freq Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.139;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(fs/2);

%% Remove VHFs
s_dt=detrend(s.v);
s_filt.v = filtfilt(AMfilter.numerator, 1, s_dt);
end

function s_filt = elim_vhfs2(s, up)
%% Filter signal to remove VHFs
% Adapted from RRest
% Adapted for single pulses

s_filt.fs = s.fs;

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Repeat pulse
s.v = repmat(s.v(:), [up.paramSet.no_pulse_repeats,1]);

%% Check to see if sampling freq is at least twice the freq of interest
if (up.paramSet.elim_vhf.Fpass/(s.fs/2)) >= 1
    % then the fs is too low to perform this filtering
    s_filt.v = s.v;
    return
end

%% Create filter
% parameters for the low-pass filter to be used
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord([up.paramSet.elim_vhf.Fstop up.paramSet.elim_vhf.Fpass]/(s.fs/2), [1 0], [up.paramSet.elim_vhf.Dstop up.paramSet.elim_vhf.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% Gives a -3 dB cutoff at cutoff_freq Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.067;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(s.fs/2);

%% Remove VHFs
s_dt=detrend(s.v);
s_filt.v = filtfilt(AMfilter.numerator, 1, s_dt);

%% Extract original pulse (from repeated pulses)
len_of_each_pulse = length(s.v)/up.paramSet.no_pulse_repeats;
start_pulse_no = floor(up.paramSet.no_pulse_repeats/2);
end_pulse_no = ceil(up.paramSet.no_pulse_repeats/2);
rel_els = (start_pulse_no*len_of_each_pulse)+1 : end_pulse_no*len_of_each_pulse;
s_filt.v = s_filt.v(rel_els);

end

function [peaks,onsets,clipp] = adaptPulseSegment(y,Fs,annot)
%ADAPTPULSESEGMENT perform adaptive pulse segmentation and artifact detection 
%in ppg signals
%   [peaks,onsets,artif] = adaptPulseSegment(y,annot)
%
% Inputs:
%       y      vector, ppg signal [Lx1] or [1xL], in which L = length(signal)
%       Fs      scalar, sampling frequency in Hz
%       annot   vector, timestamps (in samples) with location of the peaks
%
% Outputs:
%       peaks   vector, locations of peaks (in samples)
%       onsets  vector, locations of onsets of the beats (in samples)
%       artif   vector, locations of peaks classified as artefacts
% 
% References:
%       Karlen et al., Adaptive Pulse Segmentation and Artifact Detection in 
%       Photoplethysmography for Mobile Applications, 34th Internat. Conf. IEEE-EMBS 2012
%       
% Written by Marco A. Pimentel

doOptimise = 1;
doPlot = 0;
if nargin < 3
    % no annotations are provided, therefore, no optimisation will take
    % place
    doOptimise = 0; 
end


% The algorithm in the paper is applied to signals sampled at 125 Hz...
% We do not resample our signal
%ys = resample(y,125,Fs);
%Fs = 125;
% if Fs ~= 300
%     ys = resample(y,300,Fs);
%     Fs = 300;
% else
    ys = y;
% end

% The paper is not clear about the selection of the range of search for m
% We define the range of m to be between [0.005 - 0.100] secs (5ms to 100ms)
% We define "m" in terms of samples
opt.bounds = 0.005:0.005:0.100;
opt.m = unique(ceil(opt.bounds*Fs));

opt.perf = zeros(length(opt.m),4); % store results of performance for each m

if doOptimise
    % Perform optimisation
    for i = 1 : length(opt.m)
        % Determine peaks and beat onsets
        [linez,linezSig] = pulseSegment(ys,Fs,opt.m(i));
        % Calculate performance of the peak detection
        opt.perf(i,:) = evalPerf(annot,linez(:,2));
    end
    
else
    % Do not perform optimization; fix m
    opt.m = 10;
    [peaks,artifs,clipp,linez,linezSig] = pulseSegment(ys,Fs,opt.m);
end

if doPlot
    colData = {'g','y','r'};
    figure; 
    h(1) = subplot(211);
    plot(ys); hold on;
    for i = 1 : size(linez,1)
        %if linezSig(i) > -1
        plot(linez(i,:),ys(linez(i,:)),'-x','Color',colData{linezSig(i)+2});
        %end
    end
    
    h(2) = subplot(212);
    plot(ys,'g'); hold on;
    for i = 1 : size(peaks,1)
        plot(peaks(i,:),ys(peaks(i,:)),'-xr');
    end
    if ~isempty(artifs)
    for i = 1 : size(artifs,1)
        plot(artifs(i,:),ys(artifs(i,:)),'--^b');
    end
    end
    if ~isempty(clipp)
    for i = 1 : size(clipp,1)
        plot(clipp(i,:),ys(clipp(i,:)),'-om');
    end
    end
    linkaxes(h,'x');
    
end

% Correct for the downsmapling performed during the peak detection
onsets = peaks(:,1);
peaks  = peaks(:,2);
for i = 1 : size(peaks,1)
    [~,ind]  = min(ys(max([1 onsets(i)-opt.m]):min([length(ys) onsets(i)+opt.m])));
    onsets(i) = max([1 onsets(i)-opt.m]) + ind(1) - 1;
    [~,ind]  = max(ys(max([1 peaks(i)-opt.m]):min([length(ys) peaks(i)+opt.m])));
    peaks(i) = max([1 peaks(i)-opt.m]) + median(ind) - 1;
end

% Correct minimum value of onset of the beat
for i = 2 : length(onsets)
    [~,ind]   = min(ys(peaks(i-1):peaks(i)));
    onsets(i) = peaks(i-1) + ind - 1;
end

end

function [peaks,artifs,clipp,linez,linezSig] = pulseSegment(ys,Fs,m)
% Perform pulse segmentation in ys given m
% Inputs:
%       ys      vector, with ppg signal
%       m       scalar, with the length of each line (read paper for details)
% 
% Outputs:
%       linez      2-column vector, with location of beat onsets and peaks
%       linezSig   vector, with label of the slope of each line segment
%                  1 - positive slope; -1 - negative slope; 0 - constant
% 

% split signal in different segments
nseg = floor(length(ys)/m);    % number of segments
% nseg = round(length(ys)/m);    % number of segments
% intialize loop variables
seg = 1;    % segment counter
z = 1;      % line counter 
segInLine = 1;  % line controler
linez = zeros(nseg,2); linez(1,:) = [1,m];
% slope of segment/line
a = zeros(nseg,1); a(1) = slope(ys,linez(1,:));
% "classify" line segment according to the slope
linezSig = zeros(nseg,1); linezSig(1) = sign(a(1));
% Start loop over segments
z = z + 1; seg = seg + 1;
for i = 2 : nseg    % loop over segments
    linez(z,:) = [(seg-1)*m+1 seg*m];
    try
        a(z) = slope(ys,linez(z,:));
    catch
        a = 1;
    end
    linezSig(z) = sign(a(z));
    if sign(a(z)) == sign(a(z-1))
        linez(z-1,:) = [(seg-1-segInLine)*m+1 seg*m];
        seg = seg + 1;
        segInLine = segInLine + 1;
    else
        z = z + 1;
        seg = seg + 1;
        segInLine = 1;
    end
end

% remove extra spaces created in output variables
linezSig(sum(linez,2)==0,:) = [];
linez(sum(linez,2)==0,:) = [];

% Apply adaptive threshold algorithm
% For this algorithm to work, we need to first find a valide line segment 
% in order to intialize the thresholds! In order to this, we define a flag
% to control the intialization in the main loop
FOUND_L1 = 0;

% The algorithm includes the definition of 4 adaptation parameters
% We define the following adaptation parameters
% a =     | a_fast_low    a_fast_high |
%         | a_slow_low    a_slow_high |
% 
a = [0.5 1.6; ...
     0.6 2.0];
 
% Define fixed thresholds described in the paper
ThT  = 0.03 * Fs;    % Duration of the line
ThIB = 0.24 * Fs;    % Interbeat invertal (240 ms) 

% Define parameters used in the main loop
alpha = zeros(size(linez,1),1);
for i = 1 : size(linez,1)
    alpha(i) = slope(ys,linez(i,:));   % slopes of line segments
end
theta = diff(ys(linez),[],2);
durat = diff(linez,[],2);       % duration of line segments (in samples)

% remove lines that do not have the necessary duration
linez(durat<ThT,:) = [];
theta(durat<ThT,:) = [];
alpha(durat<ThT,:) = [];
horiz = horizontalLine(ys,linez,Fs);

FLAG = 0;
artifs = []; clipp = [];
% Select window for detect firs peaks!
wind = theta(theta>0);
try 
    wind = wind(1:10);
catch
    wind = wind;
end
ThAlow  = prctile(wind,95)*0.6;
ThAhigh = prctile(wind,95)*1.8;
peaks = [];
for z = 1 : size(linez,1)-1   % loop over line segments
    if FOUND_L1
        if alpha(z) > 0 && ... % slope must be positive
                alpha(z-1) ~= 0 && ...  % peaks before or after clipping are artefactual
                alpha(z+1) ~= 0
            if theta(z) >= ThAlow && theta(z) <= ThAhigh && ...
                    linez(z,2) >= peaks(end,2) + ThIB
                ThAlow  = (ThAlow + theta(z)*a(2,1))/2;
                ThAhigh = theta(z) * a(2,2);
                FLAG = 0;
                currTheta = [currTheta; theta(z)];
                peaks = [peaks; linez(z,:)];
            else
                if FLAG > 0
                    ThAlow  = (ThAlow + min(currTheta(max([1 end-4]):end))*a(1,1))/2;
                    ThAhigh = max(currTheta(max([1 end-4]):end)) * a(1,2);
                    %ThAlow  = (ThAlow + theta(z)*a(1,1))/2;
                    %ThAhigh = theta(z) * a(1,2);
                end
                FLAG = FLAG + 1;
                artifs = [artifs; linez(z,:)];
            end
        elseif theta(z) > 0 && ... 
                ((theta(z-1) ~= 0 || horiz(z-1) ~= 0) && ...
                (theta(z+1) ~= 0 || horiz(z+1) ~= 0))
            if theta(z) >= ThAlow && theta(z) <= ThAhigh && ...
                    linez(z,2) >= peaks(end,2) + ThIB
                ThAlow  = (ThAlow + theta(z)*a(2,1))/2;
                ThAhigh = theta(z) * a(2,2);
                FLAG = 0;
                currTheta = [currTheta; theta(z)];
                peaks = [peaks; linez(z,:)];
            else
                if FLAG > 0
                    %ThAlow  = (ThAlow + currTheta*a(1,1))/2;
                    %ThAhigh = currTheta * a(1,2);
                    ThAlow  = (ThAlow + min(currTheta(max([1 end-4]):end))*a(1,1))/2;
                    ThAhigh = max(currTheta(max([1 end-4]):end)) * a(1,2);
                    %ThAlow  = (ThAlow + theta(z)*a(1,1))/2;
                    %ThAhigh = theta(z) * a(1,2);
                end
                FLAG = FLAG + 1;
                artifs = [artifs; linez(z,:)];
            end
        elseif theta(z) == 0 && horiz(z) == 0
            artifs  = [artifs; linez(z,:)];
            clipp   = [clipp; linez(z,:)];
        end 
    else
        if alpha(z) > 0 && durat(z) >= ThT && ...
                theta(z) >= ThAlow && theta(z) <= ThAhigh 
            FOUND_L1 = 1;
            ThAlow  = theta(z)*0.5;
            ThAhigh = theta(z)*2.0;
            peaks = linez(z,:);    % loaction of onsets and peaks
            currTheta = theta(z);
        end
    end
end

end

function out = horizontalLine(ys,linez,Fs)
% Get horizontal lines from signal given linez
out = zeros(size(linez,1),1);
for i = 1 : size(linez,1)
    out(i) = median(abs(diff(ys(linez(i,1):linez(i,2)))));
    % check duration of the peaks
    if out(i) == 0 && diff(linez(i,:)) <= 0.200*Fs
        out(i) = 0.1;
    end
end

end

function out = slope(ys,interv)
start = interv(1); stop = interv(2);
out = sum(diff(ys([start:stop])))/(stop-start);
%out = median(gradient(ys(start:stop)));
end

function quality = assess_signal_quality(s, pulse_inds)
% ASSESS_SIGNAL_QUALITY  Assesses the signal quality of each beat of the
% pulsatile signal.
% Inputs:
%       s           -  pulsatile signal, a structure containing s.v (a
%                       vector of values), and s.fs (sampling frequency in Hz).
%       pulse_inds  -  indices of the pulse peaks
%
% Outputs:
%       quality     -  the signal quality of each beat (1 indicates high
%                       quality, 0 low quality).
%
% Adapted from RRest
%
% Reference: This function uses an adaptation of the signal quality index
% for the photoplethysmogram described in:
%     Orphanidou, C. et al., 2015. Signal-quality indices for the electrocardiogram and photoplethysmogram: derivation and applications to wireless monitoring. IEEE Journal of Biomedical and Health Informatics, 19(3), pp.832–8. Available at: http://www.ncbi.nlm.nih.gov/pubmed/25069129.

%% Setup
s.t = [0:length(s.v)-1]/s.fs;

%% Segment into windows of 10s duration
win_durn = 10;   % in secs
win.deb = s.t(1):(win_durn-2):s.t(end);
win.fin = win.deb + win_durn;

high_quality_pulse_inds = [];
for win_no = 1 : length(win.deb)
    
    % identify data for this window
    
    rel_els = s.t >= win.deb(win_no) & s.t <= win.fin(win_no);
    first_el = find(rel_els,1);
    curr_sig.v = s.v(rel_els); 
    curr_sig.t = s.t(rel_els); clear rel_els
    curr_sig.t = curr_sig.t - curr_sig.t(1);
    curr_sig.pulse_ind_inds = find(s.t(pulse_inds) >= win.deb(win_no) & s.t(pulse_inds) <= win.fin(win_no));
    curr_pulse_inds = pulse_inds(curr_sig.pulse_ind_inds) - first_el + 1;
    
    % find beat-to-beat interval
    
    beat_to_beat_interval = median(diff(curr_sig.t(curr_pulse_inds)));
    beat_to_beat_samples = round(beat_to_beat_interval*s.fs); clear beat_to_beat_interval
    
    % find a template beat
    ts = [];
    rel_els = curr_pulse_inds>beat_to_beat_samples/2 & ...
        curr_pulse_inds+floor(beat_to_beat_samples/2)<length(curr_sig.v);
    rel_pulse_inds = curr_pulse_inds(rel_els);
    curr_sig.used_pulse_ind_inds = curr_sig.pulse_ind_inds(rel_els);
    % find beat morphologies
    for rel_pulse_no = 1 : length(rel_pulse_inds)
        this_pulse = rel_pulse_inds(rel_pulse_no);
        t = curr_sig.v(this_pulse-floor(beat_to_beat_samples/2):this_pulse+floor(beat_to_beat_samples/2));
        tt = t/norm(t); tt = tt(:)';
        ts = [ts; tt]; clear tt t
    end
    clear k l j
    
    % find all templates in current window
    avtempl = mean(ts,1);
    
    % now calculate correlation for every beat in this window
    r2 = nan(size(ts,1),1);
    for k = 1:size(ts,1)
        r2(k) = corr2(avtempl,ts(k,:));
    end
    clear k
    
    high_quality_beats = r2 > 0.86;
    
    high_quality_pulse_inds = [high_quality_pulse_inds, curr_sig.used_pulse_ind_inds(high_quality_beats)];
    
    % % calculate template using only high-quality beats
    % templ = mean(ts(r2>0.86,:),1);  %threshold = 0.66 for ECG, 0.86 for PPG;    % ecg cross-correlation threshold value (for sqi)
    
end

high_quality_pulse_inds = unique(high_quality_pulse_inds);
quality = false(length(pulse_inds),1);
for ind = 1 : length(quality)
    if intersect(high_quality_pulse_inds, ind)
        quality(ind) = true;
    end
end
end

function sigs = calc_derivs2(sig, s_filt, no_of_pulses, ~)

switch no_of_pulses
    
    case 'multiple'
        
        % calculate derivatives
        derivs.first = fsg521(sig.v)*sig.fs;
        derivs.second = fsg521(derivs.first)*sig.fs;
        derivs.third = fsg521(derivs.second)*sig.fs;
        
%         % calculate derivatives of filtered signal
%         derivs.f_first = savitzky_golay(s_filt.v, 1, 9);
%         derivs.f_second = savitzky_golay(derivs.f_first, 1, 9);
%         derivs.f_third = savitzky_golay(derivs.f_second, 1, 9);

    case 'single'
        
        % repeat this single pulse several times
        no_pulses = 3; % number of times to repeat pulse
        sig.v = repmat(sig.v, [no_pulses,1]);
        s_filt.v = repmat(s_filt.v, [no_pulses,1]);
        
        % calculate derivatives of orig signal
        derivs.first = savitzky_golay(sig.v, 1, 5)*sig.fs;
        derivs.second = savitzky_golay(derivs.first, 1, 5)*sig.fs;
        derivs.third = savitzky_golay(derivs.second, 1, 5)*sig.fs;
        
%         % calculate derivatives of filtered signal
%         derivs.f_first = savitzky_golay(s_filt.v, 1, 9);
%         derivs.f_second = savitzky_golay(derivs.f_first, 1, 9);
%         derivs.f_third = savitzky_golay(derivs.f_second, 1, 9);
        
        % select derivative values corresponding to the original pulse
        orig_len = length(sig.v)/no_pulses;
        orig_els = (floor(no_pulses/2)*orig_len)+1 : ((floor(no_pulses/2)+1)*orig_len);
        derivs.first = derivs.first(orig_els);
        derivs.second = derivs.second(orig_els);
        derivs.third = derivs.third(orig_els);
%         derivs.f_first = derivs.f_first(orig_els);
%         derivs.f_second = derivs.f_second(orig_els);
%         derivs.f_third = derivs.f_third(orig_els);
end

sigs.fs = sig.fs;
sigs.s = sig.v;
sigs.s_filt = s_filt.v;
sigs.derivs = derivs;

end

function deriv = savitzky_golay(sig, deriv_no, win_size)

%% assign coefficients
% From: https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter#Tables_of_selected_convolution_coefficients
% which are calculated from: A., Gorry (1990). "General least-squares smoothing and differentiation by the convolution (Savitzky?Golay) method". Analytical Chemistry. 62 (6): 570?3. doi:10.1021/ac00205a007.

switch deriv_no
    case 0
        % - smoothing
        switch win_size
            case 5
                coeffs = [-3, 12, 17, 12, -3];
                norm_factor = 35;
            case 7
                coeffs = [-2, 3, 6, 7, 6, 3, -2];
                norm_factor = 21;
            case 9
                coeffs = [-21, 14, 39, 54, 59, 54, 39, 14, -21];
                norm_factor = 231;
            otherwise
                error('Can''t do this window size')
        end
    case 1
        % - first derivative
        switch win_size
            case 5
                coeffs = -2:2;
                norm_factor = 10;
            case 7
                coeffs = -3:3;
                norm_factor = 28;
            case 9
                coeffs = -4:4;
                norm_factor = 60;
            otherwise
                error('Can''t do this window size')
        end
        
    case 2
        % - second derivative
        switch win_size
            case 5
                coeffs = [2,-1,-2,-1,2];
                norm_factor = 7;
            case 7
                coeffs = [5,0,-3,-4,-3,0,5];
                norm_factor = 42;
            case 9
                coeffs = [28,7,-8,-17,-20,-17,-8,7,28];
                norm_factor = 462;
            otherwise
                error('Can''t do this window size')
        end
        
    case 3
        % - third derivative
        switch win_size
            case 5
                coeffs = [-1,2,0,-2,1];
                norm_factor = 2;
            case 7
                coeffs = [-1,1,1,0,-1,-1,1];
                norm_factor = 6;
            case 9
                coeffs = [-14,7,13,9,0,-9,-13,-7,14];
                norm_factor = 198;
            otherwise
                error('Can''t do this window size')
        end
        
    case 4
        % - fourth derivative
        switch win_size
            case 7
                coeffs = [3,-7,1,6,1,-7,3];
                norm_factor = 11;
            case 9 
                coeffs = [14,-21,-11,9,18,9,-11,-21,14];
                norm_factor = 143;
            otherwise
                error('Can''t do this window size')
        end
        
    otherwise
        error('Can''t do this order of derivative')        
end

if rem(deriv_no, 2) == 1
    coeffs = -1*coeffs;
end

A = [1,0];
filtered_sig = filter(coeffs, A, sig);
s=length(sig);
half_win_size = floor(win_size*0.5);
deriv=[filtered_sig(win_size)*ones(half_win_size,1);filtered_sig(win_size:s);filtered_sig(s)*ones(half_win_size,1)];
deriv = deriv/norm_factor;

end

function derivs = calc_derivs(sig, no_of_pulses, up)

switch no_of_pulses
    
    case 'multiple'
        
        % calculate derivatives
        derivs.first = fsg521(sig.v);
        derivs.second = fsg521(derivs.first);
        derivs.third = fsg521(derivs.second);

    case 'single'
        
        % repeat this single pulse several times
        no_pulses = 5; % number of times to repeat pulse
        sig.v = repmat(sig.v, [no_pulses,1]);
        
        % calculate derivatives
        derivs.first = fsg521(sig.v);
        derivs.second = fsg521(derivs.first);
        derivs.third = fsg521(derivs.second);
        
        % select derivative values corresponding to the original pulse
        orig_len = length(sig.v)/no_pulses;
        orig_els = (floor(no_pulses/2)*orig_len)+1 : ((floor(no_pulses/2)+1)*orig_len);
        derivs.first = derivs.first(orig_els);
        derivs.second = derivs.second(orig_els);
        derivs.third = derivs.third(orig_els);
        
end

derivs.fs = sig.fs;

end

function dx=fsg521(x)

% Savitsky-Golay filter function
% dx=fsg521(x)
% 	5 point SavGol filter, 2nd order polynomial, 1st derivative
%	input 	x
%	output	dx
%	corrected for time shift
%
% Adapted from Marie Willemet's code

C=[0.2,0.1];
	
for i=1:2; 
	B(i)=C(i); 
end	
B(3)=0.0;
for i=4:5;
	B(i)=-C(6-i);
end
A=[1,0];
	
s=length(x);
dx=filter(B,A,x);
dx=[dx(5)*ones(2,1);dx(5:s);dx(s)*ones(2,1)];

end

function pts = identify_fiducial_pts(sigs, onsets, up, no_of_pulses, options)
% IDENTIFY_FIDUCIAL_PTS  Identifies fiducial points on each wave of a
% pulsatile signal.
%
% Inputs:
%       sig         -  pulsatile signal, a structure containing .v (a
%                       vector of values), and .fs (sampling frequency in Hz).
%       onsets      -  indices of the pulse onsets
%       derivs      -  a structure containing .first, .second and .third
%                       derivatives of sig
%
% Outputs:
%       pts         -  a structure containing the indices of the fiducial
%                       points for each beat (a nan indicates that a fiducial
%                       point could not be identified).
%


%% find fiducial points

% setup variables
fid_pt_names = {'a', 'b', 'c', 'd', 'e', 'f', 's', 'dia', 'dic', 'p1pk', 'p2pk', 'p1in', 'p2in', 'ms', 'f1', 'f2'};
for fid_pt_no = 1 : length(fid_pt_names)
    eval(['pts.' fid_pt_names{fid_pt_no} ' = nan(length(onsets)-1,1);'])
end

% setup buffer zones
buffer_tol.deb = 0.005; % in secs
buffer_tol.fin = 0.2; % in proportion of beat
buffer_p1 = [0.1, 0.18]; % in secs

for pulse_no = 1 : length(onsets)-1
    
    % extract data for this pulse
    curr_els = onsets(pulse_no):onsets(pulse_no+1);
    curr.f_sig = sigs.s_filt(curr_els);
    curr.sig = sigs.s(curr_els);
    
    if strcmp(no_of_pulses, 'multiple')
        
        old.v = curr.sig; old.fs = sigs.fs;
        % Eliminate low frequency content
        temp = eliminate_low_freq_from_single_beat(old, up);
        
        % Ensure that signal commences at start of systole
        temp2 = align_pulse(temp, up);
        curr.sig = temp2.v;
    end
    
    % calculate derivatives
    curr.derivs.first = sigs.derivs.first(curr_els);
    curr.derivs.second = sigs.derivs.second(curr_els);
    curr.derivs.third = sigs.derivs.third(curr_els);
%     curr.derivs.f_first = sigs.derivs.f_first(curr_els);
%     curr.derivs.f_second = sigs.derivs.f_second(curr_els);
%     curr.derivs.f_third = sigs.derivs.f_third(curr_els);
    
    % find buffers
    initial_buffer = floor(buffer_tol.deb * sigs.fs);  % in samples
    end_buffer = length(curr.sig) - ceil(buffer_tol.fin * length(curr.sig));  % in samples
    
    % find f1 and f2
    temp_f1 = 1;
    temp_f2 = length(curr.sig);
    
    % find s
    temp_s = identify_s(curr);
    
    % find ms
    temp_ms = identify_ms(curr);
    
    % find a
    temp_a = identify_a(curr, initial_buffer);
    
    if isempty(temp_a)
        continue
    end
    
    % find b
    temp_b = identify_b(curr, temp_a);
    
    if isempty(temp_b)
        continue
    end
    
    % find p1
    p1_buffer = floor(buffer_p1 .* sigs.fs);  % in samples
    temp_p1 = identify_p1(curr, temp_b, temp_ms, p1_buffer, no_of_pulses, sigs.fs, up);
    
    % find e
    redo_log = 0;
    temp_e = identify_e(curr, temp_s, temp_ms, temp_b, end_buffer, redo_log);
    
    if options.manual_adjustment
        if sum(strcmp(options, 'e'))
            temp_e = options.man_e;
        else
            temp_e = 270;
        end
    end
    
    if isempty(temp_e)
        [temp_f, temp_c, temp_d, temp_dic, temp_dia, temp_p2] = deal(nan);
    else
        
        % find c
        temp_c = identify_c(curr, temp_b, temp_e);
        
        if isempty(temp_c)
            redo_log = 1;
            old_temp_e = temp_e;
            temp_e = identify_e(curr, temp_s, temp_ms, end_buffer, redo_log);
        end
        
        % find c
        temp_c = identify_c(curr, temp_b, temp_e);
        
        % find f
        temp_f = identify_f(curr, temp_e, end_buffer);
        
        % find dic
        temp_dic = identify_dic(curr,temp_e);
        
        % find dia
        temp_dia = identify_dia(curr, temp_dic, temp_e, end_buffer);
        
        if isempty(temp_c)
            [temp_d, temp_p2] = deal(nan);
        else
            % find d
            temp_d = identify_d(curr, temp_c, temp_e);
            
            % find p2
            temp_p2 = identify_p2(curr, temp_d, temp_p1, temp_dic);
            
        end
        
    end
    
    % retain timings of original p1 and p2 estimates
    temp_p1in = temp_p1; temp_p2in = temp_p2;
    
    % make p1 or p2 coincident with the systolic peak
    [~, rel_el] = min(abs(temp_s-[temp_p1,temp_p2]));
    if rel_el == 1
        temp_p1 = temp_s;
    else
        temp_p2 = temp_s;
    end
    
    if ~isnan(temp_p2) & ~isnan(temp_p1)
        % make sure p2 is at a peak if necessary
        pks = find_pks_trs(curr.sig, 'pk');
        cutoff = mean([temp_p1, temp_p2]);
        possible_pks = find(pks > cutoff & pks < temp_e & curr.sig(pks) > curr.sig(temp_p2));
        if ~isempty(possible_pks)
            [~, temp_el] = max(curr.sig(pks(possible_pks)));
            temp_p2 = pks(possible_pks(temp_el));
        end
        % make sure p1 is at a peak if necessary
        pks = find_pks_trs(curr.sig, 'pk');
        cutoff = mean([temp_p1, temp_p2]);
        possible_pks = find(pks < cutoff & pks > temp_ms & curr.sig(pks) > curr.sig(temp_p1));
        if ~isempty(possible_pks)
            [~, temp_el] = max(curr.sig(pks(possible_pks)));
            temp_p1 = pks(possible_pks(temp_el));
        end
    end
    
    % store p1pk and p2pk
    temp_p1pk = temp_p1;
    temp_p2pk = temp_p2;
    
    % store points
    for fid_pt_no = 1 : length(fid_pt_names)
        eval(['curr_temp_el = temp_' fid_pt_names{fid_pt_no} ';']);
        if ~isnan(curr_temp_el)
            eval(['pts.' fid_pt_names{fid_pt_no} '(pulse_no) = curr_temp_el + curr_els(1)-1;'])
        else
            eval(['pts.' fid_pt_names{fid_pt_no} '(pulse_no) = nan;'])
        end
    end
    
    
    clear curr temp* curr_els empty_log
    
end
clear pulse_no


%% Ensure only pts are provided for those pulses with all pts available
pt_names = fieldnames(pts);
include_pulse = true(size(pts.dia));
for pt_name_no = 1 : length(pt_names)
    eval(['curr_pt_measures = pts.' pt_names{pt_name_no} ';']);
    include_pulse(isnan(curr_pt_measures)) = false;    
end
for pt_name_no = 1 : length(pt_names)
    if strcmp(pt_names{pt_name_no}, 'f1') || strcmp(pt_names{pt_name_no}, 'f2')
        continue
    end
    eval(['pts.' pt_names{pt_name_no} '(~include_pulse) = nan;']);
end

end

function temp_s = identify_s(curr)

[~, temp_s] = max(curr.sig);

end

function temp_a = identify_a(curr, initial_buffer)

[~,filtered_ms] = max(curr.derivs.first);

pks = find_pks_trs(curr.derivs.second, 'pk');
rel_pks = pks(pks > initial_buffer & pks < filtered_ms); 
if isempty(rel_pks) && sum(pks<=initial_buffer) > 0   % Added in PulseAnalyse5
    rel_pks = pks(find(pks <= initial_buffer, 1, 'last'));
end
[~, temp_el] = max(curr.derivs.second(rel_pks));
temp_a = rel_pks(temp_el);

end

function temp_e = identify_e(curr, temp_s, temp_ms, temp_b, end_buffer, redo_log)

% Find local maxima in the second derivative
pks = find_pks_trs(curr.derivs.second, 'pk');
% Set an upper bound of 60 % of the PW duration
upper_bound = 0.6*length(curr.sig);   % const from: https://en.wikipedia.org/wiki/QT_interval#/media/File:QT_interval_corrected_for_heart_rate.png
% Set a lower bound of 'ms'
lower_bound = temp_ms;
% Identify the highest local maximum in the second derivative between these two bounds
rel_pks = pks(pks >= lower_bound & pks <= upper_bound);
[~, max_el] = max(curr.derivs.second(rel_pks));
% If this is the first local maximum in this search region ...
if max_el == 1
    % work out whether this has detected the "c" wave
    % - find out if there's an inflection point between "b" and this
    temp_trs = find_pks_trs(curr.derivs.third, 'tr');
    no_infl = sum(temp_trs > temp_b & temp_trs < rel_pks(max_el));
    % - if not then take the next peak
    if no_infl == 0
        % if there is 1 peak in this search region ...
        if length(rel_pks) < max_el+1   % Added in PulseAnalyse5
            % then take the next peak (slightly above the upper bound
            orig_el = find(pks >= lower_bound & pks <= upper_bound);
            rel_pks = pks(orig_el:orig_el+1);
        end
        rel_pk = rel_pks(max_el+1);
    else
        rel_pk = rel_pks(max_el);
    end
else
    rel_pk = rel_pks(max_el);
end
temp_e = rel_pk;

end

function temp_f = identify_f(curr, temp_e, end_buffer)

lower_bound = temp_e;
upper_bound = end_buffer;
trs = find_pks_trs(curr.derivs.second, 'tr');
possible_els = trs(trs >= lower_bound & trs <= upper_bound);

if isempty(possible_els)    % Added in PulseAnalyse5
    possible_els = trs(find(trs >=lower_bound, 1));
end

if isempty(possible_els)
    temp_f = nan;
else
    temp_f = possible_els(1);
end

end

function temp_b = identify_b(curr, temp_a)

% find b (PulseAnalyse6 onwards)

% Find local minima in second derivative
trs = find_pks_trs(curr.derivs.second, 'tr');
% define an upper bound as 25% of the duration of the signal
upper_bound = 0.25*length(curr.sig);
% find local minima between 'a' and this upper bound
temp = find(trs > temp_a & curr.derivs.second(trs) < 0 & trs < upper_bound);
% Identify the lowest of these local minima
[~, rel_el] = min(curr.derivs.second(trs(temp)));
temp = temp(rel_el);
temp_b = trs(temp); clear temp

end

function temp_d = identify_d(curr, temp_c, temp_e)

% Identify "d" as the lowest minimum of the second deriv between "c" and "e"
trs = find_pks_trs(curr.derivs.second, 'tr');
possible_trs = find(trs > temp_c & trs < temp_e);
if ~isempty(possible_trs)
    temp = trs(possible_trs);
    [~, temp_el] = min(curr.derivs.second(temp));
    temp_d = temp(temp_el); clear temp
else
    % unless there isn't a minimum, in which case it's an inflection, and
    % "d" is the same as "c"
    temp_d = temp_c;
end

end

function temp_c = identify_c(curr, temp_b, temp_e)

% Identify C as the highest local maximum on the second derivative between "b" and "e"
pks = find_pks_trs(curr.derivs.second, 'pk');
temp = find(pks > temp_b & pks < temp_e);
[~, rel_tr_el] = max(curr.derivs.second(pks(temp)));
temp_c = pks(temp(rel_tr_el)); clear temp rel_tr_el pks

% If there aren't any peaks that satisfy this criterion ...
if isempty(temp_c)
    % then identify C as the lowest local minimum on the third derivative
    % after "b" and before "e"
    trs = find_pks_trs(curr.derivs.third, 'tr');
    temp = find(trs > temp_b & trs < temp_e);
    [~, rel_tr_el] = min(curr.derivs.third(trs(temp)));
    if ~isempty(rel_tr_el)
        temp_c = trs(temp(rel_tr_el)); clear temp rel_tr_el trs
    end
end

end

function temp_dic = identify_dic(curr,temp_e)

temp_dic = temp_e;

end

function temp_dia = identify_dia(curr, temp_dic, temp_e, end_buffer)

% if there is a diastolic peak, then use that:
%  -  first peak in signal after "dic"
pks = find_pks_trs(curr.sig, 'pks');
temp_dia = pks(find(pks > temp_dic & pks < end_buffer, 1));

% if not, then ...
% % I tried (i) take first peak on first derivative after "e"
if isempty(temp_dia)
    pks = find_pks_trs(curr.derivs.first, 'pks');
    temp_dia = pks(find(pks > temp_e & pks < end_buffer, 1));
end
% % But the problem is that the first derivative isn't necessarily a peak at
% % the diastolic peak - it can be an inflection point. So:
% (ii) the soonest out of (a) first peak on first derivative after "e"
%                         (b) first min on third derivative after "e"
% if isempty(temp_dia)
%     pks = find_pks_trs(curr.derivs.first, 'pks');
%     temp_dia1 = pks(find(pks > temp_e, 1));
%     trs = find_pks_trs(curr.derivs.third, 'trs');
%     temp_dia2 = trs(find(trs > temp_e, 1));
%     temp_dia = min([temp_dia1, temp_dia2]);
% end


end

function temp_ms = identify_ms(curr)

% find max slope in DPPG
[~, temp_ms] = max(curr.derivs.first);

end

function temp_p1 = identify_p1(curr, temp_b, temp_ms, buffer_p1, no_of_pulses, fs, up)

% find p1

% find local minima in the first derivative
fd_trs = find_pks_trs(curr.derivs.first, 'tr');
% find local maxima in the second derivative
sd_pks = find_pks_trs(curr.derivs.second, 'pk');

% Find the first local minimum in the first derivative after 0.1 s
current_buffer = buffer_p1(1);
temp = find(fd_trs > current_buffer,1);
% Find the second local minimum (or failing that, the first) in the first derivative after 'b'
temp2 = find(fd_trs > temp_b, 2);
if length(temp2) > 1
    temp2 = temp2(2);
end
% Take whichever comes first:
if temp2 < temp
    temp = temp2;
end
temp_p1 = fd_trs(temp);

% If this value for p1 is after the buffer of 0.18 s ...
if temp_p1 > buffer_p1(2)
    curr.derivs.fourth = savitzky_golay(curr.derivs.third, 1, 9);
    % Then find the first local minimum in the fourth derivative after 0.1 s
    fd_trs = find_pks_trs(curr.derivs.fourth, 'tr');
    temp = find(fd_trs > current_buffer,1);
    temp_p1 = fd_trs(temp); clear temp
end

% If this value for p1 is after the buffer of 0.18 s ...
if temp_p1 > buffer_p1(2)
    % Then find the last local minimum in the first derivative before 0.18 s
    temp_p1 = fd_trs(find(fd_trs <= current_buffer,1,'last'));
end

end

function temp_p2 = identify_p2(curr, temp_d, temp_p1, temp_dic)

% Find "p2" from the minimum value of the third derivative immediately before "d"
td_trs = find_pks_trs(curr.derivs.third, 'tr');
temp = find(td_trs < temp_d,1,'last'); clear d_pks
temp_p2 = td_trs(temp);

% unless c=d, in which case p2 will now be before p1, so take the minimum
% value of the third derivative immediately after "d"
if temp_p2 < temp_p1
    temp_p2 = td_trs(find(td_trs<temp_dic,1,'last'));
end

% check to see if there is a maximum in the signal between the current
% estimate of p2, and temp_dic. If so, take this instead
pks = find_pks_trs(curr.sig, 'pk');
temp = find(pks> temp_p2 & pks < temp_dic);
if length(temp) == 1
    temp_p2 = pks(temp);
elseif length(temp) == 2
    temp_p2 = pks(temp(2));
elseif length(temp) > 1
    fprintf('\nCheck this')
end
clear pks

end

function pts = annotate_fiducial_pts(sig, onsets, derivs, up, no_of_pulses)
% IDENTIFY_FIDUCIAL_PTS  Used to annotate fiducial points on each wave of a
% pulsatile signal.
%
% Inputs:
%       sig         -  pulsatile signal, a structure containing .v (a
%                       vector of values), and .fs (sampling frequency in Hz).
%       onsets      -  indices of the pulse onsets
%       derivs      -  a structure containing .first, .second and .third
%                       derivatives of sig
%
% Outputs:
%       pts         -  a structure containing the indices of the fiducial
%                       points for each beat (a nan indicates that a fiducial
%                       point could not be identified).
%


%% annotate fiducial points

[pts.a, pts.b, pts.c, pts.d, pts.e, pts.f, pts.s, pts.dia, pts.dic, pts.p1, pts.p2, pts.ms, pts.f1, pts.f2] = deal(nan(length(onsets)-1,1));
switch no_of_pulses
    case 'multiple'
        start_pulse = 100;
        end_pulse = 100;
    case 'single'
        start_pulse = 1;
        end_pulse = length(onsets)-1;
end
for pulse_no = start_pulse : end_pulse
    
    %% extract data for this pulse
    curr_els = onsets(pulse_no):onsets(pulse_no+1);
    curr.sig.v = sig.v(curr_els);
    curr.derivs.first = derivs.first(curr_els);
    curr.derivs.second = derivs.second(curr_els);
    curr.derivs.third = derivs.third(curr_els);
    
    %% Plot waves
    
    %- setup
    paper_size = [600,1050];
    figure('Position', [50,50, paper_size])
    ftsize = 16; lwidth = 2;
    pp_annot = 0.8;
    rt_annot = -0.02;
    curr.sig.t = [0:length(curr.sig.v)-1]/sig.fs;
    
    h_b(1) = subplot('Position', [0.21,0.78,0.78,0.21]);
    
    % plot baseline curve
    plot(curr.sig.t, curr.sig.v, 'LineWidth', lwidth); hold on,
    
    % set limits
    curr_range = range(curr.sig.t);
    xlim([curr.sig.t(1)-0.1*curr_range, curr.sig.t(end)+0.1*curr_range])
    
    % set labels
    ylab = ylabel('sig', 'FontSize', ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.18, 0.5, 0]);
    set(gca, 'FontSize', ftsize -2)
    set(gca, 'XTick', 0:0.25:1, 'XTickLabel', {})
    box off
    
    %- plot first derivative
    h_b(2) = subplot('Position', [0.21,0.54,0.78,0.21]);
    
    % Plot x-axis
    plot([-10, 10], [0,0], '--k'); hold on
    
    % plot baseline curve
    plot(curr.sig.t, curr.derivs.first, 'LineWidth', lwidth); hold on,
    
    % set limits
    curr_range = range(curr.sig.t);
    xlim([curr.sig.t(1)-0.1*curr_range, curr.sig.t(end)+0.1*curr_range])
    ylims = ylim; curr_range = range(curr.derivs.first);
    ylim([min(curr.derivs.first)-0.05*curr_range, max(curr.derivs.first)+0.15*curr_range])
    
    % set labels
    ylab = ylabel('DPPG', 'FontSize', ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.18, 0.5, 0]);
    set(gca, 'FontSize', ftsize -2)
    set(gca, 'XTick', 0:0.25:1, 'XTickLabel', {})
    box off
    
    %- plot Second derivative
    h_b(3) = subplot('Position', [0.21,0.30,0.78,0.21]);
    
    % Plot x-axis
    plot([-10, 10], [0,0], '--k'); hold on
    
    % plot baseline curve
    curr_color = 0.4*0.5*[1,2,1];
    plot(curr.sig.t, curr.derivs.second, 'LineWidth', lwidth); hold on,
    
    % set limits
    curr_range = range(curr.sig.t);
    xlim([curr.sig.t(1)-0.1*curr_range, curr.sig.t(end)+0.1*curr_range])
    ylims = ylim; curr_range = range(curr.derivs.second);
    ylim([min(curr.derivs.second)-0.05*curr_range, max(curr.derivs.second)+0.15*curr_range])
    
    % set labels
    ylab = ylabel('2nd D', 'FontSize', ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.18, 0.5, 0]);
    set(gca, 'FontSize', ftsize -2, 'XTick', 0:0.25:1, 'XTickLabel', {})
    box off
    
    %- plot Third derivative
    h_b(4) = subplot('Position', [0.21,0.06,0.78,0.21]);
    
    % Plot x-axis
    plot([-10, 10], [0,0], '--k'); hold on
    
    % plot baseline curve
    curr_color = 0.4*0.5*[1,2,1];
    plot(curr.sig.t, curr.derivs.third, 'LineWidth', lwidth); hold on,
    
    % set limits
    curr_range = range(curr.sig.t);
    xlim([curr.sig.t(1)-0.1*curr_range, curr.sig.t(end)+0.1*curr_range])
    ylims = ylim; curr_range = range(curr.derivs.third);
    ylim([min(curr.derivs.third)-0.05*curr_range, max(curr.derivs.third)+0.15*curr_range])
    
    % set labels
    xlabel('Time [s]', 'FontSize', ftsize)
    ylab = ylabel('3rd D', 'FontSize', ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.18, 0.5, 0]);
    set(gca, 'FontSize', ftsize -2, 'XTick', 0:0.25:1)
    box off
    
    linkaxes(h_b, 'x')
    shg
    
    %% Annotation
    
    happy_with_annotations = false;
    while ~happy_with_annotations
        pts.a = 4;
        pts.b = 13;
        pts.c = 22;
        pts.d = 35;
        pts.e = 44;
        pts.f = 53;
        pts.s = 33;
        pts.dia = 49;
        pts.dic = 45;
        pts.p1 = 17;
        pts.p2 = 33;
        pts.ms = 5;
        pts.f1 = 1;
        pts.f2 = 57;
%         pt_names = fieldnames(pts);
%         for pt_no = 1 : length(pt_names)
%             % Update title
%             axes(h_b(1))
%             title(['Annotate pt ' pt_names{pt_no}])
%             % Perform Annotation
%             [temp,~] = ginput(1);
%             % Find closest point to annotated time
%             [~, temp2] = min(abs(curr.sig.t-temp));
%             % Store this annotation
%             eval(['pts.' pt_names{pt_no} '(pulse_no) = temp2 + curr_els(1)-1;']);
%             clear temp temp2
%         end
        
        %% Check that you are happy with these points
        
        % original pulse signal
        axes(h_b(1))
        rel_pts = {'s', 'dia', 'dic', 'f1', 'f2', 'p1', 'p2'};
        vspace0 = 0.08*range(ylim);
        for pt_name_no = 1 : length(rel_pts)
            eval(['curr_el = pts.' rel_pts{pt_name_no} ';'])
            plot(curr.sig.t(curr_el), curr.sig.v(curr_el), 'or')
            curr_text = rel_pts{pt_name_no};
            hspace0 = 0;ftsize = 12;
            text(curr.sig.t(curr_el)+hspace0, curr.sig.v(curr_el)+vspace0 , curr_text,'FontSize', ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
        end
        
        % first derivative
        axes(h_b(2))
        rel_pts = {'ms'};
        vspace0 = 0.08*range(ylim);
        for pt_name_no = 1 : length(rel_pts)
            eval(['curr_el = pts.' rel_pts{pt_name_no} ';'])
            plot(curr.sig.t(curr_el), curr.derivs.first(curr_el), 'or')
            curr_text = rel_pts{pt_name_no};
            hspace0 = 0;ftsize = 12;
            text(curr.sig.t(curr_el)+hspace0, curr.derivs.first(curr_el)+vspace0 , curr_text,'FontSize', ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
        end
        
        % second derivative
        axes(h_b(3))
        rel_pts = {'a', 'b', 'c', 'd', 'e', 'f'};
        vspace0 = 0.08*range(ylim);
        for pt_name_no = 1 : length(rel_pts)
            eval(['curr_el = pts.' rel_pts{pt_name_no} ';'])
            plot(curr.sig.t(curr_el), curr.derivs.second(curr_el), 'or')
            curr_text = rel_pts{pt_name_no};
            hspace0 = 0;ftsize = 12;
            text(curr.sig.t(curr_el)+hspace0, curr.derivs.second(curr_el)+vspace0 , curr_text,'FontSize', ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
        end
        
        % third derivative
        axes(h_b(4))
        rel_pts = {'a', 'b', 'c', 'd', 'e', 'f'};
        vspace0 = 0.08*range(ylim);
        for pt_name_no = 1 : length(rel_pts)
            eval(['curr_el = pts.' rel_pts{pt_name_no} ';'])
            plot(curr.sig.t(curr_el), curr.derivs.third(curr_el), 'or')
            curr_text = rel_pts{pt_name_no};
            hspace0 = 0;ftsize = 12;
            text(curr.sig.t(curr_el)+hspace0, curr.derivs.third(curr_el)+vspace0 , curr_text,'FontSize', ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
        end
        
        choice = questdlg('Are these annotations ok?', ...
            'Yes','No');
        % Handle response
        switch choice
            case 'Yes'
                happy_with_annotations = true;
            case 'No'
            case 'Cancel'
                error('User stopped script')
        end
        
    end
    close all
    clear curr temp* trs curr_els happy_with_annotations
    
end
clear pulse_no

end

function buttonDownCallback(o,e)
p = get(gca,'CurrentPoint');
p = p(1,1:2);
title( sprintf('(%g,%g)',p) )
end

function  Click_CallBack2(a)

switch get(ancestor(a,'figure'),'SelectionType')
    
    case 'normal' %left click
        point = get(a,'CurrentPoint');
        load(annotations_filepath);
        [pk_anns.t,inds]=unique([pk_anns.t; point(1,1)+start_time]);
        pk_anns.v=[pk_anns.v; point(1,2)]; pk_anns.v = pk_anns.v(inds);
        [pk_anns.t, inds] = sort(pk_anns.t);
        pk_anns.v = pk_anns.v(inds);
        save(annotations_filepath, 'pk_anns');
        
    case 'alt'  % right click
        point = get(a,'CurrentPoint');
        load(annotations_filepath);
        [pk_anns.t,inds]=unique([pk_anns.t; point(1,1)]);
        pk_anns.v=[pk_anns.v; point(1,2)]; pk_anns.v = pk_anns.v(inds);
        [pk_anns.t, inds] = sort(pk_anns.t);
        pk_anns.v = pk_anns.v(inds);
        save(annotations_filepath, 'pk_anns');
        
    case 'extend' % right click whilst holding down shift
        load(annotations_filepath)
        pk_anns.t = pk_anns.t(1:(end-1));
        pk_anns.v = pk_anns.v(1:(end-1));
        save(annotations_filepath, 'pk_anns');
        
end

cla(axis_h)
plot(pk_anns.t-start_time,pk_anns.v, 'ro','LineWidth',4)

end

function trs = find_pks_trs(sig,type)

if strcmp(type, 'tr')
    sig = -sig;
end

temp1 = find(sig(2:end-1) > sig(1:(end-2)) & sig(2:end-1) > sig(3:end) );

temp2 = find(sig(2:end-2) > sig(1:(end-3)) & sig(2:end-2)== sig(3:(end-1)) & sig(3:end-1) > sig(4:end) );

temp = unique([temp1; temp2]);

trs = temp+1;

end

function stiff_inds = calc_stiffness_inds(sigs, pulses, pts, options)
% CALC_STIFFNESS_INDS  Calculates stiffness indices from a pulsatile
% signal.
%
% Inputs:
%       sig         -  pulsatile signal, a structure containing .v (a
%                       vector of values), and .fs (sampling frequency in Hz).
%       derivs      -  a structure containing .first, .second and .third
%                       derivatives of sig
%       pulses      -  a structure containing pulses.quality, a logical
%                       indicating the signal quality of each pulse (true
%                       indicates high quality).
%       pts         -  a structure containing the indices of the fiducial
%                       points for each beat.
%
% Outputs:
%       stiff_inds  -  a structure containing the stiffness indices for
%                       each beat (a nan indicates that the stiffness index for that beat
%                       could not be identified).
%

%% Setup

% create time vector
sigs.t = [0:length(sigs.s)-1]'/sigs.fs;

% find out whether a height measurement has been provided
if sum(strcmp(fieldnames(sigs), 'ht'))
    ht_provided = 1;
else
    ht_provided = 0;
end

% calculate additional features
T = nan(size(pts.f2));
good_els = ~isnan(pts.s);
T(good_els) = sigs.t(pts.f2(good_els)) - sigs.t(pts.f1(good_els));

% make vectors for each stiffness index (in case there are nans in the pts indices)
si_names = {'delta_t', 'CT', 'SI', 'CT_div_h', 'prop_s', 't_sys', 't_dia', 't_ratio', 'prop_delta_t', 't_p1_dia', 't_p2_dia', 'AI', 'AP', 'RI', 'RI_p1', 'RI_p2', 'ratio_p2_p1', 'A1', 'A2', 'IPA', 'ms', 'ms_div_amp', 'b_div_a', 'c_div_a', 'd_div_a', 'e_div_a', 'a_div_amp', 'b_div_amp', 'c_div_amp', 'd_div_amp', 'e_div_amp', 'AGI', 'AGI_inf', 'AGI_mod', 't_b_c', 't_b_d', 'slope_b_c', 'slope_b_d', 'IPAD', 'k'};
for name_no = 1 : length(si_names)
    eval(['stiff_inds.' si_names{name_no} ' = nan(length(pts.f1),1);']);
end

%% Find values of each stiffness index for each beat

% calculate SIs which can be found from existing pts on pulsatile signal (timings)
stiff_inds.delta_t(good_els) = sigs.t(pts.dia(good_els)) - sigs.t(pts.s(good_els));
stiff_inds.CT(good_els) = sigs.t(pts.s(good_els)) - sigs.t(pts.f1(good_els));
if ht_provided
    stiff_inds.SI(good_els) = sigs.ht./stiff_inds.delta_t(good_els);
    stiff_inds.CT_div_h(good_els) = stiff_inds.CT(good_els)./sigs.ht;
end
stiff_inds.prop_s(good_els) = stiff_inds.CT(good_els)./T(good_els);
stiff_inds.t_sys(good_els) = sigs.t(pts.dic(good_els)) - sigs.t(pts.f1(good_els));
stiff_inds.t_dia(good_els) = sigs.t(pts.f2(good_els)) - sigs.t(pts.dic(good_els));
stiff_inds.t_ratio(good_els) = (sigs.t(pts.s(good_els)) - sigs.t(pts.f1(good_els)))./stiff_inds.t_sys(good_els);
stiff_inds.prop_delta_t(good_els) = stiff_inds.delta_t(good_els)./T(good_els);
stiff_inds.t_p1_dia(good_els) = sigs.t(pts.dia(good_els)) - sigs.t(pts.p1pk(good_els));
stiff_inds.t_p2_dia(good_els) = sigs.t(pts.dia(good_els)) - sigs.t(pts.p2pk(good_els));

% calculate SIs which can be found from existing pts on pulsatile signal (amplitudes)
stiff_inds.AI(good_els) = 100*(sigs.s(pts.p2pk(good_els)) - sigs.s(pts.p1in(good_els)))./ (sigs.s(pts.s(good_els)) - sigs.s(pts.f1(good_els)));
stiff_inds.AP(good_els) = sigs.s(pts.p2pk(good_els)) - sigs.s(pts.p1in(good_els));
stiff_inds.RI(good_els) = (sigs.s(pts.dia(good_els)) - sigs.s(pts.f1(good_els))) ./ (sigs.s(pts.s(good_els)) - sigs.s(pts.f1(good_els)));
stiff_inds.RI_p1(good_els) = (sigs.s(pts.dia(good_els)) - sigs.s(pts.f1(good_els))) ./ (sigs.s(pts.p1pk(good_els)) - sigs.s(pts.f1(good_els)));
stiff_inds.RI_p2(good_els) = (sigs.s(pts.dia(good_els)) - sigs.s(pts.f1(good_els))) ./ (sigs.s(pts.p2pk(good_els)) - sigs.s(pts.f1(good_els)));
stiff_inds.ratio_p2_p1(good_els) = (sigs.s(pts.p2pk(good_els)) - sigs.s(pts.f1(good_els))) ./ (sigs.s(pts.p1in(good_els)) - sigs.s(pts.f1(good_els)));

% calculate additional SIs which require additional pulse wave analysis
for beat_no = 1:length(pts.f1)
    
    if isnan(pts.s(beat_no))
        continue
    end
    
    % find pulse amplitude
    curr_amp = sigs.s(pts.s(beat_no)) - sigs.s(pts.f1(beat_no));
    
    % find areas
    baseline = linspace(sigs.s(pts.f1(beat_no)), sigs.s(pts.f2(beat_no)), pts.f2(beat_no) - pts.f1(beat_no)+1); baseline = baseline(:);
    rel_pts = pts.f1(beat_no) : pts.dic(beat_no);
    baseline_pts = rel_pts - pts.f1(beat_no) + 1;
    stiff_inds.A1(beat_no) = sum(sigs.s(rel_pts) - baseline(baseline_pts))/( sigs.fs*curr_amp);
    rel_pts = pts.dic(beat_no) : pts.f2(beat_no);
    baseline_pts = rel_pts - pts.f1(beat_no) + 1;
    stiff_inds.A2(beat_no) = sum(sigs.s(rel_pts) - baseline(baseline_pts))/(sigs.fs*curr_amp);

end
stiff_inds.IPA(good_els) = stiff_inds.A2(good_els) ./ stiff_inds.A1(good_els);
clear beat_no rel_beats

% calculate SIs from first derivative (amplitudes)
if sum(good_els)>0
    stiff_inds.ms(good_els) = sigs.derivs.first(pts.ms(good_els));
    stiff_inds.ms_div_amp(good_els) = sigs.derivs.first(pts.ms(good_els)) ./ (sigs.s(pts.s(good_els)) - sigs.s(pts.f1(good_els)));
end

% Calculate SIs from second derivative (amplitudes)
stiff_inds.b_div_a(good_els) = sigs.derivs.second(pts.b(good_els)) ./ sigs.derivs.second(pts.a(good_els));
stiff_inds.c_div_a(good_els) = sigs.derivs.second(pts.c(good_els)) ./ sigs.derivs.second(pts.a(good_els));
stiff_inds.d_div_a(good_els) = sigs.derivs.second(pts.d(good_els)) ./ sigs.derivs.second(pts.a(good_els));
stiff_inds.e_div_a(good_els) = sigs.derivs.second(pts.e(good_els)) ./ sigs.derivs.second(pts.a(good_els));
stiff_inds.a_div_amp(good_els) = sigs.derivs.second(pts.a(good_els)) ./ (sigs.s(pts.s(good_els)) - sigs.s(pts.f1(good_els)));
stiff_inds.b_div_amp(good_els) = sigs.derivs.second(pts.b(good_els)) ./ (sigs.s(pts.s(good_els)) - sigs.s(pts.f1(good_els)));
stiff_inds.c_div_amp(good_els) = sigs.derivs.second(pts.c(good_els)) ./ (sigs.s(pts.s(good_els)) - sigs.s(pts.f1(good_els)));
stiff_inds.d_div_amp(good_els) = sigs.derivs.second(pts.d(good_els)) ./ (sigs.s(pts.s(good_els)) - sigs.s(pts.f1(good_els)));
stiff_inds.e_div_amp(good_els) = sigs.derivs.second(pts.e(good_els)) ./ (sigs.s(pts.s(good_els)) - sigs.s(pts.f1(good_els)));

stiff_inds.AGI(good_els) = ( sigs.derivs.second(pts.b(good_els)) - sigs.derivs.second(pts.c(good_els)) - sigs.derivs.second(pts.d(good_els)) - sigs.derivs.second(pts.e(good_els)) )./sigs.derivs.second(pts.a(good_els));
stiff_inds.AGI_inf(good_els) = ( sigs.derivs.second(pts.b(good_els)) - sigs.derivs.second(pts.e(good_els)) )./sigs.derivs.second(pts.a(good_els));
stiff_inds.AGI_mod(good_els) = ( sigs.derivs.second(pts.b(good_els)) - sigs.derivs.second(pts.c(good_els)) - sigs.derivs.second(pts.d(good_els)) )./sigs.derivs.second(pts.a(good_els));

stiff_inds.t_b_c(good_els) = sigs.t(pts.c(good_els)) - sigs.t(pts.b(good_els));
stiff_inds.t_b_d(good_els) = sigs.t(pts.d(good_els)) - sigs.t(pts.b(good_els));

% Calculate SIs from second derivative (slopes)
pt1.t = sigs.t(pts.b(good_els));
pt1.v = sigs.derivs.second(pts.b(good_els));
pt2.t = sigs.t(pts.c(good_els));
pt2.v = sigs.derivs.second(pts.c(good_els));
stiff_inds.slope_b_c(good_els) = ((pt2.v - pt1.v)./sigs.derivs.second(pts.a(good_els)))./(pt2.t-pt1.t);
pt2.t = sigs.t(pts.d(good_els));
pt2.v = sigs.derivs.second(pts.d(good_els));
stiff_inds.slope_b_d(good_els) = ((pt2.v - pt1.v)./sigs.derivs.second(pts.a(good_els)))./(pt2.t-pt1.t);

% combined
stiff_inds.IPAD(good_els) = stiff_inds.IPA(good_els) + (sigs.derivs.second(pts.d(good_els)) ./ sigs.derivs.second(pts.a(good_els)));
if sum(good_els)>0
    stiff_inds.k(good_els) = (sigs.derivs.second(pts.s(good_els))) ./ (( sigs.s(pts.s(good_els)) - sigs.s(pts.ms(good_els)) ) ./ (sigs.s(pts.s(good_els)) - sigs.s(pts.f1(good_els))) );
end

%% Calculate median values of each stiffness index

stiff_ind_names = fieldnames(stiff_inds);

% if selected, then calculate the median CV indices using only high quality beats
if options.exclude_low_quality_data
    rel_beats = pulses.quality;
else
    rel_beats = true(size(pulses.onsets));
end
rel_beats = rel_beats(1:end-1);

% Find median values of each index
for stiff_ind_no = 1 : length(stiff_ind_names)
        
    % Find median value if there are multiple beats
    if length(pts.f1) > 1
        eval(['curr_stiff_ind_data.raw = stiff_inds.' stiff_ind_names{stiff_ind_no} ';'])
        curr_stiff_ind_data.v = nanmedian(curr_stiff_ind_data.raw(rel_beats));
    else
        eval(['curr_stiff_ind_data.v = stiff_inds.' stiff_ind_names{stiff_ind_no} ';'])
    end
    
    % store result
    eval(['stiff_inds.' stiff_ind_names{stiff_ind_no} ' = curr_stiff_ind_data;'])
    
end
clear curr_stiff_ind_data stiff_ind_no rel_beats stiff_ind_names

end

function S_elim_lf = eliminate_low_freq_from_single_beat(sig, up)

% Correct for low frequency baseline drift in a single beat

diff_1 = sig.v(2) - sig.v(1);
desired_val_end = sig.v(1) - diff_1;
correction_line = linspace(0, sig.v(end)-desired_val_end, length(sig.v));
S_elim_lf.v = sig.v - correction_line';
S_elim_lf.fs = sig.fs;

end

function S_aligned = align_pulse(sig, up)

% Ensure that signal commences at start of systole

[~, align_el] = min(sig.v);
S_aligned.v = sig.v([align_el:end, 1:(align_el-1)]);
S_aligned.fs = sig.fs;

% add on one additional point so that you can define a second onset
S_aligned.v(end+1) = S_aligned.v(1);

end

function pulses = generate_pulses_for_single_beat(S_filt)

% Generate "pulses" information for a single-beat input signal

pulses.onsets = [1, length(S_filt.v)];
pulses.quality = [true, false];
[~, pulses.peaks] = max(S_filt.v);
pulses.peaks(end+1) = length(S_filt.v);

end

function plot_fiducial_points(sigs, pulses, fid_pts, no_of_pulses, options)
% make plot of individual beat if needed

sig.v = sigs.s;
sig.fs = sigs.fs;
derivs = sigs.derivs;

%- setup
paper_size = [600,1050];
figure('Position', [50,50, paper_size])
ftsize = 16; lwidth = 2;
pp_annot = 0.8;
rt_annot = -0.02;
sig.t = [0:length(sig.v)-1]/sig.fs;

if options.plot_third_deriv
    y_inc = 0.23;
    n_sub = 4;
    ftsize = ftsize + 2;
else
    y_inc = 0.31;
    n_sub = 3;
end
y_offset = 0.08;

%- plot sig
h_b(1) = subplot('Position', [0.21,y_offset+(n_sub-1)*y_inc,0.78,y_inc-0.01]);

% identify relevant beat
if strcmp(no_of_pulses, 'single')
    beat_no = 1;
elseif length(fid_pts.f2) < 101
    beat_no = 2;
else
    beat_no = 100;
end
rel_els = fid_pts.f1(beat_no) : fid_pts.f2(beat_no);
curr_sig.t = sig.t(rel_els)-sig.t(rel_els(1));
curr_sig.sig = subtract_baseline(sig.v(rel_els));
[curr_sig.sig, scale_factor] = normalise(curr_sig.sig);

% plot baseline curve
plot(curr_sig.t, curr_sig.sig, 'LineWidth', lwidth); hold on,

%         % add dashed lines
%         plot(sig.t(fid_pts(1))*[1,1], [rt_annot, sig.v(fid_pts(1))], '--', 'color', 0.6*ones(1,3))
%         plot(sig.t(fid_pts(2))*[1,1], [rt_annot, sig.v(fid_pts(2))], '--', 'color', 0.6*ones(1,3))
%         plot([sig.t(fid_pts(1)),pp_annot], sig.v(fid_pts(1))*ones(1,2), '--', 'color', 0.6*ones(1,3))
%         plot([sig.t(fid_pts(2)),pp_annot], sig.v(fid_pts(2))*ones(1,2), '--', 'color', 0.6*ones(1,3))

% plot salient points
pt_names = {'dia', 'dic', 's', 'p1', 'p2'};
hspace0 = 0;
for pt_no = 1 : length(pt_names)
    eval(['curr_pt.el = fid_pts.' pt_names{pt_no} '(beat_no) - fid_pts.f1(beat_no)+1;']);
    if isnan(curr_pt.el)
        continue
    end
    
    % plot point
    curr_pt.v = curr_sig.sig(curr_pt.el);
    curr_pt.t = curr_sig.t(curr_pt.el);
    plot(curr_pt.t, curr_pt.v, 'or')
    
    % annotate point
    curr_text = pt_names{pt_no};
    vspace0 = 0.08;
    switch curr_text
        case {'dia','p2'}
            hspace0 = 0.04;
            if (strcmp(curr_text, 'p2') || strcmp(curr_text, 'p2in')) && curr_pt.el < (fid_pts.s(beat_no) - fid_pts.f1(beat_no)+1)
                hspace0 = -0.04;
            elseif (strcmp(curr_text, 'p1') || strcmp(curr_text, 'p1in'))
                hspace0 = 0;
            end
        case 'dic'
            hspace0 = -0.01;
            vspace0 = -0.08;
        case {'f1', 'p1','p1in'}
            hspace0 = -0.04;
        case 's'
            hspace0 = 0;
    end
    text(curr_sig.t(curr_pt.el)+hspace0, curr_sig.sig(curr_pt.el)+vspace0 , curr_text,'FontSize', ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
    
end

% set limits
ylim([-0.08, 1.2])
curr_range = range(curr_sig.t);
xlim([curr_sig.t(1)-0.1*curr_range, curr_sig.t(end)+0.1*curr_range])

% set labels
ylab = ylabel('PW', 'FontSize', ftsize, 'Rotation', 0);
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(gca, 'FontSize', ftsize -2, 'YTick', [])
set(gca, 'XTick', 0:0.25:1, 'XTickLabel', {})
box off

%- plot first derivative
h_b(2) = subplot('Position', [0.21,y_offset+(n_sub-2)*y_inc,0.78,y_inc - 0.01]);

curr_sig.derivs.first = derivs.first(rel_els)/scale_factor;

% Plot x-axis
plot([-10, 10], [0,0], '--k'); hold on

% plot baseline curve
plot(curr_sig.t, curr_sig.derivs.first, 'LineWidth', lwidth); hold on,

% plot salient points
pt_names = {'ms', 'dia'};
curr_range = range(curr_sig.derivs.first);
vspace0 = 0.08*curr_range;
hspace0 = 0.05;
for pt_no = 1 : length(pt_names)
    eval(['curr_fid_pt = fid_pts.' pt_names{pt_no} '(beat_no);'])
    if isnan(curr_fid_pt)
        continue
    end
    curr_pt.el = curr_fid_pt - fid_pts.f1(beat_no)+1;
    if isnan(curr_pt.el)
        continue
    end
    
    % plot point
    curr_pt.v = curr_sig.derivs.first(curr_pt.el);
    curr_pt.t = curr_sig.t(curr_pt.el);
    plot(curr_pt.t, curr_pt.v, 'or')
    
    % annotate point
    curr_text = pt_names{pt_no};
    text(curr_sig.t(curr_pt.el)+hspace0, curr_sig.derivs.first(curr_pt.el)+vspace0 , curr_text,'FontSize', ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
    
end

% set limits
curr_range = range(curr_sig.t);
xlim([curr_sig.t(1)-0.1*curr_range, curr_sig.t(end)+0.1*curr_range])
ylims = ylim; curr_range = range(curr_sig.derivs.first);
ylim([min(curr_sig.derivs.first)-0.05*curr_range, max(curr_sig.derivs.first)+0.15*curr_range])
%ylim([-0.05 0.15])

% set labels
ylab = ylabel({'1st','derivative'}, 'FontSize', ftsize, 'Rotation', 0);
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(gca, 'FontSize', ftsize -2, 'YTick', [])
set(gca, 'XTick', 0:0.25:1, 'XTickLabel', {})
set(gca, 'YTick', 0, 'YTickLabel', {'0'});
box off


%- plot Second derivative
h_b(3) = subplot('Position', [0.21,y_offset+(n_sub-3)*y_inc,0.78,y_inc - 0.01]);

curr_sig.derivs.second = derivs.second(rel_els)/scale_factor;

% Plot x-axis
plot([-10, 10], [0,0], '--k'); hold on

% plot baseline curve
curr_color = 0.4*0.5*[1,2,1];
plot(curr_sig.t, curr_sig.derivs.second, 'LineWidth', lwidth); hold on,

% plot salient points
pt_names = {'a', 'b', 'c', 'd', 'e', 'f'};
curr_range = range(curr_sig.derivs.second);
vspace_const = 0.08;
hspace0 = 0;
for pt_no = 1 : length(pt_names)
    eval(['curr_pt.el = fid_pts.' pt_names{pt_no} '(beat_no) - fid_pts.f1(beat_no)+1;']);
    if isnan(curr_pt.el)
        continue
    end
    
    % plot point
    curr_pt.v = curr_sig.derivs.second(curr_pt.el);
    curr_pt.t = curr_sig.t(curr_pt.el);
    plot(curr_pt.t, curr_pt.v, 'or')
    
    % annotate point
    curr_text = pt_names{pt_no};
    switch curr_text
        case {'a','c', 'e'}
            vspace0 = (1.3*vspace_const)*curr_range;
        case {'b', 'd', 'f'}
            vspace0 = -1*vspace_const*curr_range;
    end
    text(curr_sig.t(curr_pt.el), curr_sig.derivs.second(curr_pt.el)+vspace0 , curr_text,'FontSize', ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
    
end

% set limits
curr_range = range(curr_sig.t);
xlim([curr_sig.t(1)-0.1*curr_range, curr_sig.t(end)+0.1*curr_range])
ylims = ylim; curr_range = range(curr_sig.derivs.second);
ylim([min(curr_sig.derivs.second)-0.05*curr_range, max(curr_sig.derivs.second)+0.15*curr_range])
%ylim([-0.025 0.025])

% set labels
ylab = ylabel({'2nd','derivative'}, 'FontSize', ftsize, 'Rotation', 0);
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(gca, 'FontSize', ftsize -2, 'XTick', 0:0.25:1, 'YTick', [])
box off
if ~options.plot_third_deriv
    xlabel('Time [s]', 'FontSize', ftsize)
else
    set(gca, 'XTickLabel', {})
end
set(gca, 'YTick', 0, 'YTickLabel', {'0'});

if options.plot_third_deriv
%- plot Third derivative
h_b(4) = subplot('Position', [0.21,y_offset+(n_sub-4)*y_inc,0.78,y_inc-0.01]);

curr_sig.derivs.third = derivs.third(rel_els)/scale_factor;

% Plot x-axis
plot([-10, 10], [0,0], '--k'); hold on

% plot baseline curve
curr_color = 0.4*0.5*[1,2,1];
plot(curr_sig.t, curr_sig.derivs.third, 'LineWidth', lwidth); hold on,

% plot salient points
pt_names = {'p1', 'p2'};
curr_range = range(curr_sig.derivs.third);
vspace_const = 0.08;
hspace0 = 0;
for pt_no = 1 : length(pt_names)
    eval(['curr_pt.el = fid_pts.' pt_names{pt_no} '(beat_no) - fid_pts.f1(beat_no)+1;']);
    if isnan(curr_pt.el)
        continue
    end
    
    % plot point
    curr_pt.v = curr_sig.derivs.third(curr_pt.el);
    curr_pt.t = curr_sig.t(curr_pt.el);
    plot(curr_pt.t, curr_pt.v, 'or')
    
    % annotate point
    curr_text = pt_names{pt_no};
    switch curr_text
        case {'p1'}
            vspace0 = vspace_const*curr_range;
        case {'p2'}
            vspace0 = -1*vspace_const*curr_range;
    end
    text(curr_sig.t(curr_pt.el), curr_sig.derivs.third(curr_pt.el)+vspace0 , curr_text,'FontSize', ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
    
end

% set limits
curr_range = range(curr_sig.t);
xlim([curr_sig.t(1)-0.1*curr_range, curr_sig.t(end)+0.1*curr_range])
ylims = ylim; curr_range = range(curr_sig.derivs.third);
ylim([min(curr_sig.derivs.third)-0.05*curr_range, max(curr_sig.derivs.third)+0.15*curr_range])
%ylim([-0.025 0.025])

% set labels
xlabel('Time [s]', 'FontSize', ftsize)
ylab = ylabel({'3rd','derivative'}, 'FontSize', ftsize, 'Rotation', 0);
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(gca, 'FontSize', ftsize -2, 'XTick', 0:0.25:1, 'YTick', [])
box off
end


linkaxes(h_b, 'x')
shg

if ~isempty(options.save_folder)
    savepath = [options.save_folder, options.save_file, 'fid_pts'];
    set(gcf,'color','w');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize', [paper_size(1), paper_size(2)]./40);
    set(gcf,'PaperPosition',[0 0 paper_size(1) paper_size(2)]./40);
    print(gcf,'-dpdf',savepath)
end
shg

end

function plot_cv_inds(sigs, pulses, fid_pts, cv_inds, no_of_pulses, options)
% make plot of individual beat


sig.v = sigs.s;
sig.fs = sigs.fs;
derivs = sigs.derivs;

%- setup figure
paper_size = [600,1050];
figure('Position', [550,50, paper_size])
fig_props.ftsize = 16; fig_props.lwidth = 2;
pp_annot = 0.8;
rt_annot = -0.02;
sig.t = [0:length(sig.v)-1]/sig.fs;

if options.plot_third_deriv
    fig_props.y_inc = 0.23;
    fig_props.n_sub = 4;
    fig_props.ftsize = fig_props.ftsize + 2;
else
    fig_props.y_inc = 0.31;
    fig_props.n_sub = 3;
end
fig_props.y_offset = 0.08;
fig_props.plot_third_deriv = options.plot_third_deriv;

%- identify relevant beat
if strcmp(no_of_pulses, 'single')
    beat_no = 1;
else
    beat_no = 4;
end
rel_els = fid_pts.f1(beat_no) : fid_pts.f2(beat_no);
curr_sig.t = sig.t(rel_els)-sig.t(rel_els(1));
curr_sig.sig = subtract_baseline(sig.v(rel_els));
[curr_sig.sig, scale_factor] = normalise(curr_sig.sig);
curr_sig.derivs.first = derivs.first(rel_els)/scale_factor;
curr_sig.derivs.second = derivs.second(rel_els)/scale_factor;
curr_sig.derivs.third = derivs.third(rel_els)/scale_factor;

% Plot signal

pts_to_plot = {'s', 'dic', 'dia'};
include_time_axis = 1;
plot_sig(curr_sig, fid_pts, cv_inds, fig_props, beat_no, pts_to_plot, include_time_axis)
% add dashed lines
% delta T
plot(curr_sig.t(fid_pts.s)*[1,1], [rt_annot, curr_sig.sig(fid_pts.s)], '--', 'color', 0.6*ones(1,3))
plot(curr_sig.t(fid_pts.dia)*[1,1], [rt_annot, curr_sig.sig(fid_pts.dia)], '--', 'color', 0.6*ones(1,3))
% RI
plot([curr_sig.t(fid_pts.s),curr_sig.t(fid_pts.f2)], 0*ones(1,2), '--', 'color', 0.6*ones(1,3))
plot([curr_sig.t(fid_pts.dia),curr_sig.t(fid_pts.f2)], curr_sig.sig(fid_pts.dia)*ones(1,2), '--', 'color', 0.6*ones(1,3))
% t_dia
plot(curr_sig.t(fid_pts.dic)*ones(1,2), [curr_sig.sig(fid_pts.dic), 0.98], '--', 'color', 0.6*ones(1,3))
plot(curr_sig.t(fid_pts.f2)*ones(1,2), [0, 0.98], '--', 'color', 0.6*ones(1,3))

% - DeltaT
normalised1  = coords_to_pos(curr_sig.t(fid_pts.s), 0);
normalised2  = coords_to_pos(curr_sig.t(fid_pts.dia), 0);
ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
text(mean([curr_sig.t(fid_pts.s),curr_sig.t(fid_pts.s),curr_sig.t(fid_pts.dia)]), 0.05, '\DeltaT','FontSize', fig_props.ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
% - RI
t_val = curr_sig.t(end) + 0.05*range(curr_sig.t);
normalised1  = coords_to_pos(curr_sig.t(end), 0);
normalised2  = coords_to_pos(curr_sig.t(end), curr_sig.sig(fid_pts.dia));
ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
text(t_val+0.02, curr_sig.sig(fid_pts.dia)/2, 'RI','FontSize', fig_props.ftsize, 'Color', 'k','HorizontalAlignment','center','VerticalAlignment', 'middle');
% - CT
normalised1  = coords_to_pos(0, 0);
normalised2  = coords_to_pos(curr_sig.t(fid_pts.s), 0);
ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
text(curr_sig.t(fid_pts.s)/1.6, 0.0+0.05, 'CT','FontSize', fig_props.ftsize, 'Color', 'k','HorizontalAlignment','center','VerticalAlignment', 'bottom');
% - t_dia
t_val = mean([curr_sig.t(fid_pts.dic), curr_sig.t(fid_pts.f2)]);
normalised1  = coords_to_pos(curr_sig.t(fid_pts.dic), 0.98);
normalised2  = coords_to_pos(curr_sig.t(fid_pts.f2), 0.98);
ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
text(t_val, 0.98+0.02, 't_{dia}','FontSize', fig_props.ftsize, 'Color', 'k','HorizontalAlignment','center','VerticalAlignment', 'bottom');

% Plot first derivative

pts_to_plot = {'ms'};
plot_first_deriv(curr_sig, fid_pts, cv_inds, fig_props, beat_no, pts_to_plot)

% add dashed lines
plot([-0.01, curr_sig.t(fid_pts.ms)], curr_sig.derivs.first(fid_pts.ms)*ones(1,2), '--', 'color', 0.6*ones(1,3))

% - ms
t_val = -0.015;
normalised1  = coords_to_pos(t_val, 0);
normalised2  = coords_to_pos(t_val, curr_sig.derivs.first(fid_pts.ms));
ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
text(1.2*t_val, 0.5*curr_sig.derivs.first(fid_pts.ms), 'ms','FontSize', fig_props.ftsize, 'Color', 'k','HorizontalAlignment','right','VerticalAlignment', 'bottom');

% Plot second derivative

pts_to_plot = {'b', 'd'};
plot_second_deriv(curr_sig, fid_pts, cv_inds, fig_props, beat_no, pts_to_plot, options)

% add line
plot([curr_sig.t(fid_pts.b), curr_sig.t(fid_pts.d)], [curr_sig.derivs.second(fid_pts.b), curr_sig.derivs.second(fid_pts.d)], '-k', 'LineWidth', 2)
text(curr_sig.t(fid_pts.b), 0.75*max(curr_sig.derivs.second), 'slope_{b-d}','FontSize', fig_props.ftsize, 'Color', 'k','HorizontalAlignment','left');

if options.plot_third_deriv
    pts_to_plot = {''};
    plot_third_deriv(curr_sig, fid_pts, cv_inds, fig_props, beat_no, pts_to_plot)
end
 
if ~isempty(options.save_folder)
    savepath = [options.save_folder, options.save_file, 'cv_inds'];
    set(gcf,'color','w');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize', [paper_size(1), paper_size(2)]./40);
    set(gcf,'PaperPosition',[0 0 paper_size(1) paper_size(2)]./40);
    print(gcf,'-dpdf',savepath)
end
shg

end

function make_demo_plot(sig, derivs, pulses, fid_pts, cv_inds, no_of_pulses, options)
% make demo plot of identifying fiducial points and taking feature
% measurements

%- setup
paper_size = [1000,350];
figure('Position', [50,50, paper_size])
ftsize = 13; lwidth = 2;
pp_annot = 0.8;
rt_annot = -0.02;
sig.t = [0:length(sig.v)-1]/sig.fs;

%% Plot Fiducial Points

subplot(1,2,1)

% identify relevant beat
if strcmp(no_of_pulses, 'single')
    beat_no = 1;
elseif length(fid_pts.f2) < 101
    beat_no = 2;
else
    beat_no = 100;
end
rel_els = fid_pts.f1(beat_no) : fid_pts.f2(beat_no);
curr_sig.t = sig.t(rel_els)-sig.t(rel_els(1));
curr_sig.sig = subtract_baseline(sig.v(rel_els));
[curr_sig.sig, scale_factor] = normalise(curr_sig.sig);

% plot baseline curve
if strcmp(options.demo_plot_wave, '2nd_deriv')
    rel_sig = derivs.second;
    init_vspace = 0.05*range(rel_sig);
else
    rel_sig = curr_sig.sig;
    init_vspace = 0.08;
end
plot([-10, 10], [0,0], '--k'), hold on
plot(curr_sig.t, rel_sig, 'LineWidth', lwidth); hold on,

% plot salient points
if strcmp(options.demo_plot_wave, '2nd_deriv')
    pt_names = {'a', 'b', 'c', 'd', 'e'};
else
    pt_names = {'dia', 'dic', 's', 'f1', 'f2'};
end
hspace0 = 0;
for pt_no = 1 : length(pt_names)
    eval(['curr_pt.el = fid_pts.' pt_names{pt_no} '(beat_no) - fid_pts.f1(beat_no)+1;']);
    if isnan(curr_pt.el)
        continue
    end
    
    % plot point
    curr_pt.v = rel_sig(curr_pt.el);
    curr_pt.t = curr_sig.t(curr_pt.el);
    plot(curr_pt.t, curr_pt.v, 'or')
    
    % annotate point
    curr_text = pt_names{pt_no};
    vspace0 = init_vspace;
    switch curr_text
        case {'a', 'c', 'e'}
            hspace0 = 0;
        case {'b', 'd'}
            hspace0 = 0;
            vspace0 = -1*init_vspace;            
        case {'dia','p2'}
            hspace0 = 0.04;
            if strcmp(curr_text, 'p2') && curr_pt.el < (fid_pts.s(beat_no) - fid_pts.f1(beat_no)+1)
                hspace0 = -0.04;
            elseif strcmp(curr_text, 'p2')
                hspace0 = 0;
            end
        case 'dic'
            hspace0 = -0.01;
            vspace0 = -1*init_vspace; 
        case {'f1', 'p1'}
            hspace0 = -0.04;
        case 'f2'
            hspace0 = 0.04;
        case 's'
            hspace0 = 0;
    end
    text(curr_sig.t(curr_pt.el)+hspace0, rel_sig(curr_pt.el)+vspace0 , curr_text,'FontSize', ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
    
end

% set limits
if strcmp(options.demo_plot_wave, '2nd_deriv')
    temp = [min(rel_sig), max(rel_sig)];
    ylim([temp(1) - 0.15*range(temp), temp(2)+0.15*range(temp)])
else
    ylim([-0.08, 1.2])
end
curr_range = range(curr_sig.t);
xlim([curr_sig.t(1)-0.1*curr_range, curr_sig.t(end)+0.1*curr_range])

% set labels
if strcmp(options.demo_plot_wave, '2nd_deriv')
    ylab = ylabel('PW''''', 'FontSize', ftsize, 'Rotation', 0);
else
    ylab = ylabel({'PW', '[au]'}, 'FontSize', ftsize, 'Rotation', 0);
end
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(gca, 'FontSize', ftsize -2, 'YTick', [])
xlabel('Time [s]', 'FontSize', ftsize)
box off

%% Plot CV Indices

subplot(1,2,2)

if strcmp(options.demo_plot_wave, '2nd_deriv')
    rel_sig = rel_sig/rel_sig(fid_pts.a);
end

% Plot signal
plot([-10, 10], [0,0], '--k'), hold on
plot(curr_sig.t, rel_sig, 'LineWidth', lwidth); hold on,

% plot salient points
hspace0 = 0;
for pt_no = 1 : length(pt_names)
    eval(['curr_pt.el = fid_pts.' pt_names{pt_no} '(beat_no) - fid_pts.f1(beat_no)+1;']);
    if isnan(curr_pt.el)
        continue
    end
    
    % plot point
    curr_pt.v = rel_sig(curr_pt.el);
    curr_pt.t = curr_sig.t(curr_pt.el);
    plot(curr_pt.t, curr_pt.v, 'or')
    
end

% set limits
if strcmp(options.demo_plot_wave, '2nd_deriv')
    temp = [min(rel_sig), max(rel_sig)];
    ylims = [temp(1) - 0.15*range(temp), temp(2)+0.15*range(temp)];
else
    ylims = [-0.08, 1.2];
end
ylim(ylims)
curr_range = range(curr_sig.t);
xlims = [curr_sig.t(1)-0.1*curr_range, curr_sig.t(end)+0.1*curr_range];
xlim(xlims)

% add dashed lines
% - a
%plot(curr_sig.t(fid_pts.a)*[1,1], [0, rel_sig(fid_pts.a)], '--', 'color', 0.6*ones(1,3))

y_offset = 0.07*range(ylim);
y_temp = 0.1*y_offset;
% - b
normalised1  = coords_to_pos(curr_sig.t(fid_pts.b), 0);
normalised2  = coords_to_pos(curr_sig.t(fid_pts.b), rel_sig(fid_pts.b));
ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]+y_temp);
ylim(ylims), xlim(xlims)
text(curr_sig.t(fid_pts.b), -y_offset+rel_sig(fid_pts.b), 'b/a','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center');
ylim(ylims), xlim(xlims)
% - d
normalised1  = coords_to_pos(curr_sig.t(fid_pts.d), 0);
normalised2  = coords_to_pos(curr_sig.t(fid_pts.d), rel_sig(fid_pts.d));
ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]+y_temp);
ylim(ylims), xlim(xlims)
text(curr_sig.t(fid_pts.d), -y_offset+rel_sig(fid_pts.d), 'd/a','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center');
ylim(ylims), xlim(xlims)
% - e
normalised1  = coords_to_pos(curr_sig.t(fid_pts.e), 0);
normalised2  = coords_to_pos(curr_sig.t(fid_pts.e), rel_sig(fid_pts.e));
ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]+y_temp);
ylim(ylims), xlim(xlims)
text(curr_sig.t(fid_pts.e), y_offset+rel_sig(fid_pts.e), 'e/a','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center');
ylim(ylims), xlim(xlims)

% set labels
ylab = ylabel('PW''''/a', 'FontSize', ftsize, 'Rotation', 0);
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(gca, 'FontSize', ftsize -2, 'YTick', [])
xlabel('Time [s]', 'FontSize', ftsize)
box off

%% Save Figure

if ~isempty(options.save_folder)
    savepath = [options.save_folder, options.save_file, 'demo_plot'];
    set(gcf,'color','w');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize', [paper_size(1), paper_size(2)]./40);
    set(gcf,'PaperPosition',[0 0 paper_size(1) paper_size(2)]./40);
    print(gcf,'-dpdf',savepath)
    print(gcf,'-depsc',savepath)
end

end

function normalised  = coords_to_pos(x_coord, y_coord)

pos = get(gca, 'Position');
normalised(1) = (x_coord - min(xlim))/diff(xlim) * pos(3) + pos(1);
normalised(2) = (y_coord - min(ylim))/diff(ylim) * pos(4) + pos(2);
 
end

function plot_sig(curr_sig, fid_pts, cv_inds, fig_props, beat_no, pts_to_plot, include_time_axis)
%- plot sig

if nargin < 6
    pts_to_plot = {'dia', 'dic', 's', 'p1', 'p2'};
end

h_b(1) = subplot('Position', [0.21,fig_props.y_offset+(fig_props.n_sub-1)*fig_props.y_inc,0.78,fig_props.y_inc-0.01]);

plot(curr_sig.t, curr_sig.sig, 'LineWidth', fig_props.lwidth); hold on,

% plot salient points
pt_names = intersect({'dia', 'dic', 's', 'p1', 'p2'}, pts_to_plot);
for pt_no = 1 : length(pt_names)
    eval(['curr_pt.el = fid_pts.' pt_names{pt_no} '(beat_no) - fid_pts.f1(beat_no)+1;']);
    if isnan(curr_pt.el)
        continue
    end
    
    % plot point
    curr_pt.v = curr_sig.sig(curr_pt.el);
    curr_pt.t = curr_sig.t(curr_pt.el);
    plot(curr_pt.t, curr_pt.v, 'or')
    
    % annotate point
    curr_text = pt_names{pt_no};
    vspace0 = 0.08;
    hspace0 = 0;
    switch curr_text
        case {'dia','p2'}
            hspace0 = 0.04;
            if strcmp(curr_text, 'p2') && curr_pt.el < (fid_pts.s(beat_no) - fid_pts.f1(beat_no)+1)
                hspace0 = -0.04;
            end
            vspace0 = 0.11;
        case 'dic'
            hspace0 = 0;
            vspace0 = -0.11;
        case {'f1', 'p1'}
            hspace0 = -0.04;
        case 's'
            hspace0 = 0;
            vspace0 = 0.11;
    end
    text(curr_sig.t(curr_pt.el)+hspace0, curr_sig.sig(curr_pt.el)+vspace0 , curr_text,'FontSize', fig_props.ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
    
end

% set limits
ylim([-0.08, 1.2])
curr_range = range(curr_sig.t);
xlim([curr_sig.t(1)-0.1*curr_range, curr_sig.t(end)+0.1*curr_range])

% set labels
ylab = ylabel('PW', 'FontSize', fig_props.ftsize, 'Rotation', 0);
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(gca, 'FontSize', fig_props.ftsize -2, 'YTick', [])
if include_time_axis
    xlabel('Time [s]', 'FontSize', fig_props.ftsize)
    set(gca, 'XTick', 0:0.25:1)
else
    set(gca, 'XTick', 0:0.25:1, 'XTickLabel', {})
end
box off

end

function plot_first_deriv(curr_sig, fid_pts, cv_inds, fig_props, beat_no, pts_to_plot)
%- plot first derivative

if nargin < 6
    pts_to_plot = {'ms' , 'dia'};
end

h_b(2) = subplot('Position', [0.21,fig_props.y_offset+(fig_props.n_sub-2)*fig_props.y_inc,0.78,fig_props.y_inc-0.01]);


% Plot x-axis
plot([-10, 10], [0,0], '--k'); hold on

% plot first derivative
plot(curr_sig.t, curr_sig.derivs.first, 'LineWidth', fig_props.lwidth); hold on,

% plot salient points
pt_names = intersect({'ms', 'dia'}, pts_to_plot);
curr_range = range(curr_sig.derivs.first);
vspace0 = 0.08*curr_range;
hspace0 = 0.05;
for pt_no = 1 : length(pt_names)
    eval(['curr_pt.el = fid_pts.' pt_names{pt_no} '(beat_no) - fid_pts.f1(beat_no)+1;']);
    if isnan(curr_pt.el)
        continue
    end
    
    % plot point
    curr_pt.v = curr_sig.derivs.first(curr_pt.el);
    curr_pt.t = curr_sig.t(curr_pt.el);
    plot(curr_pt.t, curr_pt.v, 'or')
    
    % annotate point
    curr_text = pt_names{pt_no};
    text(curr_sig.t(curr_pt.el)+hspace0, curr_sig.derivs.first(curr_pt.el)+vspace0 , curr_text,'FontSize', fig_props.ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
    
end

% set limits
curr_range = range(curr_sig.t);
xlim([curr_sig.t(1)-0.1*curr_range, curr_sig.t(end)+0.1*curr_range])
ylims = ylim; curr_range = range(curr_sig.derivs.first);
ylim([min(curr_sig.derivs.first)-0.05*curr_range, max(curr_sig.derivs.first)+0.15*curr_range])

% set labels
ylab = ylabel({'1st','derivative'}, 'FontSize', fig_props.ftsize, 'Rotation', 0);
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(gca, 'FontSize', fig_props.ftsize -2, 'YTick', 0, 'YTickLabel', {'0'})
set(gca, 'XTick', 0:0.25:1, 'XTickLabel', {})

box off

end

function plot_second_deriv(curr_sig, fid_pts, cv_inds, fig_props, beat_no, pts_to_plot, options)

%- plot Second derivative

if nargin < 6
    pts_to_plot = {'a', 'b', 'c', 'd', 'e', 'f'};
end

h_b(3) = subplot('Position', [0.21,fig_props.y_offset+(fig_props.n_sub-3)*fig_props.y_inc,0.78,fig_props.y_inc-0.01]);


% Plot x-axis
plot([-10, 10], [0,0], '--k'); hold on

% plot second derivative
plot(curr_sig.t, curr_sig.derivs.second, 'LineWidth', fig_props.lwidth); hold on,

% plot salient points
pt_names = intersect({'a', 'b', 'c', 'd', 'e', 'f'}, pts_to_plot);
curr_range = range(curr_sig.derivs.second);
vspace_const = 0.08;
hspace0 = 0;
for pt_no = 1 : length(pt_names)
    eval(['curr_pt.el = fid_pts.' pt_names{pt_no} '(beat_no) - fid_pts.f1(beat_no)+1;']);
    if isnan(curr_pt.el)
        continue
    end
    
    % plot point
    curr_pt.v = curr_sig.derivs.second(curr_pt.el);
    curr_pt.t = curr_sig.t(curr_pt.el);
    plot(curr_pt.t, curr_pt.v, 'or')
    
    % annotate point
    curr_text = pt_names{pt_no};
    switch curr_text
        case {'a','c', 'e'}
            vspace0 = (1.35*vspace_const)*curr_range;
        case {'b', 'd', 'f'}
            vspace0 = -1.2*vspace_const*curr_range;
    end
    text(curr_sig.t(curr_pt.el), curr_sig.derivs.second(curr_pt.el)+vspace0 , curr_text,'FontSize', fig_props.ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
    
end

% set limits
curr_range = range(curr_sig.t);
xlim([curr_sig.t(1)-0.1*curr_range, curr_sig.t(end)+0.1*curr_range])
ylims = ylim; curr_range = range(curr_sig.derivs.second);
ylim([min(curr_sig.derivs.second)-0.05*curr_range, max(curr_sig.derivs.second)+0.15*curr_range])

% set labels
ylab = ylabel({'2nd', 'derivative'}, 'FontSize', fig_props.ftsize, 'Rotation', 0);
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(gca, 'FontSize', fig_props.ftsize -2, 'XTick', 0:0.25:1, 'YTick', 0, 'YTickLabel', {'0'})
if options.plot_third_deriv
    set(gca, 'XTickLabel', {})
else
    xlabel('Time [s]', 'FontSize', fig_props.ftsize)
end
box off

end

function plot_third_deriv(curr_sig, fid_pts, cv_inds, fig_props, beat_no, pts_to_plot)

%- plot Third derivative

if nargin < 6
    pts_to_plot = {'p1', 'p2'};
end

h_b(4) = subplot('Position', [0.21,fig_props.y_offset+(fig_props.n_sub-4)*fig_props.y_inc,0.78,fig_props.y_inc-0.01]);

% Plot x-axis
plot([-10, 10], [0,0], '--k'); hold on

% plot third derivative
curr_sig.derivs.third(1:8) = nan;
curr_sig.derivs.third(end-7:end) = nan;
plot(curr_sig.t, curr_sig.derivs.third, 'LineWidth', fig_props.lwidth); hold on,

% plot salient points
pt_names = intersect({'p1', 'p2'}, pts_to_plot);
curr_range = range(curr_sig.derivs.third);
vspace_const = 0.11;
hspace0 = 0;
for pt_no = 1 : length(pt_names)
    eval(['curr_pt.el = fid_pts.' pt_names{pt_no} '(beat_no) - fid_pts.f1(beat_no)+1;']);
    if isnan(curr_pt.el)
        continue
    end
    
    % plot point
    curr_pt.v = curr_sig.derivs.third(curr_pt.el);
    curr_pt.t = curr_sig.t(curr_pt.el);
    plot(curr_pt.t, curr_pt.v, 'or')
    
    % annotate point
    curr_text = pt_names{pt_no};
    switch curr_text
        case {'p1'}
            vspace0 = vspace_const*curr_range;
        case {'p2'}
            vspace0 = -1*vspace_const*curr_range;
    end
    text(curr_sig.t(curr_pt.el), curr_sig.derivs.third(curr_pt.el)+vspace0 , curr_text,'FontSize', fig_props.ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
    
end

% set limits
curr_range = range(curr_sig.t);
xlim([curr_sig.t(1)-0.1*curr_range, curr_sig.t(end)+0.1*curr_range])
ylims = ylim; curr_range = range(curr_sig.derivs.third);
ylim([min(curr_sig.derivs.third)-0.05*curr_range, max(curr_sig.derivs.third)+0.15*curr_range])

% set labels
xlabel('Time [s]', 'FontSize', fig_props.ftsize)
ylab = ylabel({'3rd', 'derivative'}, 'FontSize', fig_props.ftsize, 'Rotation', 0);
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(gca, 'FontSize', fig_props.ftsize -2, 'XTick', 0:0.25:1, 'YTick', [])
box off

end

function subtracted = subtract_baseline(sig)

baseline = linspace(sig(1),sig(end),length(sig));

subtracted = sig(:) - baseline(:); 

end

function [norm, scale_factor] = normalise(sig)

norm = sig - min(sig); 
scale_factor = max(norm);
norm = norm / scale_factor;

end