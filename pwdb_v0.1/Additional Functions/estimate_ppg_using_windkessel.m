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
ppg.v = normwave(Volvb.v);

end

function norm_wave = normwave(orig_wave)

norm_wave = orig_wave;

norm_wave = (norm_wave-min(norm_wave))./range(norm_wave);

end