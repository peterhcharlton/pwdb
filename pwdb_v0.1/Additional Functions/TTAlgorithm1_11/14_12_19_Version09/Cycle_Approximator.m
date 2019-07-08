function indmin = Cycle_Approximator(t,signal1)
% The primary use of this is as a pre-routine to find the cardiac time
% period T.  With this, we can adjust the MAX/MIN threshold to be per
% cardiac cycle, so that we do not miss any pulses.

signal = signal1(1,:);

%%%%% The Identification of Peaks %%%%%
% indmax are the maxima positions, indmin are the minima positions
indmax = [];
indmin = [];

threshold = 2.5.*sqrt(var( signal ));

% Threshold should be set to the height of the smallest peak that should be
% detected
sa = signal(end);  % sa is the maximal element since the last minimum
sb = signal(end);  % sb is the minimum element since the last maximum
sc = 0;

i=length(signal);
a=length(signal);
d=0;   % d is the direction of the signal 0-indetermined direction,
%        1-from minimum to maximum, 2-from maximum to minimum
sigLen=length(signal);

while (i > 1)
    i=i-1;
    if (d==0) % direction not yet determined
        if (sa >= (signal(i)+threshold))
            d=2; % FALL
            if (i < sigLen-3)
            else
            end
        elseif (signal(i) >= (sb+threshold))
            d=1; % RISE
        end;
        if (sa <= signal(i))
            sa = signal(i);
            a = i;
        elseif (signal(i) <= sb)
            sb = signal(i);
            b = i;
        end;
    elseif (d==1) % direction from min to max (RISE)
        if (sa <= signal(i))
            sa = signal(i);
            a = i;
        elseif (sa >= (signal(i)+threshold))% && ~isempty(indmin)
            indmax = [indmax a];
            sb = signal(i);
            b = i;
            d=2;
        end;
    elseif (d==2) % direction from max to min (FALL)
        if (signal(i) <= sb)
            sb = signal(i);
            b = i;
        elseif (signal(i) >= (sb + threshold)) || (t(a) - t(i)) > 0.4
            indmin = [indmin b];
            sa = signal(i);
            a = i;
            d=1;
        end;
    end
end
