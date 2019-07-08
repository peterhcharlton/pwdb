%:::::::::::::::::::::::::::::::::%
% FILTER WAVEFORM USING HARMONICS %
%:::::::::::::::::::::::::::::::::%
%
% Jordi Alastruey
% King's College London
% January 2014

function [tfunc,funcf] = Harmonics_filter(tfunc,func,SR,fmt,OutputFile)

T = tfunc(end) + tfunc(2)-tfunc(1); % Cardiac period in s
NHarm = round((length(func)-1.1)/2);

if (nargin == 5)
    [FT,H_amp,H_ph] = Harmonics(func(1:end),NHarm,fmt,OutputFile,T);
else
    [FT,H_amp,H_ph] = Harmonics(func(1:end),NHarm,fmt);
end

tfunc = 0:1/SR:T;
funcf = FT(1);
for j=1:NHarm
    funcf = funcf + H_amp(j)*sin(2*j*pi*tfunc/T + H_ph(j));
end



%-- Freq file:
% NHarm T FT(1)
% [H_amp(:) H_ph(:)]
