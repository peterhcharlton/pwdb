%::::::::::::::::::::::::::::::::%
% HARMONICS OF AN INPUT FUNCTION %
%::::::::::::::::::::::::::::::::%
%
% Jordi Alastruey
% King's College London
% January 2014

function [FT,H_amp,H_ph] = Harmonics(func,NHarm,fmt,OutputFile,T)

%NHarm: Number of harmonics plotted
%func: Input function
%SF: Sampling frequency

format compact

n_pts = length(func);
FT=0;
func_cal=0;

%% Fourier transform
clear k j i
for k=1:n_pts,
    omega(k)=2*pi*(k-1)/n_pts; % omega
    FT(k)=0;
    for j=1:n_pts,
	    FT(k)=FT(k)+func(j)*exp(-i*omega(k)*(j-1))/n_pts;
    end;
end;
%% Inverse Fourier transform
clear j k i
for j=1:n_pts,
    func_cal(j)=0;
    for k=1:n_pts,
	    func_cal(j)=func_cal(j)+FT(k)*exp(i*omega(j)*(k-1));
    end;
end;
clear j k i
%% Magnitude and phase of FT
mFT= abs(FT);  % modulus of FT
phFT = angle(FT);  % phase angle of FT
%% Plot of the intensity
if (fmt==1)
figure
plot (omega(1:round(n_pts/2))/2/pi*n_pts, mFT(1:round(n_pts/2)))
grid
xlabel('Harmonic');
ylabel('Intensity');
end
%% Check if the initial function is obtained from the harmonics calculated
if (fmt==1)
figure
plot (omega/2/pi*n_pts, func_cal,'b', omega/2/pi*n_pts, func, 'r');
xlabel('Sample');
ylabel('Value');
legend('calculated','initial function')
grid
end
%% Calculation of the harmonics
H = zeros(round(n_pts/2-1),n_pts);
H(1,:) = FT(1); % DC frequency
for j=2:n_pts/2-1,
    for k=1:n_pts,
	    H(j,k)=2*FT(j)*exp(i*omega(k)*(j-1));
    end;
end;
%% Plot the harmonic indicated in H(?,:) (H(1,:) corresponds to DC, H(2,:) to 1st harmonic, and so on)
if (fmt==1)
figure
plot (omega/2/pi*n_pts, H(2,:), omega/2/pi*n_pts, func);
legend('harmonic','initial function')
grid
xlabel('Sample');
ylabel('Value');
end
%% Check if the harmonics are well-calculated
HT = H(1,:);
for j=2:round(n_pts/2)-1,
    HT(1,:) = HT(1,:) + H(j,:); % The summation of all the harmonics should yield the original sample vector
end
if (fmt==1)
figure
plot (omega/2/pi*n_pts, HT,'b', omega/2/pi*n_pts, func, 'r');
xlabel('Sample');
ylabel('Value');
legend('calculated from harmonics','initial function') 
grid
end
%% Check the function generated by the first NHarm harmonics
HT = H(1,:);
for j=2:NHarm,
    HT(1,:) = HT(1,:) + H(j,:); % The summation of all the harmonics should yield the original sample vector
end
if (fmt==1)
figure
plot (omega/2/pi*n_pts, HT,'b', omega/2/pi*n_pts, func, 'r');
xlabel('Sample');
ylabel('Value');
legend('calculated from NHarm harmonics','initial function') 
grid
end
%% Amplitude and phase angle of each harmonic
% func = H_amp(j)*sin(2*pi*j/T + H_ph(j))
H_ph = zeros(1,round(n_pts/2-1));
H_amp = zeros(1,round(n_pts/2-1));
for j=1:round(n_pts/2)-1,
    H_ph(j) = atan(-real(FT(j+1))/imag(FT(j+1))); % Phase angle
    H_amp(j) = 2*real(FT(j+1))/sin(H_ph(j)); % Amplitude
end;
%% Display the function for 1D Nektar for the first NHarm harmonics
%clc
if (fmt==1)
fprintf('bc=%0.5g\n',FT(1))
for j=1:NHarm
if(H_amp(j)>0)
    if(H_ph(j)>0)
        fprintf('+%0.5g*sin(%d*PI*t/T+%0.5g)',H_amp(j),2*j,H_ph(j))
    else
        fprintf('+%0.5g*sin(%d*PI*t/T%0.5g)',H_amp(j),2*j,H_ph(j))
    end
else
    if(H_ph(j)>0)
        fprintf('%0.5g*sin(%d*PI*t/T+%0.5g)',H_amp(j),2*j,H_ph(j))
    else
        fprintf('%0.5g*sin(%d*PI*t/T%0.5g)',H_amp(j),2*j,H_ph(j))
    end
end
end
fprintf('\n\n')
FT = FT(1);
end

%% Display the amplitude and phase for the first NHarm harmonics
% in the right format for an input .bcs file for Nektar 1-D
fid = 1; % Print results on the screen by default
if (nargin > 3)
    fid = fopen(OutputFile,'wt');
end

if (fmt==0) && (fid==1)
    fprintf(fid,'Harmonics calculated\n');
else
    if (nargin == 5)
        fprintf(fid,'%d %0.5g %0.5g\n',NHarm, T, FT(1));
    else
        fprintf(fid,'%d [To be replaced with T] %0.5g\n',NHarm, FT(1));
    end
 
    for j=1:NHarm
        fprintf(fid,'%0.5g %0.5g\n',H_amp(j),H_ph(j));
    end
end

if (nargin > 3)
    fclose(fid);
end


