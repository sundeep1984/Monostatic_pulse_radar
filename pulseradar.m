% Basic design and simulation of a monostatic pulse radar 

% Copyright (C) 2017  Sundeep Venkatraman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

clear all;
close all;
%clc;
pd = 0.9;            % Probability of detection
pfa = 1e-6;          % Probability of false alarm
max_range = 5000;    % Maximum unambiguous range
range_res = 50;      % Required range resolution
tgt_rcs = 1;         % Required target radar cross section


%prop_speed = 2.99792e8;
prop_speed=3e8;                         % Propagation speed
pulse_bw = prop_speed/(2*range_res);    % Pulse bandwidth
pulse_width = 1/pulse_bw;               % Pulse width
prf = prop_speed/(2*max_range);         % Pulse repetition frequency
fs = 2*pulse_bw;                        % Sampling rate

noise_bw = pulse_bw;

num_pulse_int = 10;

% minimum required SNR in dB calculated using Albersheim's equation
A = log(0.62/pfa);
B=log(pd/(1-pd));
snr_min = -5*log10(num_pulse_int) + (6.2 +4.54/sqrt(num_pulse_int+0.44))*log10(A+0.12*A*B+1.7*B);

tx_gain = 20;%gains are in dB
txg=10^(tx_gain/10);

rx_gain = 20;
rxg=10^(rx_gain/10);

fc = 10e9;
lambda = prop_speed/fc;

k=1.3806e-23;
T=300;
N0=k*T*noise_bw;


P_r=N0/2*10^(snr_min/10);

peak_power = P_r*(4*pi)^2*max_range^4/(txg*rxg*lambda^2*tgt_rcs);



tgtpos = [[2024.66;0;0],[3518.63;0;0],[3845.04;0;0],[4900;0;0]];

tgtvel = [[0;0;0],[0;0;0],[0;0;0]];

sampling_interval = 10e-9;% this is assumed for the simulation - signal is not physically sampled at this rate 

tot_time = 2*max_range/prop_speed;

tot_samples = floor(tot_time/sampling_interval);
pulse_samples = floor(pulse_width/sampling_interval);


pulsewave = zeros(1,tot_samples);
pulsewave(1:pulse_samples) = sqrt(peak_power);


ret_pulse = randn(num_pulse_int,length(pulsewave))*sqrt(N0/2);


s=size(tgtpos);
distance=sqrt(tgtpos(1,:).^2 + tgtpos(2,:).^2 + tgtpos(3,:).^2);
delays=floor((2*distance/prop_speed)/sampling_interval);
    
    
for i=1:num_pulse_int    
for i1=1:s(2)
    % use the reflected signal equation at https://www.mathworks.com/help/phased/ref/phased.radartarget-class.html
    % in addition to 1/(4*pi*R^2 )
    if delays(i1)<tot_samples
        ret_pulse(i,:) = ret_pulse(i,:) + [zeros(1,delays(i1)) pulsewave(1:(tot_samples-delays(i1)))]*sqrt(txg*rxg)*(lambda)/((4*pi)*(distance(i1))^2);
        
    end
    

end


    
end


% find threshold for coherent integration/detection and use approx
% threshold(calculated offline in terms of inverse incomplete gamma
% function) for noncoherent integration
c_thresh = sqrt(N0/2*num_pulse_int)*erfcinv(2*pfa);
nc_thresh = sqrt(N0/2)*5.7193; 



% need to do matched filtering before noncoherent or coherent integration
mf = fliplr(conj(pulsewave(1:pulse_samples))/sqrt(peak_power));
for i=1:num_pulse_int

    temp=conv(mf,ret_pulse(i,:));    
    mfout(i,:)=temp(floor(pulse_samples/2)+1:end-floor(pulse_samples/2));
end




noncoh = sqrt(sum(abs(mfout).^2))/num_pulse_int;
coh = sum(mfout)/num_pulse_int;

figure
plot((1:tot_samples)*sampling_interval,coh);
hold on
plot((1:tot_samples)*sampling_interval,c_thresh*ones(1,length(coh)),'r-');


% the code with det1 and det2 is used to detect when the output exceeds
% threshold and marks the beginning of each time it exceeds threshold as a
% detected object
det1=find(coh>c_thresh);
det1 = [det1(1) det1];
det2 = diff(det1);
stops = find(det2>1);

if length(det1)>0
    locs = det1([1 stops+1]);
    coh_detected_dist=locs*sampling_interval*prop_speed/2
    
end

figure
plot((1:tot_samples)*sampling_interval,noncoh);
hold on
plot((1:tot_samples)*sampling_interval,nc_thresh*ones(1,length(coh)),'r-');

det1=find(noncoh>nc_thresh);
det1 = [det1(1) det1];
det2 = diff(det1);
stops = find(det2>1);

if length(det1)>0
    locs = det1([1 stops+1]);
    noncoh_detected_dist=locs*sampling_interval*prop_speed/2
    
end


