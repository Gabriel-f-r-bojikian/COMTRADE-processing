clear all;
close all;
clc;

% Choose the name of the file to be analysed
filename = 'AD_GR_Rf25_F2_BCT_0';

% Data input
COMTRADE_data = csvread([filename '.dat']);
config = csvread([filename '.cfg']);

n = COMTRADE_data(:,1);
ta = COMTRADE_data(:,2)*1E-6;
fa_comtrade = 1/(ta(2) - ta(1));

% Separating the channels, recording them in engineering values
VANs = COMTRADE_data(:, 3)*config(3, 6) + config(3, 7);
VBNs = COMTRADE_data(:, 4)*config(4, 6) + config(4, 7);
VCNs = COMTRADE_data(:, 5)*config(5, 6) + config(5, 7);
IALp = (COMTRADE_data(:, 6)*config(6, 6) + config(6, 7))*1E3;
IBLp = (COMTRADE_data(:, 7)*config(7, 6) + config(7, 7))*1E3;
ICLp = (COMTRADE_data(:, 8)*config(8, 6) + config(8, 7))*1E3;

% Converting primary currents to secondary currents
IALs = IALp / 1000;
IBLs = IBLp / 1000;
ICLs = ICLp / 1000;

% Plotting the three phase voltages
figure;
plot(ta, VANs, ta, VBNs, ta, VCNs);
grid;
legend('VANs', 'VBNs', 'VCNs');
title('Secondary voltages');
xlabel('Time [s]');
ylabel('Voltages [V]');

% Plotting the three secondary line currents
figure;
plot(ta, IALs, ta, IBLs, ta, ICLs);
grid;
legend('IALs', 'IBLs', 'ICLs');
title('Secondary currents');
xlabel('Time [s]');
ylabel('Currents [A]');

% Conditioning the signals
FMV = 3.0/(115*sqrt(2)/sqrt(3)); % V/V
FMI = 3.0/15; % A/V

% Generate an analog signal that will be processed by the lowpass filter
VANs_cond = VANs * FMV;
VBNs_cond = VBNs * FMV;
VCNs_cond = VCNs * FMV;
% Plotting the A phase current after conditioning
figure;
plot(ta, VBNs, ta, VBNs_cond); 
grid;
legend('VBNs', 'VBNs_cond');
title('Voltages inside the IED before and after the FMV');
xlabel('Time [s]');
ylabel('Voltage [V]');

% Generate an analog signal that will be processed by the lowpass filter
VIBLs_cond = IBLs * FMI;
% Plotting the A phase current after conditioning
figure;
plot(ta, IBLs, ta, VIBLs_cond); 
grid;
legend('VIBLs_cond', 'VIBLs_fil');
title('Current voltages inside the IED before and after FMI');
xlabel('Time [s]');
ylabel('Voltage [V]');


% Simulating the filter
pkg load signal;

[nf, wc] = buttord(2*pi*90, 2*pi*60*8, 0.1, 36, 's');
[num, den] = butter(nf, wc, 's');

% Esta funcao gera o filtro de fato
lpf = tf(num, den);
VIALs_fil = lsim(lpf, VIALs_cond, ta);