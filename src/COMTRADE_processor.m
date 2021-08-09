clear all;
close all;
clc;

% Choose the name of the file to be analyzed
filename = 'AD_GR_Rf25_F2_BCT_0';

% Key variables
nominalCurrent = 1; % Ampère
samplesPerCycle = 32;

% Data input
COMTRADE_data = csvread([filename '.dat']);
config = csvread([filename '.cfg']);

n = COMTRADE_data(:,1);
ta = COMTRADE_data(:,2)*1E-6;
fa_comtrade = 1/(ta(2) - ta(1));

% Functions we will use for this script

function firstHarmonicAvg = fourierFilter(inputSignal, windowSize)
  
  % Applying discreet fourier transform to see energy levels in the first harmonic

  FIR_Cosseno = (sqrt(2)/windowSize)*cos(2*pi*(1:windowSize)/windowSize);
  FIR_Seno    = (sqrt(2)/windowSize)*sin(2*pi*(1:windowSize)/windowSize);

  pos = 1;
  buffer = zeros(1,windowSize);
  while(pos<=size(inputSignal,2))
    bufferPosition = windowSize;

    while(bufferPosition > 1)
      buffer(bufferPosition) = buffer(bufferPosition - 1);
      bufferPosition = bufferPosition - 1;
    endwhile
    buffer(1) = inputSignal(pos);
    
    YR = 0;
    auxCounter = 1;
    while(auxCounter <= windowSize)
      YR = YR + FIR_Cosseno(auxCounter)*buffer(auxCounter);
      auxCounter = auxCounter + 1;
    endwhile

    YI = 0;
    auxCounter = 1;
    while(auxCounter <= windowSize)
      YI = YI + FIR_Seno(auxCounter)*buffer(auxCounter);
      auxCounter = auxCounter + 1;
    endwhile



    firstHarmonicAvg(pos) = abs(YR+i*YI);
    
    pos = pos+1;
        
  endwhile
endfunction

function electricCurrentSpikePosition = findElectricCurrentSpike(inputSignal, nominalCurrent)
  electricCurrentSpikePosition = Inf; % If not found, return impossible value
  for pos = 1:size(inputSignal,2)
    if (inputSignal(pos) > nominalCurrent)
      electricCurrentSpikePosition = pos;
      break;
    endif
  endfor
endfunction

function shutdownPosition = findShutdownPosition(inputSignal, timeoutLength)
  epsilon = 1E-3;
  shutdownPosition = Inf; % If not found, return impossible value

  timeoutAccumulator = 0;
  for pos = 1:size(inputSignal,2)
    if (inputSignal(pos) < epsilon)
      timeoutAccumulator = timeoutAccumulator + 1;
    else
      timeoutAccumulator = 0;
    endif

    if (timeoutAccumulator == timeoutLength)
      shutdownPosition = pos - timeoutLength;
      break;
    endif
  endfor
endfunction

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

VIALs_cond = IALs * FMI;
VIBLs_cond = IBLs * FMI;
VICLs_cond = ICLs * FMI;


% Plotting the A phase current after conditioning
figure;
plot(ta, VBNs, ta, VBNs_cond); 
grid;
legend('VBNs', 'VBNs_cond');
title('Voltages inside the IED before and after the FMV');
xlabel('Time [s]');
ylabel('Voltage [V]');

% Generate an analog signal that will be processed by the lowpass filter
% Plotting the A phase current after conditioning
figure;
plot(ta, IBLs, ta, VIBLs_cond); 
grid;
legend('VIBLs_cond', 'VIBLs_fil');
title('Current voltages inside the IED before and after FMI');
xlabel('Time [s]');
ylabel('Voltage [V]');


% Step 3 - Clamping the voltages internal to the IED
max_ADC_voltage = 3.3;

VANs_lim= VANs_cond;
VANs_lim(VANs_lim>max_ADC_voltage)= max_ADC_voltage;
VANs_lim(VANs_lim<-max_ADC_voltage)= -max_ADC_voltage;
VBNs_lim= VBNs_cond;
VBNs_lim(VBNs_lim>max_ADC_voltage)= max_ADC_voltage;
VBNs_lim(VBNs_lim<-max_ADC_voltage)= -max_ADC_voltage;
VCNs_lim= VCNs_cond;
VCNs_lim(VCNs_lim>max_ADC_voltage)= max_ADC_voltage;
VCNs_lim(VCNs_lim<-max_ADC_voltage)= -max_ADC_voltage;

VIALs_lim= VIALs_cond;
VIALs_lim(VIALs_lim>max_ADC_voltage)= max_ADC_voltage;
VIALs_lim(VIALs_lim<-max_ADC_voltage)= -max_ADC_voltage;
VIBLs_lim= VIBLs_cond;
VIBLs_lim(VIBLs_lim>max_ADC_voltage)= max_ADC_voltage;
VIBLs_lim(VIBLs_lim<-max_ADC_voltage)= -max_ADC_voltage;
VICLs_lim= VICLs_cond;
VICLs_lim(VICLs_lim>max_ADC_voltage)= max_ADC_voltage;
VICLs_lim(VICLs_lim<-max_ADC_voltage)= -max_ADC_voltage;


% Plotting the A phase current after Clamping
figure;
plot(ta, VBNs_cond, ta, VBNs_lim); 
grid;
legend('VBNs_cond', 'VBNs_lim');
title('Clamped voltages inside the IED');
xlabel('Time [s]');
ylabel('Voltage [V]');

% Plotting the A phase current after clamping
figure;
plot(ta, VIBLs_cond, ta, VIBLs_lim); 
grid;
legend('VIBLs_cond', 'VIBLs_lim');
title('Clamped Current voltages inside the IED');
xlabel('Time [s]');
ylabel('Voltage [V]');


% Simulating the filter
pkg load signal;

resolucao_bits = floor(11+1/3.2);

% Frequencies
fp = 90; % Hz
fs = 12 * fp; % Hz

wp = 2*pi*fp; % rad/s
ws = 2*pi*fs; % rad/s

% Attenuations
Amax = 0.1; %dB
Amin = 0.6*20*log10(2^resolucao_bits);

% Designing the Butterworth lowpass filter
[nf, wc] = buttord(wp, ws, Amax, Amin, 's');
[num, den] = butter(nf, wc, 's');

% Generating the transfer function
lpf = tf(num, den);

% Applying the filtering
VANs_fil = lsim(lpf, VANs_lim, ta);
VBNs_fil = lsim(lpf, VBNs_lim, ta);
VCNs_fil = lsim(lpf, VCNs_lim, ta);
VIALs_fil = lsim(lpf, VIALs_lim, ta);
VIBLs_fil = lsim(lpf, VIBLs_lim, ta);
VICLs_fil = lsim(lpf, VICLs_lim, ta);

% Plotting the B phase voltage after filtering
figure;
plot(ta, VIBLs_lim, ta, VIBLs_fil); 
grid;
legend('VIBLs_lim', 'VIBLs_fil');
title('Current voltages inside the IED before and after LPF');
xlabel('Time [s]');
ylabel('Voltage [V]');

% Digitalization through the ADC
% Decimation
fa_final = samplesPerCycle*60;

decim_fact = round(fa_comtrade/fa_final);

samp=1;
for ccomt = 1:decim_fact:size(VIALs_fil,1),
  VANs_samp(samp) = VANs_fil(ccomt,1);
  VBNs_samp(samp) = VBNs_fil(ccomt,1);
  VCNs_samp(samp) = VCNs_fil(ccomt,1);
  VIALs_samp(samp) = VIALs_fil(ccomt,1);
  VIBLs_samp(samp) = VIBLs_fil(ccomt,1);
  VICLs_samp(samp) = VICLs_fil(ccomt,1);
  
  ta_samp(samp)    = ta(ccomt);
  samp = samp + 1;
endfor

figure;
plot(ta,VIBLs_fil, ta_samp,VIBLs_samp, 'o');
grid;
legend('VIB filtrado', 'VIB amostrado');
title('Amostragem das informações');
xlabel('Tempo [s]');
ylabel ('Tensao [V]');


% Quantization
q = (max_ADC_voltage - (-max_ADC_voltage))/(2^resolucao_bits);

VANs_dig = VANs_samp / q;
VBNs_dig = VBNs_samp / q;
VCNs_dig = VCNs_samp / q;
IALs_dig = VIALs_samp / q; %Cada valor individual é quantizado no número de bits, simbolos, niveis
% resultantes para sua representação na saída do AD
IBLs_dig = VIBLs_samp / q;
ICLs_dig = VICLs_samp / q;

figure;
stem(IBLs_dig);
grid;
legend('IB digitalizado');
title('Final da digitalização');
xlabel('amostra k');
ylabel ('IB(k)');

YModA = fourierFilter(IALs_dig, samplesPerCycle);
YModB = fourierFilter(IBLs_dig, samplesPerCycle);
YModC = fourierFilter(ICLs_dig, samplesPerCycle);

spikePositionA = findElectricCurrentSpike(YModA, 2.5*nominalCurrent*FMI/q);
spikePositionB = findElectricCurrentSpike(YModB, 2.5*nominalCurrent*FMI/q);
spikePositionC = findElectricCurrentSpike(YModC, 2.5*nominalCurrent*FMI/q);

currentShutdownPositionA = findShutdownPosition(IALs_dig, samplesPerCycle);
currentShutdownPositionB = findShutdownPosition(IBLs_dig, samplesPerCycle);
currentShutdownPositionC = findShutdownPosition(ICLs_dig, samplesPerCycle);

voltageShutdownPositionA = findShutdownPosition(VANs_dig, samplesPerCycle);
voltageShutdownPositionB = findShutdownPosition(VBNs_dig, samplesPerCycle);
voltageShutdownPositionC = findShutdownPosition(VCNs_dig, samplesPerCycle);

shortCircuitInstant = ta_samp(min([spikePositionA, spikePositionB, spikePositionC]))
currentShutdownInstant = ta_samp(min([currentShutdownPositionA, currentShutdownPositionB, currentShutdownPositionC]))
voltageShutdownInstant = ta_samp(min([voltageShutdownPositionA, voltageShutdownPositionB, voltageShutdownPositionC]))


figure;
plot(ta_samp,IALs_dig,ta_samp,YModA,'r*');
legend('IALs_dig', 'Valor eficaz');
title('Amostragem das informações - Fase A');
xlabel('Tempo [s]');
ylabel ('Tensao [V]');
grid;

figure;
plot(ta_samp,IBLs_dig,ta_samp,YModB,'g*');
legend('IBLs_dig', 'Valor eficaz');
title('Amostragem das informações - Fase B');
xlabel('Tempo [s]');
ylabel ('Tensao [V]');
grid;

figure;
plot(ta_samp,ICLs_dig,ta_samp,YModC,'b*');
legend('ICLs_dig', 'Valor eficaz');
title('Amostragem das informações - Fase C');
xlabel('Tempo [s]');
ylabel ('Tensao [V]');
grid;
