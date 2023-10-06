% ECE 345 Fall 2023
% Author: Shrey Dobariya (smd403)
%Exercise 1
sample_rate = 1e5;
set(groot, 'defaultAxesFontSize',16);
set(groot, 'defaultLineLineWidth', 2);

% Part A Cosine DT PLOT of a(n)
time = -4:11;
sigfunc_a = @(n) (1.5*cos((3*pi/8) * n));
stem(time, sigfunc_a(time) );
title('(a) Plot of a[n]');
xlabel('time (n)');
ylabel('amplitude');

% Part B Power DT plot of b(n)
time = 0:8;
sigfunc_b = @(n) (3*0.8.^(n-3));
stem(time, sigfunc_b(time) );
title('(b) Plot of b[n]');
xlabel('time (n)');
ylabel('amplitude');

% Part C DT Plot of c(n)
time = 0:15;
sigfun_c = @(n) (sigfunc_a(n) + sigfunc_b(n));
stem(time, sigfun_c(time) );
title('(c) Plot of c[n]');
xlabel('time (n)');
ylabel('amplitude');


% Part D
sim_sample_rate = 4e4;
time = 0:(1/sim_sample_rate):0.005;
sigfun_x = @(t) (1000*cos(880*pi * t));
plot(time, sigfun_x(time) );
title('(d) Plot of x[t]');
xlabel('time (t)');
ylabel('amplitude');

% Part E
sim_sample_rate = 50;
time = 0:(1/sim_sample_rate):4;
signal_function_y = @(t) (1 <= t & t <= 2) .* t + (t < 1 | t > 2) .* 0;
plot(time, signal_function_y(time) );
title('(e) Plot of y[t]');
xlabel('time (t)');
ylabel('amplitude');

% Part F
sim_sample_rate = 1.2e3;
time = -3:(1/sim_sample_rate):3;
signal_function_z = @(t) (1.2.^t .* sin(10*pi*t));
plot(time, signal_function_z(time) );
title('(f) Plot of z[t]');
xlabel('time (t)');
ylabel('amplitude');
%Exercise 2
% Part A
time = 0:15;
window_2a = (time >= 5) .* (time <= 10);
clip = sigfun_c(time) .* window_2a;
stem(clip)
title('Window plot ofc[n]');
xlabel('time (t)');
ylabel('amplitude');

% Part B
sim_sample_rate = 4e4;
time = 0:(1/sim_sample_rate):0.005;
window_x = (time >= 0.002) .* (time <= 0.003);
clip = sigfun_x(time) .* window_x;
plot(clip)
title('Window plot ofx[t]');
xlabel('time (t)');
ylabel('amplitude');

% Part C
sim_sample_rate = 1.2e3;
time = -3:(1/sim_sample_rate):3;
window_z = (time >= 0) .* (time <= 2);
clip = signal_function_z(time) .* window_z;
plot(clip)
title('Window plot ofz[t]');
xlabel('time (t)');
ylabel('amplitude');

% Exercise 3
%%%
% Part A
sim_sample_rate = 5e3;
time = 0:1/(sim_sample_rate):25;
part_a_signal = @(t) ((0.9).^(t) .* sin(10*pi.*t));
maxValue = max(part_a_signal(time));
maxIndices = find(round(part_a_signal(time), 4) == round(maxValue,4));
minValue = min(part_a_signal(time));
minIndices = find(round(part_a_signal(time), 4) == round(minValue, 4));
% We will do a rough approximation here  because some points are 0 but don't equal 0 exactly
zeroes = find(part_a_signal(time) < 0.01 & part_a_signal(time) > -0.01);
plot(time, part_a_signal(time));
hold on;
plot(time(maxIndices), maxValue*ones(size(maxIndices)), 'xr');
plot(time(minIndices), minValue *ones(size(minIndices)), 'xg');
plot(time(zeroes), zeros(size(zeroes)), 'o');
grid on;
hold off;
legend('Signal', 'Max', 'Min', 'Zeroes');
title('Finding Peaks On Plot A');
xlabel('time (t)');
ylabel('a(t)');


% Part B
sim_sample_rate = 5e3;
time = -2:1/(sim_sample_rate):2;
part_b_signal = @(t) (sin(10*pi.*t)+cos(8.*t));

maxValue = max(part_b_signal(time));
%rough approximation due to some rounding issues
maxIndices = find(round(part_b_signal(time), 4) == round(maxValue,4));

minValue = min(part_b_signal(time));
%rough approximation due to some rounding issues
minIndices = find(round(part_b_signal(time), 4) == round(minValue, 4));

% rough approximate because some points are 0 but don't equal 0 exactly
zeroes = find(part_b_signal(time) < 0.01 & part_b_signal(time) > -0.01);

plot(time, part_b_signal(time));
hold on;
plot(time(maxIndices), maxValue*ones(size(maxIndices)), 'xr');
plot(time(minIndices), minValue *ones(size(minIndices)), 'xg');
plot(time(zeroes), zeros(size(zeroes)), 'o');
grid on;
hold off;
legend('Signal', 'Max', 'Min', 'Zeroes');
title('Finding Peaks On Plot B');
xlabel('time (t)');
ylabel('b(t)');


% Part C
time = -5:5;
part_c_signal = @(t) (-1.1).^t;

maxValue = max(part_c_signal(time));
%rough approximation due to some rounding issues
maxIndices = find(round(part_c_signal(time), 4) == round(maxValue,4));

minValue = min(part_c_signal(time));
%rough approximation due to some rounding issues
minIndices = find(round(part_c_signal(time), 4) == round(minValue, 4));

stem(time, part_c_signal(time));
hold on;
plot(time(maxIndices), maxValue*ones(size(maxIndices)), 'xr');
plot(time(minIndices), minValue *ones(size(minIndices)), 'xg');
%this plot can't have zeroes so we just skip that code
grid on;
hold off;
legend('Signal', 'Max', 'Min');
title('Finding Peaks On Plot C');
xlabel('timet)');
ylabel('c[n]');


% Part D
time = -100:100;
part_d_signal = @(t) cos(3*pi/5*t) + cos(2 * pi/5 *(t-3));

maxValue = max(part_d_signal(time));
%rough approximation due to some rounding issues
maxIndices = find(round(part_d_signal(time), 4) == round(maxValue,4));

minValue = min(part_d_signal(time));
%rough approximation due to some rounding issues
minIndices = find(round(part_d_signal(time), 4) == round(minValue, 4));

% rough approximate because some points are 0 but don't equal 0 exactly
zeroes = find(part_d_signal(time) < 0.0001 & part_d_signal(time) > -0.0001);

stem(time, part_d_signal(time));
hold on;
plot(time(maxIndices), maxValue*ones(size(maxIndices)), 'xr');
plot(time(minIndices), minValue *ones(size(minIndices)), 'xg');
plot(time(zeroes), zeros(size(zeroes)), 'o');
grid on;
hold off;
legend('Signal', 'Max', 'Min', 'Zeroes');
title('Finding Peaks On Plot D');
xlabel('time (t)');
ylabel('d[n]');



% Exercise 4
% for each of them I'll show both the methods , the asked one and the
% riemann sum so that we can see how accuate the values are 
%Part A
partA = @(t) cos(120*pi*t).^2;
n = 1000;
a = 0;
b = 1;
width = (b - a) / n;
ri_sum = 0;

for i = 1:n
    m = (i - 0.5)*width;
    ri_sum = ri_sum + partA(m)*width;
end
fprintf('The value of the integral for part A from %d to %d : %.3f\n',a,b,ri_sum);
% partA method2 
t = 0:Deltasim:1;
x = cos(120*pi*t).^2;
A = sum(x) * Deltasim;
fprintf('The approximate value of integral A is: %f\n', A);

%Part B
partB = @(t) (cos(120*pi*t) .* sin(120*pi*t));
n = 1000;
a = 0;
b = 1;
width = (b - a) / n;
ri_sum = 0;

for i = 1:n
    m = (i - 0.5)*width;
    ri_sum = ri_sum + partB(m)*width;
end
fprintf('The value of the integral for part B from %d to %d: %.8f\n',a,b,ri_sum);
%partB method2
fsim = 1e5;
Deltasim = 1 / fsim;
t_b = 0:Deltasim:1;
x_b = cos(120 * pi * t_b) .* sin(120 * pi * t_b);
area_b = sum(x_b) * Deltasim;
fprintf('The approximate value of integral B is: %f\n', area_b);

%Part C
e = exp(1);
partC = @(t) (e.^(-0.5*t)+ 2*e.^(-t)).*cos(60*pi*t);
n = 1000;
a = 3;
b = 5;
width = (b - a) / n;
ri_sum = 0;
for i = 1:n
    m = (i - 0.5)*width;
    ri_sum = ri_sum + partC(m)*width;
end
fprintf('The value of the integral for part C from %d to %d: %.8f\n',a,b,ri_sum);
% partC method2
fsim = 1e5;
Deltasim = 1 / fsim;
t_c = 3:Deltasim:5;
x_c = (exp(-0.5*t_c) + 2 * exp(-t_c)) .* cos(60 * pi * t_c);
area_c = sum(x_c) * Deltasim;
fprintf('The value of the integral for part C: %.8f\n',area_c);

%Part D
partD = @(t) cos(6*pi*t).*cos(6*pi*(t - 1/5));
n = 1000;
a = 0;
b = 2;
width = (b - a) / n;
ri_sum = 0;

for i = 1:n
    m = (i - 0.5)*width;
    ri_sum = ri_sum + partD(m)*width;
end
fprintf('The value of the integral for part D from  %d to %d: %.8f\n',a,b,ri_sum);
%partD method2
fsim = 1e5;
Deltasim = 1 / fsim;
t_d = 0:Deltasim:2;
x_d = cos(6 * pi * t_d) .* cos(6 * pi * (t_d - 1/5));
area_d = sum(x_d) * Deltasim;
fprintf('The value of the integral for part d: %.8f\n',area_d);


%Exercise 5
file= 'whalesong.mp3';
[audio, fs] = audioread(file); 
fs_kHz = fs / 1000;
clip_length = length(audio) / fs_kHz;
fprintf('Sampling Rate in kHz is:%0.6f\n',fs_kHz);
fprintf('length of the clip  is:%0.6f\n',clip_length);

t = (0:length(audio) - 1) /fs;
plot(t, audio);
xlabel('Time (s)');
ylabel('Amplitude');
title('Audio Signal');
[max_magnitude, max_index] = max(abs(audio));
max_time = t(max_index);
fprintf('maximum magnitude  is:%0.6f\n',max_magnitude);
fprintf('Timestamp of max  is:%0.6f\n',max_time);
% extracting and plotting whale call between 31 and 33 seconds
startInd = round(31* fs)+1;
endInd = round(33*fs);
whcallAUDIO = audio(startInd:endInd);
whaleCallTime = (0:length(whcallAUDIO) - 1) / fs;
soundsc(whcallAUDIO);
plot(whaleCallTime, whcallAUDIO);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('31s to 33s Whale Call');
% Calculate the energy of the 2-second clip
totalEnergy = trapz(abs(whcallAUDIO).^2);
fprintf('Energy in 2s Clip  is:%0.6f\n',totalEnergy );







