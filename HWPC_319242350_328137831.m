close all;
clc;

%% Question 1

%% Utility functions
function signal = createGaussianEnvelope(t_vector, t_peak, fwhm)
    % Calculate sigma from FWHM
    % Formula: FWHM = 2 * sqrt(2 * ln(2)) * sigma
    sigma = fwhm / (2 * sqrt(2 * log(2))); 
    
    % Calculate Gaussian Envelope
    envelope = exp(-(t_vector - t_peak).^2 / (2 * sigma^2));
    signal = envelope;
end

function signal = createPulse(t_vector, t_peak, fwhm, f_carrier_frequency)
    % Calculate Gaussian Envelope
    envelope = createGaussianEnvelope(t_vector, t_peak, fwhm); % a(t)
    carrier = cos(2 * pi * f_carrier_frequency * t_vector); % Re{e^(j*w*t)}
    signal = envelope .* carrier; % Re{a(t) * e^(j*w*t)}
end

function signal_at_z = propagateWaveInZ(frequencies, fourier_signal, n, z)
    c0 = 3e8; % velocity of light in vacum (m/s)
    omega = 2 * pi * frequencies;
    k_0 = omega / c0;
    beta = k_0 .* n;
    
    % 2. Enforce Physical Antisymmetry (beta(-f) = -beta(f))
    % For real-valued time signals, the phase shift must be an odd function.
    % We take the beta calculated for positive frequencies and mirror it
    % for the negative side with a sign flip.
    beta_antisymmetric = sign(frequencies) .* interp1(frequencies, beta, abs(frequencies), 'linear', 'extrap');
    
    % 3. Apply the Phase Shift (Propagation)
    % H(f) = exp(-j * beta * z)
    propagation_factor = exp(-1j * beta_antisymmetric * z);

    f_signal_propagated = fourier_signal .* propagation_factor;
    
    % Return to time domain using Inverse FFT
    signal_at_z = ifft(ifftshift(f_signal_propagated));
end

function f_filtered_signal = SingleSidebandFilter(frequencies, fourier_signal, cutoff)
    % SINGLESIDEBANDFILTER Filters the spectrum to keep only frequencies above a specific value.
    % This creates a single-sided spectrum by removing negative frequencies 
    % and components below the cutoff, often used for analytic signal processing.
    %
    % frequencies    : The frequency vector [Hz] (from fftshift)
    % fourier_signal : The shifted FFT of the input signal
    % cutoff         : The lower frequency boundary [Hz]. Components below this are zeroed.

    % 1. Create a logical mask to isolate frequencies strictly greater than or equal to the cutoff.
    mask = frequencies >= cutoff;
    
    % 2. Apply the mask to the fourier signal
    f_filtered_signal = fourier_signal .* mask;
end

function plotPropagatedSignalsInZ(graph_name, t_vector, z_vals, frequencies, f_filtered_signal, n)
    figure('Name', graph_name);
    total_signal = zeros(size(t_vector)); % Initialize vector for superposition
    colors = {'b', [0 0.5 0], 'm'}; % Blue, Dark Green, Magenta
    
    plots_num = length(z_vals)+1;
    for i = 1:plots_num-1
        z_curr = z_vals(i);
        % Propagate the filtered signal to the current distance
        t_signal_at_z = propagateWaveInZ(frequencies, f_filtered_signal, n, z_curr);
        
        % Accumulate for superposition
        current_real_sig = real(t_signal_at_z);
        total_signal = total_signal + current_real_sig;
        
        % Plot individual signals with unique colors
        subplot(plots_num, 1, i);
        plot(t_vector * 1e9, current_real_sig, 'Color', colors{i}, 'LineWidth', 1);
        title(['Individual Signal at z = ', num2str(z_curr*1e3, '%.2f'), ' mm']);
        xlabel('Time [nsec]'); ylabel('Amp');
        grid on;
    end
    
    % Plot the Superposition (The sum of all three waves) in Red
    subplot(plots_num, 1, plots_num);
    plot(t_vector * 1e9, total_signal, 'Color', 'r', 'LineWidth', 1.5);
    title('Superposition of all signals');
    xlabel('Time [nsec]'); ylabel('Amp');
    grid on;
end

%% Part A: Signal Generation
dt = 100e-12;               % 100 psec to seconds
t = 0:dt:0.1e-6;            % Time vector from 0 to 0.1 microsec

% Gaussian Pulse Parameters
t_center = 50e-9;           % Center at 50 nsec
fwhm = 4.2e-9;              % 4.2 nsec Full Width at Half Maximum

% Carrier Parameters
f0 = 1e9;                   % 1 GHz carrier frequency

t_signal = createPulse(t,t_center,fwhm,f0);

%% Part B: Fourier Transform
f_signal = fftshift(fft(t_signal));
f = linspace(-1/(2*dt), 1/(2*dt), length(t_signal));

% Plotting the Spectrum
figure;
stem(f / 1e9, abs(f_signal));
title('Spectrum of Gaussian Pulse');
xlabel('Frequency [GHz]'); ylabel('|X(f)|');
grid on;

% Main Frequencies of the signal are 1[GHz] and -1[GHz].
% We shall choose positive frequencies only and use the real part of the
% result (symetry)
% we might get the same result without passing through positive filter
% with the constraint beta(-omega) = -beta(omega)
f_filtered_signal = SingleSidebandFilter(f, f_signal, 0);

%% Part C: Phase Velocity and Propagation
c0 = 3e8;               % velocity of light in vacum (m/s)

A = 29; B = 1;           % Example constants (A+B > 1)
n_op = A+B;
v_p = c0 / n_op;        % Phase velocity
v_g = v_p;              % Linear Disperssion (Group Velocity = Phase velocity)

% Distances
z = [0, v_g * 10e-9, v_g * 20e-9]; % z1, z2, z3

% Expected result:
% signal(z,t) = Re{a(t-beta_1*z) * e^(j*(w*t-beta_0*z))}
% beta = k_0 * n = omega / c0 * n
% beta_1 = dbeta / domega = n / c0

% the delay should be:
% beta_1 * z = n_op / c0 * v_p * 10^-8 * (i-1) = 10^-8 * (i-1)
% the delay expected to be 10 [nsec].

plotPropagatedSignalsInZ('Part C: Constant Refractive Index', t, z, f, f_filtered_signal, n_op)

%% Part D: Quadratic Dispersion
delta_n = 10 * A * (f - f0) * dt;
n_f = n_op + delta_n;

beta_1 = (n_op + 10 * A * dt * f0) / c0;
v_g = 1 / beta_1;        % Group Velocity

% Distances
z = [0, v_g * 10e-9, v_g * 20e-9]; % z1, z2, z3

% Expected result:
% beta(omega) = k_0 * n = omega / c0 * n
% beta is polynomial of order 2 and thus the expected result is chage
% the width of the pulse and also phase-chage:

% tau(z) = tau(0)*sqrt{1+(beta_2*z/tau_0^2)^2}
plotPropagatedSignalsInZ('Part D: Quadratic Dispersion', t, z, f, f_filtered_signal, n_f)

%% Part E: Quadratic Refractive Index
delta_n = 3000 * B * (f - f0).^2 * dt^2;
n_f = n_op + delta_n;

beta_1 = n_op / c0;
v_g = 1 / beta_1;        % Group Velocity

% Distances
z = [0, v_g * 10e-9, v_g * 20e-9]; % z1, z2, z3

% Expected result:
% beta(omega) = k_0 * n = omega / c0 * n
% beta is polynomial of order 3 and we did not calculate that in class.

plotPropagatedSignalsInZ('Part E: Quadratic Refractive Index', t, z, f, f_filtered_signal, n_f)

%% Question 2e1

lambda = 1.3:0.005:1.8;  
n1 = 1;
n3 = 1.45;
n2 = sqrt(n1*n3);
h0 = 1.55/(4*n2);

r1_2 = (n1-n2)/(n1+n2);
r2_3 = (n2-n3)/(n2+n3);

delta = 2*(2*pi./lambda)*n2*h0;
r_tot = (r2_3*exp(-1i*delta/2)+r1_2*exp(1i*delta/2))./(exp(1i*delta/2)+r2_3*r1_2*exp(-1i*delta/2));  
R_tot = ((abs(r_tot)).^2)*100;  

figure()
plot(lambda,R_tot,linewidth=1.5);
grid minor
xlabel('Wavelength [micron]');
ylabel('Power Reflection Coefficient R [%]');
title({'Power Reflection Coefficient vs. Wavelength'});


%% Question 2e2

lambda0 = 1.55;
degree1 = 0:0.1:90;  
degree2 = asind(n1*(sind(degree1))/n2);
degree3 = asind(n1*(sind(degree1))/n3);

delta_v2 = 2*(2*pi/lambda0)*n2*h0*cosd(degree2);

r1_2_v2 = (n1*cosd(degree1)-n2*cosd(degree2))./(n1*cosd(degree1)+n2*cosd(degree2));
r2_3_v2 = (n2*cosd(degree2)-n3*cosd(degree3))./(n2*cosd(degree2)+n3*cosd(degree3));

r_tot_v2 = (r2_3_v2.*exp(-1i*delta_v2/2)+r1_2_v2.*exp(1i*delta_v2/2))./(exp(1i*delta_v2/2)+r2_3_v2.*r1_2_v2.*exp(-1i*delta_v2/2)); 
R_tot_v2 = ((abs(r_tot_v2)).^2)*100; 

figure()
plot(degree1,R_tot_v2,linewidth=1.5);
grid minor
xlabel('Incidence Angle [degrees]');
ylabel('Power Reflection Coefficient R [%]');
title({'Power Reflection Coefficient vs Incidence Angle'});


%% Question 2f

figure()
plot(lambda,R_tot,linewidth=1.5);
grid minor
xlabel('Wavelength [micron]');
ylabel('Power Reflection Coefficient R [%]');
title({'Power Reflection Coefficient vs. Wavelength'});
hold on

new_h0 = h0;
num_of_local_min = 1;

while (~(num_of_local_min>1))
    new_h0 = new_h0 + 2*h0;
    delta_v3 = 2*(2*pi./lambda)*n2*new_h0;
    r_tot_v3 = (r2_3*exp(-1i*delta_v3/2)+r1_2*exp(1i*delta_v3/2))./(exp(1i*delta_v3/2)+r2_3*r1_2*exp(-1i*delta_v3/2));  
    R_tot_v3 = ((abs(r_tot_v3)).^2)*100; 
    plot(lambda,R_tot_v3,linewidth=1.5);
    indices_of_local_min = islocalmin(R_tot_v3);
    num_of_local_min = sum(indices_of_local_min);
end


%% Question 3c1

lambda0 = 1.55;   
n_L = 1.45;      
n_H = 1.01*n_L;  

r_L_to_H = (n_L - n_H)/(n_L + n_H); 
N_periods = floor(1/(abs(r_L_to_H)));
T_tot = (1/(1-(r_L_to_H)^2))*[-(1+(r_L_to_H)^2), 2*r_L_to_H; 2*r_L_to_H, -(1+(r_L_to_H)^2)]; 
T_to_power_N = (T_tot)^(N_periods);  
amplitude_ref_coeff = -(T_to_power_N(2,1))/(T_to_power_N(2,2)); 
power_ref_coeff = (abs(amplitude_ref_coeff))^2; 

while power_ref_coeff <= 0.99
    N_periods = N_periods + 1;
    T_to_power_N = (T_tot)*(T_to_power_N);
    amplitude_ref_coeff = -(T_to_power_N(2,1))/(T_to_power_N(2,2));
    power_ref_coeff = (abs(amplitude_ref_coeff))^2;
end

lambda_range = 1.45:0.000005:1.65;  

phase_amp_ref_coeff = []; 
power_ref_coeff = [];     

for current_lambda=lambda_range
    delta_L = pi*lambda0/current_lambda;  
    delta_H = pi*lambda0/current_lambda;  
    % calculating transfer matrix elements
    A = (1/(1-(r_L_to_H)^2))*exp(-1i*delta_L/2)*(exp(-1i*delta_H/2)-(exp(1i*delta_H/2))*((r_L_to_H)^2));
    B = (1/(1-(r_L_to_H)^2))*2*1i*r_L_to_H*(exp(-1i*delta_L/2))*sin(delta_H/2);
    C = (1/(1-(r_L_to_H)^2))*(-2)*1i*r_L_to_H*(exp(1i*delta_L/2))*sin(delta_H/2);
    D = (1/(1-(r_L_to_H)^2))*exp(1i*delta_L/2)*(exp(1i*delta_H/2)-(exp(-1i*delta_H/2))*((r_L_to_H)^2));
    T_tot = [A, B; C, D];              
    T_to_power_N = (T_tot)^(N_periods);   
    amplitude_ref_coeff = -(T_to_power_N(2,1))/(T_to_power_N(2,2));
    current_phase = angle(amplitude_ref_coeff);
    phase_amp_ref_coeff = [phase_amp_ref_coeff current_phase];
    current_power_ref_coeff = (abs(amplitude_ref_coeff))^2;
    power_ref_coeff = [power_ref_coeff current_power_ref_coeff];
end

figure() 
plot(lambda_range, unwrap(phase_amp_ref_coeff), linewidth=1.5);
grid minor
xlabel('Wavelength [\mum]');
ylabel('Phase of Amplitude Reflection Coefficient r-DBR [rad]');
title({'Phase of Amplitude Reflection Coefficient r-DBR vs Wavelength'});
ylim([-0.5 7]);

figure()
plot(lambda_range, power_ref_coeff*100, linewidth=1.5);
grid minor
xlabel('Wavelength [\mum]');
ylabel('Power Reflection Coefficient R [%]');
title({'Power Reflection Coefficient R vs Wavelength'});


%% Question 3c2

lambda0 = 1.55;  
n_L = 1.45;     
n_H = 1.1*n_L;   

r_L_to_H = (n_L - n_H)/(n_L + n_H);
N_periods = floor(1/(abs(r_L_to_H)));
T_tot = (1/(1-(r_L_to_H)^2))*[-(1+(r_L_to_H)^2), 2*r_L_to_H; 2*r_L_to_H, -(1+(r_L_to_H)^2)];
T_to_power_N = (T_tot)^(N_periods); 
amplitude_ref_coeff = -(T_to_power_N(2,1))/(T_to_power_N(2,2)); 
power_ref_coeff = (abs(amplitude_ref_coeff))^2;

while power_ref_coeff <= 0.99
    N_periods = N_periods + 1;
    T_to_power_N = (T_tot)*(T_to_power_N);
    amplitude_ref_coeff = -(T_to_power_N(2,1))/(T_to_power_N(2,2));
    power_ref_coeff = (abs(amplitude_ref_coeff))^2;
end

lambda_range = 1.45:0.000005:1.65; 

phase_amp_ref_coeff = [];
power_ref_coeff = [];    

for current_lambda=lambda_range
    delta_L = pi*lambda0/current_lambda; 
    delta_H = pi*lambda0/current_lambda;  
    % calculating transfer matrix elements
    A = (1/(1-(r_L_to_H)^2))*exp(-1i*delta_L/2)*(exp(-1i*delta_H/2)-(exp(1i*delta_H/2))*((r_L_to_H)^2));
    B = (1/(1-(r_L_to_H)^2))*2*1i*r_L_to_H*(exp(-1i*delta_L/2))*sin(delta_H/2);
    C = (1/(1-(r_L_to_H)^2))*(-2)*1i*r_L_to_H*(exp(1i*delta_L/2))*sin(delta_H/2);
    D = (1/(1-(r_L_to_H)^2))*exp(1i*delta_L/2)*(exp(1i*delta_H/2)-(exp(-1i*delta_H/2))*((r_L_to_H)^2));
    T_tot = [A, B; C, D];               
    T_to_power_N = (T_tot)^(N_periods);  
    amplitude_ref_coeff = -(T_to_power_N(2,1))/(T_to_power_N(2,2));
    current_phase = angle(amplitude_ref_coeff);
    phase_amp_ref_coeff = [phase_amp_ref_coeff current_phase];
    current_power_ref_coeff = (abs(amplitude_ref_coeff))^2;
    power_ref_coeff = [power_ref_coeff current_power_ref_coeff];
end

figure() 
plot(lambda_range, unwrap(phase_amp_ref_coeff), linewidth=1.5);
grid minor
xlabel('Wavelength [\mum]');
ylabel('Phase of Amplitude Reflection Coefficient r-DBR [rad]');
title({'Phase of Amplitude Reflection Coefficient r-DBR vs Wavelength'});
ylim([-0.5 7]);

figure()
plot(lambda_range, power_ref_coeff*100, linewidth=1.5);
grid minor
xlabel('Wavelength [\mum]');
ylabel('Power Reflection Coefficient R [%]');
title({'Power Reflection Coefficient R vs Wavelength'});


%% Question 3d

lambda0 = 1.55;  
n_L = 1.45;      
n_H = 1.01*n_L;  

N_periods = 400; 
theta_L = deg2rad(5);  
theta_H = asin((n_L/n_H)*sin(theta_L)); 
r_L_to_H = (n_L*cos(theta_L) - n_H*cos(theta_H))/(n_L*cos(theta_L) + n_H*cos(theta_H)); 

lambda_range = 1.45:0.000005:1.65;
phase_amp_ref_coeff = []; 
power_ref_coeff = [];    

for current_lambda=lambda_range
    delta_L = (pi*lambda0/current_lambda)*cos(theta_L);
    delta_H = (pi*lambda0/current_lambda)*cos(theta_H);
    A = (1/(1-(r_L_to_H)^2))*exp(-1i*delta_L/2)*(exp(-1i*delta_H/2)-(exp(1i*delta_H/2))*((r_L_to_H)^2));
    B = (1/(1-(r_L_to_H)^2))*2*1i*r_L_to_H*(exp(-1i*delta_L/2))*sin(delta_H/2);
    C = (1/(1-(r_L_to_H)^2))*(-2)*1i*r_L_to_H*(exp(1i*delta_L/2))*sin(delta_H/2);
    D = (1/(1-(r_L_to_H)^2))*exp(1i*delta_L/2)*(exp(1i*delta_H/2)-(exp(-1i*delta_H/2))*((r_L_to_H)^2));
    T_tot = [A, B; C, D];                
    T_to_power_N = (T_tot)^(N_periods);  
    amplitude_ref_coeff = -(T_to_power_N(2,1))/(T_to_power_N(2,2));
    current_phase = angle(amplitude_ref_coeff);
    phase_amp_ref_coeff = [phase_amp_ref_coeff current_phase];
    current_power_ref_coeff = (abs(amplitude_ref_coeff))^2;
    power_ref_coeff = [power_ref_coeff current_power_ref_coeff];
end

figure() 
plot(lambda_range, unwrap(phase_amp_ref_coeff), linewidth=1.5);
grid minor
xlabel('Wavelength [\mum]');
ylabel('Phase of Amplitude Reflection Coefficient r-DBR [rad]');
title({'Phase of Amplitude Reflection Coefficient r-DBR vs Wavelength'});
ylim([-0.5 7]);

figure()
plot(lambda_range, power_ref_coeff*100, linewidth=1.5);
grid minor
xlabel('Wavelength [\mum]');
ylabel('Power Reflection Coefficient R [%]');
title({'Power Reflection Coefficient R vs Wavelength'});

