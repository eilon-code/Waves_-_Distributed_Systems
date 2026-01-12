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
