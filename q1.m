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


%% Part C: Phase Velocity and Propagation
c0 = 3e8;               % velocity of light in vacum (m/s)

A = 2; B = 1;           % Example constants (A+B > 1)
n_op = A+B;
v_p = c0 / n_op;        % Phase velocity

% Distances
z = [0, v_p * 10e-9, v_p * 20e-9]; % z1, z2, z3

% Expected result:
% signal(z,t) = Re{a(t-beta_1*z) * e^(j*(w*t-beta_0*z))}
% beta_0 is 0 so we should get delay (of the Gaussian envelope) only.
% beta = k_0 * n = omega / c0 * n
% beta_1 = beta / omega = n / c0

% the delay should be:
% beta_1 * z = n_op / c0 * v_p * 10^-8 * (i-1) = 10^-8 * (i-1)
figure('Name', 'Part C: Constant Refractive Index');
for i = 1:3
    z_curr = z(i);
    t_signal_at_z = propagateWaveInZ(f, f_signal, n_op, z_curr);

    % Plotting (taking the real part since the signal is physical)
    subplot(3, 1, i);
    plot(t * 1e9, real(t_signal_at_z));
    title(['Signal at z = ', num2str(z_curr, '%.2f'), ' meters']);
    xlabel('Time [nsec]'); ylabel('Amplitude');
    grid on;
end

%% Part D: Quadratic Dispersion
delta_n = 10 * A * (f - f0) * dt;
n_f = n_op + delta_n;

% Expected result:
% beta(omega) = k_0 * n = omega / c0 * n
% beta is polynomial of order 2 and thus the expected result is chage
% the width of the pulse and also phase-chage:

% tau(z) = tau(0)*sqrt{1+(beta_2*z/tau_0^2)^2}

% beta_2 = d^2beta / domega^2 = 10 * A * dt * 2 *pi
% beta_2 * z = 10 * A * dt * 2 * pi * v_p * 10^-8 * (i-1)
% beta_2 * z = 10^-17 * 4 * pi * (i-1) * v_p
figure('Name', 'Part D: Quadratic Dispersion');
for i = 1:3
    z_curr = z(i);
    t_signal_at_z = propagateWaveInZ(f, f_signal, n_f, z_curr);

    % Plotting (taking the real part since the signal is physical)
    subplot(3, 1, i);
    plot(t * 1e9, real(t_signal_at_z));
    title(['Signal at z = ', num2str(z_curr, '%.2f'), ' meters']);
    xlabel('Time [nsec]'); ylabel('Amplitude');
    grid on;
end

%% Part E: Quadratic Refractive Index
delta_n = 3000 * B * (f - f0).^2 * dt^2;
n_f = n_op + delta_n;

% Expected result:
% beta(omega) = k_0 * n = omega / c0 * n
% beta is polynomial of order 3 and we did not calculate that in class.

figure('Name', 'Part E: Quadratic Refractive Index');
for i = 1:3
    z_curr = z(i);
    t_signal_at_z = propagateWaveInZ(f, f_signal, n_f, z_curr);

    % Plotting (taking the real part since the signal is physical)
    subplot(3, 1, i);
    plot(t * 1e9, real(t_signal_at_z));
    title(['Signal at z = ', num2str(z_curr, '%.2f'), ' meters']);
    xlabel('Time [nsec]'); ylabel('Amplitude');
    grid on;
end
