% This is the main file for Mie Theory
% We want to compute extinction cross-section for infinite cylinders

%% Clear workspace and command window
clear; clc;

%% Set parameters
omega_p = 8.65e15;          % Plasma frequency
gamma = 0.01 * omega_p;     % Damping constant
v_F = 1.07e6;               % Fermi velocity
epsilon_m = 1;              % Relative permittivity of the outside
a = 2e-9;                   % Radius of the cylinder
c = 3e8;                    % Speed of light in vacuum

%% Calculate dielectric function
% transverse dielectric function
epsilon_T_val = @(omega) (1 - (omega_p^2)./(omega.^2 + 1i*gamma*omega)); 
% longitudinal dielectric function
epsilon_L_val = @(k, omega) (1 - (omega_p^2)./ ...
                       (omega.^2 + 1i*gamma*omega - (0.6*v_F^2)*k.^2)); 

%% Frequency and wavevector
SamplingPoints = 1000;
omega = linspace(0.4*omega_p, 1.4*omega_p, SamplingPoints).';           % Frequency range
% omega = 0.8 * omega_p;                                                  % Single frequency for testing
epsilon_T = epsilon_T_val(omega);
k_L = sqrt((omega.^2 + 1i*gamma*omega - omega_p^2) / (0.6*v_F^2));      % Longitudinal wavevector inside
k_T = sqrt(epsilon_T) .* omega / c;                                     % Transverse wavevector inside
k_o = sqrt(epsilon_m) * omega / c;                                      % Wavevector of the outside

%% Calculate scattering coefficients
N = 20;                     % Maximum order of scattering coefficients
a_n_local = ScattingCoefficients_Local(N, k_o, k_T, a, epsilon_T, epsilon_m);
a_n_nonlocal = ScattingCoefficients_Nonlocal(N, k_o, k_T, k_L, a, epsilon_T, epsilon_m);

%% Calulate extinction cross-section
C_ext_local = ExtinctionCrossSection(k_o, a, a_n_local);
C_ext_nonlocal = ExtinctionCrossSection(k_o, a, a_n_nonlocal);

%% Draw results
figure;
x_axis = linspace(0.4, 1.4, SamplingPoints).';
plot(x_axis, log(C_ext_local), 'LineWidth', 2, 'LineStyle', '--');
hold on;
plot(x_axis, log(C_ext_nonlocal), 'LineWidth', 2);
legend('Local Response', 'Nonlocal Response', 'FontSize', 12);
xlabel('$\omega/\omega_p$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\log(C_{ext})$', 'Interpreter', 'latex', 'FontSize', 14);
title(['Na cylinder of radius ' num2str(a*1e9) ' nm, in a medium with $\varepsilon_m =$ ' num2str(epsilon_m)], 'Interpreter', 'latex', 'FontSize', 16);
grid on;