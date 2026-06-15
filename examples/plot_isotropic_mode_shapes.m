%PLOT_ISOTROPIC_MODE_SHAPES  Plot 1D and 3D mode shapes for one isotropic case.

clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));

Nx = 33;
R2 = 1.0;
alpha = deg2rad(30);
L = 0.75 / sin(alpha);
R1 = R2 - L*sin(alpha);
h = 0.01 * R2;
n = 3;

lam = laminate_isotropic(1, 0.3, 1, h);
params = struct('L',L, 'R1',R1, 'alpha',alpha, 'n',n, ...
                'rho',lam.rho, 'h',lam.h, ...
                'A',lam.A, 'B',lam.B, 'D',lam.D, 'Nx',Nx);

mdl = build_conical_shell_gdq_model(params);
sol = solve_conical_shell_modes(mdl, 'SS', 'SS', 3);

fprintf('First three lambda values for n=%d: ', n);
fprintf('%.6f ', sol.lambda);
fprintf('\n');

figure('Name', '1D mode shapes');
plot_mode_shapes_1d(mdl, sol, 1:3);

plot_mode_shape_3d(mdl, sol, 1);

