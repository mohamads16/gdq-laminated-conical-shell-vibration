clear; clc;
Nx = 100;

% --- geometry from Table 1 ---
R2    = 1.0;
alpha = deg2rad(30);
L     = 0.75 / sin(alpha);
R1    = R2 - L*sin(alpha);   % = 0.25
h     = 0.01 * R2;

% --- isotropic laminate (E, rho cancel in lambda) ---
E   = 1;
nu  = 0.3;
rho = 1;
lam = laminate_isotropic(E, nu, rho, h);

geom = struct('L',L, 'R1',R1, 'alpha',alpha);

% --- reference λ1 from Shu, Table 1, SS–SS ---
lambda_ref = [0.1493 0.6946 0.3593 0.2117 0.1907 ...
              0.2161 0.2542 0.2977 0.3516 0.4138]';

% --- compute your λ1 for n=0..9 ---
lambda_num = zeros(10,1);

for n = 0:9
    params = struct('L',    geom.L, ...
                    'R1',   geom.R1, ...
                    'alpha',geom.alpha, ...
                    'n',    n, ...
                    'rho',  lam.rho, ...
                    'h',    lam.h, ...
                    'A',    lam.A, ...
                    'B',    lam.B, ...
                    'D',    lam.D, ...
                    'Nx',   Nx);

    mdl = step1(params);
    sol = step2(mdl, 'SS', 'SS', 3);   % compute first 3 modes

    lambda_num(n+1) = sol.lambda(1);

    % 1D mode shapes for mode 1
    figure('Name', sprintf('1D mode shapes, n=%d', n));
    shape_modes_1d(mdl, sol, 1);

    % 3D view for a few selected n
    if n == 0 || n == 2
        shape_mode_3d(mdl, sol, 1);
    end
end

% --- percent error (signed) ---
pct_err = 100*(lambda_num - lambda_ref)./lambda_ref;

% --- table with row names n=0..9 ---
T3 = table(lambda_num, lambda_ref, pct_err, ...
    'VariableNames', {'lambda1_num','lambda1_ref','pct_error'});
T3.Properties.RowNames = compose('n=%d', 0:9);

disp(T3);
