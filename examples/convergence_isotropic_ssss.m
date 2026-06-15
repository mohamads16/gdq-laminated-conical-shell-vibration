%CONVERGENCE_ISOTROPIC_SSSS  Grid convergence for selected Table 1 cases.

clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));

R2 = 1.0;
alpha = deg2rad(30);
L = 0.75 / sin(alpha);
R1 = R2 - L*sin(alpha);
h = 0.01 * R2;
lam = laminate_isotropic(1, 0.3, 1, h);

Ns = [5 7 9 11 13 17 21 33 65 99];
nCases = [0 6];
lambdaRef = [0.1493 0.2542];

modeN = zeros(numel(Ns)*numel(nCases),1);
gridN = zeros(size(modeN));
lambda = zeros(size(modeN));
reference = zeros(size(modeN));

row = 0;
for q = 1:numel(nCases)
    for j = 1:numel(Ns)
        row = row + 1;
        params = struct('L',L, 'R1',R1, 'alpha',alpha, 'n',nCases(q), ...
                        'rho',lam.rho, 'h',lam.h, ...
                        'A',lam.A, 'B',lam.B, 'D',lam.D, 'Nx',Ns(j));
        sol = solve_conical_shell_modes(build_conical_shell_gdq_model(params), 'SS', 'SS', 3);
        modeN(row) = nCases(q);
        gridN(row) = Ns(j);
        lambda(row) = sol.lambda(1);
        reference(row) = lambdaRef(q);
    end
end

absError = abs(lambda - reference);
results = table(modeN, gridN, lambda, reference, absError, ...
    'VariableNames', {'n','Nx','lambda','lambda_ref','abs_error'});
disp(results);

