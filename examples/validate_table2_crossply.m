%VALIDATE_TABLE2_CROSSPLY  Reproduce Shu (1996) Table 2 convergence cases.

clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));

R2 = 1.0;
alpha = deg2rad(30);
L = 0.25 / sin(alpha);
R1 = R2 - L*sin(alpha);
h = 0.01 * R2;
Ns = [5 7 9 11 13 21];
layerCases = [Inf 2];
caseNames = ["without coupling"; "maximum coupling"];
refs = [
    0.3158 0.2240 0.1961 0.1971 0.1976 0.1976
    0.3155 0.2226 0.1820 0.1787 0.1799 0.1799
];

caseName = strings(numel(layerCases)*numel(Ns),1);
layers = zeros(numel(caseName),1);
N = zeros(numel(caseName),1);
lambda = zeros(numel(caseName),1);
lambdaRef = zeros(numel(caseName),1);

row = 0;
for q = 1:numel(layerCases)
    lam = laminate_shu_crossply(h, layerCases(q));
    for j = 1:numel(Ns)
        row = row + 1;
        params = struct('L',L, 'R1',R1, 'alpha',alpha, 'n',0, ...
                        'rho',lam.rho, 'h',lam.h, ...
                        'A',lam.A, 'B',lam.B, 'D',lam.D, 'Nx',Ns(j));

        mdl = build_conical_shell_gdq_model(params);
        sol = solve_conical_shell_modes(mdl, 'SS', 'SS', 3);

        caseName(row) = caseNames(q);
        layers(row) = layerCases(q);
        N(row) = Ns(j);
        lambda(row) = sol.lambda(1);
        lambdaRef(row) = refs(q,j);
    end
end

absError = abs(lambda - lambdaRef);
pctError = 100*(lambda - lambdaRef)./lambdaRef;
results = table(caseName, layers, N, lambda, lambdaRef, absError, pctError, ...
    'VariableNames', {'case','layers','N','lambda','lambda_ref','abs_error','pct_error'});

disp(results);
fprintf('\nMax absolute error:\n');
for q = 1:numel(layerCases)
    mask = caseName == caseNames(q);
    fprintf('  %-18s %.6g\n', caseNames(q), max(absError(mask)));
end

