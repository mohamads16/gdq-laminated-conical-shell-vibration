%VALIDATE_TABLE1_ISOTROPIC  Reproduce Shu (1996) Table 1 isotropic cases.

clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));

Nx = 99;
R2 = 1.0;
alpha = deg2rad(30);
L = 0.75 / sin(alpha);
R1 = R2 - L*sin(alpha);
h = 0.01 * R2;

lam = laminate_isotropic(1, 0.3, 1, h);
bcs = {'SS','SS'; 'SS','C'; 'C','SS'; 'C','C'};
refs = table1_reference();

bcName = strings(40,1);
nValue = zeros(40,1);
lambda = zeros(40,1);
lambdaRef = zeros(40,1);

row = 0;
for q = 1:size(bcs,1)
    for n = 0:9
        row = row + 1;
        params = struct('L',L, 'R1',R1, 'alpha',alpha, 'n',n, ...
                        'rho',lam.rho, 'h',lam.h, ...
                        'A',lam.A, 'B',lam.B, 'D',lam.D, 'Nx',Nx);

        mdl = build_conical_shell_gdq_model(params);
        sol = solve_conical_shell_modes(mdl, bcs{q,1}, bcs{q,2}, 3);

        bcName(row) = bcs{q,1} + "-" + bcs{q,2};
        nValue(row) = n;
        lambda(row) = sol.lambda(1);
        lambdaRef(row) = refs(q,n+1);
    end
end

absError = abs(lambda - lambdaRef);
pctError = 100*(lambda - lambdaRef)./lambdaRef;
results = table(bcName, nValue, lambda, lambdaRef, absError, pctError, ...
    'VariableNames', {'BC','n','lambda','lambda_ref','abs_error','pct_error'});

disp(results);
fprintf('\nMax absolute error by boundary condition (Nx=%d):\n', Nx);
for q = 1:size(bcs,1)
    name = bcs{q,1} + "-" + bcs{q,2};
    mask = bcName == name;
    fprintf('  %-5s %.6g\n', name, max(absError(mask)));
end

function refs = table1_reference()
refs = [
    0.1493 0.6946 0.3593 0.2117 0.1907 0.2161 0.2542 0.2977 0.3516 0.4138
    0.9749 0.7451 0.4506 0.2976 0.2522 0.2640 0.2935 0.3349 0.3869 0.4485
    0.6018 0.7185 0.3763 0.2324 0.1977 0.2166 0.2525 0.2977 0.3516 0.4138
    0.9886 0.7525 0.4659 0.3153 0.2606 0.2647 0.2935 0.3349 0.3869 0.4485
];
end

