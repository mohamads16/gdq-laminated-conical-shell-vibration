function tests = test_paper_validation
tests = functiontests(localfunctions);
end

function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));
end

function testTable1IsotropicAtHighResolution(testCase)
[computed, reference] = compute_table1_isotropic(99);
verifyLessThanOrEqual(testCase, max(abs(computed - reference), [], 'all'), 2.0e-3);

% Regression check for the corrected mixed-boundary labels.
verifyEqual(testCase, computed(2,1), 0.9749, 'AbsTol', 1.0e-3); % SS-C
verifyEqual(testCase, computed(3,1), 0.6018, 'AbsTol', 1.0e-3); % C-SS
end

function testTable2CrossplyConvergence(testCase)
[computed, reference] = compute_table2_crossply();
verifyLessThanOrEqual(testCase, max(abs(computed(1,:) - reference(1,:))), 7.0e-4);
verifyLessThanOrEqual(testCase, max(abs(computed(2,:) - reference(2,:))), 7.0e-3);
verifyEqual(testCase, computed(2,5), 0.1799, 'AbsTol', 3.0e-3); % N=13, Nl=2
end

function [computed, reference] = compute_table1_isotropic(Nx)
R2 = 1.0;
alpha = deg2rad(30);
L = 0.75 / sin(alpha);
R1 = R2 - L*sin(alpha);
h = 0.01 * R2;
lam = laminate_isotropic(1, 0.3, 1, h);
bcs = {'SS','SS'; 'SS','C'; 'C','SS'; 'C','C'};
reference = [
    0.1493 0.6946 0.3593 0.2117 0.1907 0.2161 0.2542 0.2977 0.3516 0.4138
    0.9749 0.7451 0.4506 0.2976 0.2522 0.2640 0.2935 0.3349 0.3869 0.4485
    0.6018 0.7185 0.3763 0.2324 0.1977 0.2166 0.2525 0.2977 0.3516 0.4138
    0.9886 0.7525 0.4659 0.3153 0.2606 0.2647 0.2935 0.3349 0.3869 0.4485
];
computed = zeros(size(reference));

for q = 1:size(bcs,1)
    for n = 0:9
        params = struct('L',L, 'R1',R1, 'alpha',alpha, 'n',n, ...
                        'rho',lam.rho, 'h',lam.h, ...
                        'A',lam.A, 'B',lam.B, 'D',lam.D, 'Nx',Nx);
        sol = solve_conical_shell_modes(build_conical_shell_gdq_model(params), ...
                                        bcs{q,1}, bcs{q,2}, 3);
        computed(q,n+1) = sol.lambda(1);
    end
end
end

function [computed, reference] = compute_table2_crossply()
R2 = 1.0;
alpha = deg2rad(30);
L = 0.25 / sin(alpha);
R1 = R2 - L*sin(alpha);
h = 0.01 * R2;
Ns = [5 7 9 11 13 21];
layerCases = [Inf 2];
reference = [
    0.3158 0.2240 0.1961 0.1971 0.1976 0.1976
    0.3155 0.2226 0.1820 0.1787 0.1799 0.1799
];
computed = zeros(size(reference));

for q = 1:numel(layerCases)
    lam = laminate_shu_crossply(h, layerCases(q));
    for j = 1:numel(Ns)
        params = struct('L',L, 'R1',R1, 'alpha',alpha, 'n',0, ...
                        'rho',lam.rho, 'h',lam.h, ...
                        'A',lam.A, 'B',lam.B, 'D',lam.D, 'Nx',Ns(j));
        sol = solve_conical_shell_modes(build_conical_shell_gdq_model(params), ...
                                        'SS', 'SS', 3);
        computed(q,j) = sol.lambda(1);
    end
end
end
