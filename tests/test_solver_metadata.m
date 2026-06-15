function tests = test_solver_metadata
tests = functiontests(localfunctions);
end

function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));
end

function testEigenFilterAndBoundaryDiagnosticsExist(testCase)
R2 = 1.0;
alpha = deg2rad(30);
L = 0.75 / sin(alpha);
R1 = R2 - L*sin(alpha);
h = 0.01 * R2;

lam = laminate_isotropic(1, 0.3, 1, h);
params = struct('L',L, 'R1',R1, 'alpha',alpha, 'n',3, ...
                'rho',lam.rho, 'h',lam.h, ...
                'A',lam.A, 'B',lam.B, 'D',lam.D, 'Nx',33);

sol = solve_conical_shell_modes(build_conical_shell_gdq_model(params), 'SS', 'SS', 3);

verifyTrue(testCase, isfield(sol.meta, 'boundaryBlockRcond'));
verifyTrue(testCase, isfield(sol.meta, 'boundaryBlockRcondWarningThreshold'));
verifyGreaterThan(testCase, sol.meta.boundaryBlockRcond, 0);

verifyTrue(testCase, isfield(sol.meta, 'eigenFilter'));
verifyTrue(testCase, isfield(sol.meta.eigenFilter, 'imagAbsTol'));
verifyTrue(testCase, isfield(sol.meta.eigenFilter, 'imagRelTol'));
verifyTrue(testCase, isfield(sol.meta.eigenFilter, 'positiveTol'));
verifyTrue(testCase, isfield(sol.meta.eigenFilter, 'rejectedEigenvalues'));
verifyTrue(testCase, isfield(sol.meta.eigenFilter, 'rejectedReasons'));
verifyTrue(testCase, all(sol.mu > 0));
verifyTrue(testCase, all(abs(imag(sol.meta.rawSelectedEigenvalues)) <= ...
    max(sol.meta.eigenFilter.imagAbsTol, ...
        sol.meta.eigenFilter.imagRelTol * max(1, abs(sol.meta.rawSelectedEigenvalues)))));
end

