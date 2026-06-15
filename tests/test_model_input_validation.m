function tests = test_model_input_validation
tests = functiontests(localfunctions);
end

function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));
end

function testMissingRequiredField(testCase)
params = valid_params();
params = rmfield(params, 'A');

verifyError(testCase, @() build_conical_shell_gdq_model(params), ...
    'build_conical_shell_gdq_model:MissingField');
end

function testInvalidNx(testCase)
params = valid_params();
params.Nx = 4;
verifyError(testCase, @() build_conical_shell_gdq_model(params), ...
    'build_conical_shell_gdq_model:InvalidNx');

params = valid_params();
params.Nx = 9.5;
verifyError(testCase, @() build_conical_shell_gdq_model(params), ...
    'build_conical_shell_gdq_model:InvalidNx');
end

function testInvalidPositiveScalars(testCase)
for fieldName = ["L","rho","h","R1"]
    params = valid_params();
    params.(fieldName) = 0;
    verifyError(testCase, @() build_conical_shell_gdq_model(params), ...
        'build_conical_shell_gdq_model:InvalidPositiveScalar');
end
end

function testInvalidWaveNumber(testCase)
params = valid_params();
params.n = -1;
verifyError(testCase, @() build_conical_shell_gdq_model(params), ...
    'build_conical_shell_gdq_model:InvalidWaveNumber');

params = valid_params();
params.n = 1.5;
verifyError(testCase, @() build_conical_shell_gdq_model(params), ...
    'build_conical_shell_gdq_model:InvalidWaveNumber');
end

function testInvalidMatricesAndA11(testCase)
params = valid_params();
params.A = eye(2);
verifyError(testCase, @() build_conical_shell_gdq_model(params), ...
    'build_conical_shell_gdq_model:InvalidMatrix');

params = valid_params();
params.B(1,1) = NaN;
verifyError(testCase, @() build_conical_shell_gdq_model(params), ...
    'build_conical_shell_gdq_model:InvalidMatrix');

params = valid_params();
params.A(1,1) = 0;
verifyError(testCase, @() build_conical_shell_gdq_model(params), ...
    'build_conical_shell_gdq_model:InvalidA11');
end

function testInvalidRadiusOverGrid(testCase)
params = valid_params();
params.alpha = -pi/2;
verifyError(testCase, @() build_conical_shell_gdq_model(params), ...
    'build_conical_shell_gdq_model:InvalidRadius');
end

function params = valid_params()
R2 = 1.0;
alpha = deg2rad(30);
L = 0.75 / sin(alpha);
R1 = R2 - L*sin(alpha);
h = 0.01 * R2;
lam = laminate_isotropic(1, 0.3, 1, h);

params = struct('L',L, 'R1',R1, 'alpha',alpha, 'n',3, ...
                'rho',lam.rho, 'h',lam.h, ...
                'A',lam.A, 'B',lam.B, 'D',lam.D, 'Nx',13);
end

