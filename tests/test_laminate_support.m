function tests = test_laminate_support
tests = functiontests(localfunctions);
end

function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));
end

function testCrossPlyLaminateIsAccepted(testCase)
lam = laminate_from_unidirectional_plies(15, 1, 0.25, 0.5, 1, 0.005, [0 90 90 0]);

verifySize(testCase, lam.A, [3 3]);
verifyLessThanOrEqual(testCase, lam.maxCoupling16_26, lam.coupling16_26Tolerance);
end

function testAnglePlyLaminateIsRejected(testCase)
verifyError(testCase, @() laminate_from_unidirectional_plies(15, 1, 0.25, 0.5, 1, 0.01, 45), ...
    'laminate_abd_clt:UnsupportedCoupling16_26');
end
