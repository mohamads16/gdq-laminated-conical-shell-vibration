function tests = test_eigenvalue_filter
tests = functiontests(localfunctions);
end

function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));
end

function testPositiveRealEigenvaluesAccepted(testCase)
eigsIn = [1; 2.5; 10];
[keep, meta] = filter_physical_eigenvalues(eigsIn);

verifyEqual(testCase, keep, true(size(eigsIn)));
verifyEmpty(testCase, meta.rejectedEigenvalues);
end

function testNegativeZeroAndComplexRejected(testCase)
eigsIn = [-1; 0; 2 + 1e-3i; 4];
[keep, meta] = filter_physical_eigenvalues(eigsIn);

verifyEqual(testCase, keep, [false; false; false; true]);
verifyEqual(testCase, meta.rejectedIndices, [1; 2; 3]);
verifyEqual(testCase, meta.rejectedEigenvalues, eigsIn(1:3));
end

function testTinyImaginaryRoundoffAccepted(testCase)
eigsIn = [1 + 1e-12i; 5 - 1e-12i];
[keep, meta] = filter_physical_eigenvalues(eigsIn);

verifyEqual(testCase, keep, [true; true]);
verifyEmpty(testCase, meta.rejectedEigenvalues);
end
