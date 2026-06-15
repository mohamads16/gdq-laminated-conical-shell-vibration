function tests = test_gdq_derivatives
tests = functiontests(localfunctions);
end

function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));
end

function testPolynomialDerivatives(testCase)
[x, D] = gdq_chebyshev_lobatto(13, 0, 1, 4);

for p = 0:4
    f = x.^p;

    if p == 0
        exact1 = zeros(size(x));
    else
        exact1 = p*x.^(p-1);
    end
    verifyLessThanOrEqual(testCase, max(abs(D{1}*f - exact1)), 1e-10);

    if p < 2
        exact2 = zeros(size(x));
    else
        exact2 = p*(p-1)*x.^(p-2);
    end
    verifyLessThanOrEqual(testCase, max(abs(D{2}*f - exact2)), 1e-9);
end
end

function testNodeOrdering(testCase)
[x, ~] = gdq_chebyshev_lobatto(7, 0, 2, 4);
verifyEqual(testCase, x(1), 0, 'AbsTol', 1e-14);
verifyEqual(testCase, x(end), 2, 'AbsTol', 1e-14);
verifyTrue(testCase, all(diff(x) > 0));
end
