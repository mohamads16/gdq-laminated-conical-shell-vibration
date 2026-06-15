function [keep, meta] = filter_physical_eigenvalues(eigenvalues)
%FILTER_PHYSICAL_EIGENVALUES  Identify positive real shell eigenvalues.
%
%   [keep, meta] = filter_physical_eigenvalues(eigenvalues) returns a
%   logical mask for eigenvalues that can represent mu = rho*h*omega^2.
%   Values with non-negligible imaginary parts or nonpositive real parts
%   are rejected. Tolerances and rejected values are returned in meta.

eigenvalues = eigenvalues(:);
imagAbsTol = 1e-10;
imagRelTol = 1e-8;
imagLimit = max(imagAbsTol, imagRelTol * max(1, abs(eigenvalues)));
imagReject = abs(imag(eigenvalues)) > imagLimit;

positiveTol = 1e-12 * max(1, max(abs(real(eigenvalues))));
nonpositiveReject = real(eigenvalues) <= positiveTol;

keep = ~(imagReject | nonpositiveReject);
rejectedIdx = find(~keep);
reasons = strings(numel(rejectedIdx), 1);
for j = 1:numel(rejectedIdx)
    idx = rejectedIdx(j);
    parts = strings(0,1);
    if imagReject(idx)
        parts(end+1,1) = "non-negligible imaginary part"; %#ok<AGROW>
    end
    if nonpositiveReject(idx)
        parts(end+1,1) = "nonpositive physical eigenvalue"; %#ok<AGROW>
    end
    reasons(j) = strjoin(parts, "; ");
end

meta = struct();
meta.imagAbsTol = imagAbsTol;
meta.imagRelTol = imagRelTol;
meta.positiveTol = positiveTol;
meta.rawEigenvalues = eigenvalues;
meta.acceptedIndices = find(keep);
meta.rejectedIndices = rejectedIdx;
meta.rejectedEigenvalues = eigenvalues(rejectedIdx);
meta.rejectedReasons = reasons;
end
