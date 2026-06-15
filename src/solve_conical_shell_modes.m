function sol = solve_conical_shell_modes(mdl, bcSmall, bcLarge, nmodes)

% Inputs:
%   mdl     - model struct from build_conical_shell_gdq_model, containing:
%               x, D1..D4, S, BC, params, geom
%   bcSmall - boundary type at small radius: 'SS' or 'C'
%   bcLarge - boundary type at large radius: 'SS' or 'C'
%   nmodes  - (optional) number of modes to compute (default 12)

% Output:
%     ''sol'' struct with fields:
%     lambda                     - non-dimensional frequencies
%     omega                      - dimensional frequencies
%     mu                         - eigenvalues mu = rho*h*omega^2
%     XU,XV,XW                   - mode shapes for U,V,W
%     idx.I, idx.B               - interior and boundary DOF indices
%     Kred, T                    - reduced operator and boundary map
%     meta.bcSmall, meta.bcLarge - boundary types used

% If nmodes not given or empty, default to 12 modes
if nargin < 4 || isempty(nmodes), nmodes = 12; end

% Shorthands
N  = numel(mdl.x);
D1 = mdl.D1; D2 = mdl.D2; D3 = mdl.D3; D4 = mdl.D4;
S  = mdl.S;

% ---------------- Blocks (vectorized: diag(S)*Dk) ---------------------------------
sd = S.spdiag;

KU = sd(S.S110) + sd(S.S111)*D1 + sd(S.S112)*D2;
KV = sd(S.S120) + sd(S.S121)*D1;
KW = sd(S.S130) + sd(S.S131)*D1 + sd(S.S132)*D2 + sd(S.S133)*D3;

LU = sd(S.S210) + sd(S.S211)*D1;
LV = sd(S.S220) + sd(S.S221)*D1 + sd(S.S222)*D2;
LW = sd(S.S230) + sd(S.S231)*D1 + sd(S.S232)*D2;

MU = sd(S.S310) + sd(S.S311)*D1 + sd(S.S312)*D2 + sd(S.S313)*D3;
MV = sd(S.S320) + sd(S.S321)*D1 + sd(S.S322)*D2;
MW = sd(S.S330) + sd(S.S331)*D1 + sd(S.S332)*D2 + sd(S.S333)*D3 + sd(S.S334)*D4;


% Assemble full 3N x 3N matrix A for the coupled system:
% Rows are (U-equation; V-equation; W-equation).
% Columns correspond to stacked unknowns [U; V; W].
A = sparse([KU, KV, KW;
            LU, LV, LW;
            MU, MV, MW]);   % LHS operator in (29)

%% ---------------- Interior rows (GDQ locations) -----------------------------------
% For U and V, use interior nodes i = 2..N-1
% For W, use interior nodes i = 3..N-2 (leaving two nodes near each edge for BCs)
rowU_I = (2:N-1);                     % U rows to keep as interior equations
rowV_I = N + (2:N-1);                 % V rows in A: offset by N
rowW_I = 2*N + (3:N-2);               % W rows in A: offset by 2N
rows_I = [rowU_I, rowV_I, rowW_I];    % all interior equation rows

% Column index helpers for the stacked unknown vector X = [U;V;W]
colU = 1:N; colV = N+(1:N); colW = 2*N+(1:N);

% Interior unknown indices: U(2..N-1), V(2..N-1), W(3..N-2)
idxU_I = colU(2:N-1); idxV_I = colV(2:N-1); idxW_I = colW(3:N-2);
idxI   = [idxU_I, idxV_I, idxW_I];    % overall interior DOF indices

% Boundary unknown indices:
%   U(1), U(N),
%   V(1), V(N),
%   W(1), W(2), W(N-1), W(N)
idxB = [colU(1), colU(N), colV(1), colV(N), colW(1), colW(2), colW(N-1), colW(N)];

%% ---------------- Boundary equations B * X = 0  -----------------------------------
% Build 8 boundary equations (4 at small edge, 4 at large edge)
% in terms of the full stacked unknowns X = [U;V;W].
B = build_boundary_rows(mdl, bcSmall, bcLarge);  % size 8 x 3N

% Partition A into interior blocks corresponding to rows_I
A_I  = A(rows_I, :);
A_II = A_I(:, idxI);
A_IB = A_I(:, idxB);

% Likewise partition B into interior- and boundary-unknown parts
B_I  = B(:, idxI);
B_B  = B(:, idxB);

boundaryRcond = rcond(full(B_B));
boundaryRcondWarnThreshold = 1e-12;
if boundaryRcond < boundaryRcondWarnThreshold
    warning('solve_conical_shell_modes:IllConditionedBoundaryBlock', ...
        'Boundary block rcond is %.3e; mode results may be unreliable.', ...
        boundaryRcond);
end

% From B * X = 0:
%   B_I * X_I + B_B * X_B = 0  =>  X_B = -B_B^{-1} * B_I * X_I
% Boundary map X_B = T * X_I
T    = -(B_B \ B_I);

% *** SIGN FIX ***
% The continuous equation (29) is written as LHS = -rho*h*omega^2 * X.
% We built +LHS here, so to match eigenproblem form (34)-(36) we negate.
% With condensation, we have:
%   (A_II + A_IB*T) X_I = -rho*h*omega^2 X_I
% So define Kred so that:
%   Kred * X_I = rho*h*omega^2 * X_I

%Kred: reduced stiffness-like matrix for the interior unknowns
Kred = -(A_II + A_IB*T);      % == -(A_II - A_IB*(B_B\B_I))


% ---------------- Solve ------------------------------------------------------------

% Number of modes to compute: cannot exceed (dimension of Kred - 2) for safety
requestedModes = min(nmodes, size(Kred,1)-2);
if requestedModes < 1
    error('solve_conical_shell_modes:TooFewDofs', ...
        'Reduced eigenproblem is too small to compute modes.');
end

% Options for eigs (iterative eigen-solver)
opts.tol = 1e-10;
opts.maxit = 2e4;
opts.isreal = true;            % <-- keep; drop opts.issym

rho   = mdl.params.rho;  h = mdl.params.h;
A11   = mdl.params.A(1,1);  R2 = mdl.geom.R2;

eigMeta = struct();
candidateModes = min(size(Kred,1)-2, max(requestedModes + 8, 2*requestedModes));
if candidateModes < requestedModes
    candidateModes = requestedModes;
end

try
    [Vraw, Draw] = eigs(Kred, candidateModes, 0, opts);
    eigMeta.solver = 'eigs';
catch
    [Vraw, Draw] = eig(full(Kred));
    eigMeta.solver = 'eig';
end

rawEigenvalues = diag(Draw);
[keep, filterMeta] = filter_physical_eigenvalues(rawEigenvalues);

if nnz(keep) < requestedModes && ~strcmp(eigMeta.solver, 'eig')
    [Vraw, Draw] = eig(full(Kred));
    rawEigenvalues = diag(Draw);
    [keep, filterMeta] = filter_physical_eigenvalues(rawEigenvalues);
    eigMeta.solver = 'eig';
    eigMeta.fallbackReason = 'Too few physical eigenvalues from eigs candidate set.';
end

acceptedIdx = find(keep);
if isempty(acceptedIdx)
    error('solve_conical_shell_modes:NoPhysicalEigenvalues', ...
        'No positive real eigenvalues passed the physical-mode filter.');
end

lambdaCandidates = R2 * sqrt(real(rawEigenvalues(acceptedIdx)) / A11);
[~, acceptedOrder] = sort(lambdaCandidates(:));
selectedIdx = acceptedIdx(acceptedOrder(1:min(requestedModes, numel(acceptedIdx))));

Vred = Vraw(:, selectedIdx);
rawSelectedEigenvalues = rawEigenvalues(selectedIdx);
mu = real(rawSelectedEigenvalues);         % = rho*h*omega^2
k = numel(mu);

% Dimensional frequency: omega = sqrt(mu / (rho*h))
omega  = sqrt(mu / (rho*h));

% Non-dimensional frequency parameter:
%   lambda = R2 * sqrt(mu / A11) = R2 * sqrt(rho*h/A11) * omega
lambda = R2 * sqrt(mu / A11);

%% ---------------- Reconstruct full modes at all nodes ------------------------------
% X_I -> X_B -> X, then split into U, V, W
XU = zeros(N, k); XV = zeros(N, k); XW = zeros(N, k);

for j = 1:k
    XI = Vred(:,j);         % interior DOFs for mode j
    XB = T * XI;            % boundary DOFs from condensation

    X = zeros(3*N,1);       % full DOF vector [U;V;W]
    X(idxI) = XI;           % fill interior entries
    X(idxB) = XB;           % fill boundary entries

    % Split X back into U, V, W at all nodes
    XU(:,j) = X(colU);
    XV(:,j) = X(colV);
    XW(:,j) = X(colW);
end

% Normalize (peak W = 1 for readability)
for j = 1:k
    sc = max(abs(XW(:,j))); if sc>0, XU(:,j)=XU(:,j)/sc; XV(:,j)=XV(:,j)/sc; XW(:,j)=XW(:,j)/sc; end
end

% ---------------- Package ----------------------------------------------------------
[lambda, ord] = sort(lambda(:));
omega = omega(ord); mu = mu(ord);
XU = XU(:,ord); XV = XV(:,ord); XW = XW(:,ord);
rawSelectedEigenvalues = rawSelectedEigenvalues(ord);

sol.lambda = lambda; sol.omega = omega; sol.mu = mu;
sol.XU = XU; sol.XV = XV; sol.XW = XW;
sol.idx.I = idxI; sol.idx.B = idxB;
sol.Kred  = Kred; sol.T = T;
sol.meta.bcSmall = bcSmall; sol.meta.bcLarge = bcLarge;
sol.meta.boundaryBlockRcond = boundaryRcond;
sol.meta.boundaryBlockRcondWarningThreshold = boundaryRcondWarnThreshold;
sol.meta.eigenFilter = filterMeta;
sol.meta.eigenSolver = eigMeta;
sol.meta.rawSelectedEigenvalues = rawSelectedEigenvalues;
end

%===============================================================================
function [keep, meta] = filter_physical_eigenvalues(eigenvalues)
% Keep eigenvalues that represent positive physical values of rho*h*omega^2.

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

%===============================================================================
function B = build_boundary_rows(mdl, bcSmall, bcLarge)
% Build 8 boundary equations (4 per edge) in the stacked unknowns [U;V;W].

N  = numel(mdl.x);
D1 = mdl.D1; D2 = mdl.D2;
colU = 1:N; colV = N+(1:N); colW = 2*N+(1:N);

B = sparse(8, 3*N); r = 0;

% small-edge row index in GDQ matrices
iS = 1;
% large-edge row index
iL = N;

% --- small end ---
switch upper(strtrim(bcSmall))
    case 'SS'
        [B, r] = add_SS_rows(B, r, iS, mdl.BC.small, colU, colV, colW, D1, D2, N);
    case 'C'
        [B, r] = add_C_rows(B, r, iS, colU, colV, colW, D1);
    otherwise
        error('bcSmall must be ''SS'' or ''C''.');
end

% --- large end ---
switch upper(strtrim(bcLarge))
    case 'SS'
        [B, ~] = add_SS_rows(B, r, iL, mdl.BC.large, colU, colV, colW, D1, D2, N);
    case 'C'
        [B, ~] = add_C_rows(B, r, iL, colU, colV, colW, D1);
    otherwise
        error('bcLarge must be ''SS'' or ''C''.');
end
end

%===============================================================================
function [B, r] = add_SS_rows(B, r, i, BCedge, colU, colV, colW, D1, D2, N)
% Simply-supported: V=0, W=0, plus mixed U/W conditions with {a_i},{b_i}.
a = BCedge.a;  % [a1..a5]
b = BCedge.b;  % [b1..b5]

% 1) V(i)=0
r = r+1; B(r, colV(i)) = 1;

% 2) W(i)=0
r = r+1; B(r, colW(i)) = 1;

% 3) U_x + a1 U + a2 V + a3 W + a4 W_x + a5 W_xx = 0
r = r+1;
Bu = sparse(1, 3*N);
Bu(1, colU) = D1(i,:);          % U_x(i)
Bu(1, colU(i)) = Bu(1, colU(i)) + a(1);                    % + a1 U(i)
Bu(1, colV(i)) = Bu(1, colV(i)) + a(2);                    % + a2 V(i)
Bu(1, colW(i)) = Bu(1, colW(i)) + a(3);                    % + a3 W(i)
Bu(1, colW)    = Bu(1, colW)    + a(4)*D1(i,:) + a(5)*D2(i,:); % + a4 W_x + a5 W_xx
B(r,:) = Bu;

% 4) W_xx + b1 U + b2 U_x + b3 V + b4 W + b5 W_x = 0
r = r+1;
Bw = sparse(1, 3*N);
Bw(1, colW)    = D2(i,:);       % W_xx(i)
Bw(1, colU(i)) = Bw(1, colU(i)) + b(1);                    % + b1 U(i)
Bw(1, colU)    = Bw(1, colU)    + b(2)*D1(i,:);            % + b2 U_x(i)
Bw(1, colV(i)) = Bw(1, colV(i)) + b(3);                    % + b3 V(i)
Bw(1, colW(i)) = Bw(1, colW(i)) + b(4);                    % + b4 W(i)
Bw(1, colW)    = Bw(1, colW)    + b(5)*D1(i,:);            % + b5 W_x(i)
B(r,:) = Bw;
end

%===============================================================================
function [B, r] = add_C_rows(B, r, i, colU, colV, colW, D1)
% Clamped: U=0, V=0, W=0, W_x=0
% 1) U(i)=0
r = r+1; B(r, colU(i)) = 1;
% 2) V(i)=0
r = r+1; B(r, colV(i)) = 1;
% 3) W(i)=0
r = r+1; B(r, colW(i)) = 1;
% 4) W_x(i)=0
r = r+1; B(r, colW) = D1(i,:);
end
