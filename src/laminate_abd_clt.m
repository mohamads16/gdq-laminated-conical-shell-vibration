function lam = laminate_abd_clt(plies)
% LAMINATE_ABD_CLT  Classical Lamination Theory for ABD from ply stack.
%   lam = laminate_abd_clt(plies)

%   Input:
%       plies              - 1xN array of structs with fields:
%       E1, E2, nu12, G12  - ply orthotropic properties (plane-stress)
%       rho                - ply density
%       t                  - ply thickness
%       theta_deg          - ply fiber angle (degrees, measured from x-axis)

%   Output struct 'lam':
%     A, B, D              - 3x3 extensional, coupling, and bending stiffness matrices
%     h                    - total laminate thickness
%     rho                  - equivalent density (thickness-averaged)
%     Qbar_per_ply         - diagnostics: per-ply Qbar and z-interval
%
%   Note:
%     The current conical-shell solver follows Shu's equations using only
%     the 11, 12, 22, and 66 ABD terms. Laminates with non-negligible
%     16/26 coupling terms are rejected rather than silently truncated.


arguments
    plies (1,:) struct    % Enforce that plies is a 1xN struct array
end

np = numel(plies);              % number of plies
t  = [plies.t];                 % thickness per ply
h  = sum(t);                    % total laminate thickness
z  = linspace(-h/2, h/2, np+1); % interfaces from bottom to top

% Initialize A, B, D matrices and running mass integral
A = zeros(3); B = zeros(3); D = zeros(3);
rho_m = 0;

for k = 1:np
    p = plies(k);
    % Qb = transformed reduced stiffness matrix [Qbar] for ply k,
    % computed from its orthotropic properties and fiber angle

    [Qb, ~] = qbar_from_lamina(p.E1, p.E2, p.nu12, p.G12, deg2rad(p.theta_deg));

    % Thickness interval for ply k (from interface z(k) to z(k+1))
    dz  = z(k+1) - z(k);
    dz2 = z(k+1)^2 - z(k)^2;
    dz3 = z(k+1)^3 - z(k)^3;

    % Classical Lamination Theory thickness integrals:
    A = A + Qb * dz;
    B = B + 0.5 * Qb * dz2;
    D = D + (1/3) * Qb * dz3;

    % Integrate density through thickness for equivalent mass per area
    rho_m = rho_m + p.rho * dz;
end

% Equivalent density: average over total thickness
rho_eq = rho_m / h;

couplingTerms = [A(1,3) A(2,3) B(1,3) B(2,3) D(1,3) D(2,3) ...
                 A(3,1) A(3,2) B(3,1) B(3,2) D(3,1) D(3,2)];
couplingScale = max([1, max(abs(A),[],'all'), max(abs(B),[],'all'), max(abs(D),[],'all')]);
couplingTol = 1e-10 * couplingScale;
maxCoupling16_26 = max(abs(couplingTerms));
if maxCoupling16_26 > couplingTol
    error('laminate_abd_clt:UnsupportedCoupling16_26', ...
        ['The current conical-shell solver does not support nonzero ', ...
         'ABD 16/26 coupling terms. Maximum 16/26 term is %.3e; ', ...
         'tolerance is %.3e. Use isotropic, cross-ply, or specially ', ...
         'orthotropic laminates for this solver.'], ...
         maxCoupling16_26, couplingTol);
end

% Keep only the components used by the cone model (11,12,22,66 positions).
lam.A = [A(1,1) A(1,2) 0; A(2,1) A(2,2) 0; 0 0 A(3,3)];
lam.B = [B(1,1) B(1,2) 0; B(2,1) B(2,2) 0; 0 0 B(3,3)];
lam.D = [D(1,1) D(1,2) 0; D(2,1) D(2,2) 0; 0 0 D(3,3)];
lam.h = h;                                                 % total thickness
lam.rho = rho_eq;                                          % equivalent density
lam.maxCoupling16_26 = maxCoupling16_26;
lam.coupling16_26Tolerance = couplingTol;

% Optional diagnostics: store per-ply Qbar and z-intervals
lam.Qbar_per_ply = Qbar_per_ply(plies, z);
end

% ---------------------- helpers ----------------------
function [Qbar, Q] = qbar_from_lamina(E1,E2,nu12,G12,theta)
% plane-stress reduced stiffness and transformed stiffness
%   Outputs:
%     Qbar - 3x3 transformed reduced stiffness matrix in plate coordinates
%     Q    - 3x3 local reduced stiffness matrix in material axes

nu21 = nu12 * E2 / E1;
den  = 1 - nu12*nu21;    % denominator for plane-stress reduced stiffness

% Local reduced stiffness in material axes (plane-stress, in 1-2 coordinates)
Q11 =  E1/den;         Q22 =  E2/den;
Q12 =  nu12*E2/den;    Q66 =  G12;
Q   = [Q11 Q12 0; Q12 Q22 0; 0 0 Q66];

m = cos(theta); n = sin(theta);
m2=m*m; n2=n*n; 

% Standard CLT transformation formulas for Qbar (in 1-2 to x-y axes):
Qbar11 = Q11*m2*m2 + 2*(Q12+2*Q66)*m2*n2 + Q22*n2*n2;
Qbar22 = Q11*n2*n2 + 2*(Q12+2*Q66)*m2*n2 + Q22*m2*m2;
Qbar12 = (Q11+Q22-4*Q66)*m2*n2 + Q12*(m2*m2 + n2*n2);
Qbar16 = (Q11 - Q12 - 2*Q66)*m2*m* n - (Q22 - Q12 - 2*Q66)*m* n2*n;
Qbar26 = (Q11 - Q12 - 2*Q66)*m* n2*n - (Q22 - Q12 - 2*Q66)*m2*m* n;
Qbar66 = (Q11 + Q22 - 2*Q12 - 2*Q66)*m2*n2 + Q66*(m2*m2 + n2*n2);

% Assemble transformed reduced stiffness matrix [Qbar]
Qbar = [Qbar11 Qbar12 Qbar16;
        Qbar12 Qbar22 Qbar26;
        Qbar16 Qbar26 Qbar66];
end

function out = Qbar_per_ply(plies, z)
% Store per-ply Qbar and z-intervals (for diagnostics).
%   Output:
%     out - 1xN struct array with fields:
%         theta_deg - ply angle in degrees
%         z_bot     - bottom interface coordinate
%         z_top     - top interface coordinate
%         Qbar      - 3x3 transformed stiffness for that ply

np = numel(plies);
out(np) = struct('theta_deg',[],'z_bot',[],'z_top',[],'Qbar',[]);
zt = z;    % local copy

for k=1:np
    p = plies(k);
    % Compute Qbar for each ply (again) for reporting purposes
    [Qb,~] = qbar_from_lamina(p.E1,p.E2,p.nu12,p.G12,deg2rad(p.theta_deg));
    out(k).theta_deg = p.theta_deg;
    out(k).z_bot = zt(k); out(k).z_top = zt(k+1);
    out(k).Qbar = Qb;
end
end
