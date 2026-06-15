function lam = laminate_shu_crossply(h, Nlayers)
% laminate_shu_crossply
% ABD for antisymmetric cross-plied laminated conical shell
% using Shu's simplified formulas (eq. 39–40).
%
% Inputs:
%   h        shell thickness
%   Nlayers  number of layers N_l (2,4,6,10,20 or Inf for "without coupling")
%
% Output:
%   lam struct with fields: A(3x3), B(3x3), D(3x3), h, rho

    % -------- material data from eq. (39) --------
    Eu      = 1.0;               % circumferential Young's modulus (scale)
    Ex_over_Eu  = 15.0;
    Gxu_over_Eu = 0.5;
    nu_xu       = 0.25;

    Ex  = Ex_over_Eu * Eu;
    Gxu = Gxu_over_Eu * Eu;

    % Poisson in the other direction (classical reciprocity)
    nu_ux = nu_xu * Eu / Ex;

    % Reduced stiffnesses Q_ij in shell coordinates (eq. 17)
    den = 1 - nu_xu * nu_ux;
    Q11 = Ex / den;
    Q22 = Eu / den;
    Q12 = nu_xu * Eu / den;
    Q66 = Gxu;

    % -------- extensional stiffness A (eq. 40) --------
    A11 = 0.5 * h * (Q11 + Q22);
    A22 = A11;
    A12 = Q12 * h;
    A66 = Q66 * h;

    % -------- coupling stiffness B (eq. 40) --------
    if isinf(Nlayers)
        % "without coupling" case N_l = ∞  -> B = 0
        B11 = 0;   B22 = 0;
    else
        % Shu's text after eq. (40) states that the maximum-coupling
        % two-layer case has B11 = -B22 = (Q11-Q22)*h^2/8.
        % That corresponds to B11 = (h^2/(4*Nlayers))*(Q11-Q22).
        B11 = (h^2 / (4 * Nlayers)) * (Q11 - Q22);
        B22 = -B11;
    end
    B12 = 0;  B66 = 0;

    % -------- bending stiffness D (eq. 40) --------
    D11 = (h^3 / 24) * (Q11 + Q22);
    D22 = D11;
    D12 = (1/12) * Q12 * h^3;
    D66 = (1/12) * Q66 * h^3;

    % pack into 3x3 ABD, matching your cone model layout
    A = [A11 A12 0;
         A12 A22 0;
         0   0   A66];

    B = [B11 B12 0;
         B12 B22 0;
         0   0   B66];

    D = [D11 D12 0;
         D12 D22 0;
         0   0   D66];

    rho = 1.0;   % density scale (only rho*h/A11 matters in lambda)

    lam = struct('A',A, 'B',B, 'D',D, 'h',h, 'rho',rho);
end
