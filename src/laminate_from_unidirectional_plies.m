function lam = laminate_from_unidirectional_plies(E1,E2,nu12,G12,rho,ply_t,angles_deg)
%LAMINATE_FROM_UNIDIRECTIONAL_PLIES  Build ABD data from identical UD plies.
%
%   lam = laminate_from_unidirectional_plies(E1,E2,nu12,G12,rho,ply_t,angles_deg)
%   creates a ply stack with identical material properties and one angle per
%   ply, then calls laminate_abd_clt. Angle-ply stacks that produce
%   non-negligible 16/26 ABD coupling terms are rejected because the current
%   conical-shell solver does not include those governing-equation terms.

angles_deg = angles_deg(:).';
np = numel(angles_deg);

plies = repmat(struct('E1',E1,'E2',E2,'nu12',nu12,'G12',G12, ...
                      'rho',rho,'t',ply_t,'theta_deg',0), 1, np);
for k = 1:np
    plies(k).theta_deg = angles_deg(k);
end

lam = laminate_abd_clt(plies);
end
