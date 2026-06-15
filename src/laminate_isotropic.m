function lam = laminate_isotropic(E,nu,rho,h)
%LAMINATE_ISOTROPIC  Equivalent ABD data for an isotropic shell.

G = E/(2*(1+nu));    % Shear modulus G for isotropic material

% Plane-stress reduced stiffness matrix [Q] for isotropic lamina:
Q = [E/(1-nu^2), nu*E/(1-nu^2), 0;
     nu*E/(1-nu^2), E/(1-nu^2), 0;
     0, 0, G];

% Classical Lamination Theory: A = Q*h, B = 0, D = Q*h^3/12 for a single layer
A = Q*h; 
B = zeros(3); 
D = Q*h^3/12;

% Build lam struct:
lam = struct('A',[A(1,1) A(1,2) 0; A(2,1) A(2,2) 0; 0 0 A(3,3)], ...
             'B',B, ...
             'D',[D(1,1) D(1,2) 0; D(2,1) D(2,2) 0; 0 0 D(3,3)], ...
             'h',h, 'rho',rho);
end
