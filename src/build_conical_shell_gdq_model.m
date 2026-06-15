function mdl = build_conical_shell_gdq_model(params)

arguments
    params struct
end

%% --- Parameters
L     = params.L;    
R1    = params.R1;   
alpha = params.alpha;
n     = params.n;    
A     = params.A;    
B     = params.B;   
D     = params.D;
Nx    = params.Nx;

assert(isscalar(L) && L>0 && Nx>=5, 'Use L>0 and Nx>=5 (W needs up to 4th order).');

% --- GDQ nodes and derivative matrices on [0,L]
[x, Dx] = gdq_chebyshev_lobatto(Nx, 0, L, 4);
D1 = Dx{1}; D2 = Dx{2}; D3 = Dx{3}; D4 = Dx{4};

% --- geometry at nodes
s = sin(alpha); c = cos(alpha);
R = R1 + x*s;
invR  = 1./R; invR2 = invR.^2; invR3 = invR.^3; invR4 = invR.^4;
R2_big = R1 + L*sin(alpha);   % exact large-edge radius


% --- ABD entries (consistent with Shu's setup)
A11=A(1,1); A12=A(1,2); A22=A(2,2); A66=A(3,3);
B11=B(1,1); B12=B(1,2); B22=B(2,2); B66=B(3,3);
D11=D(1,1); D12=D(1,2); D22=D(2,2); D66=D(3,3);

%% ======================= Appendix 1: S_{ijk}(x) at each node =======================
% U-equation (25a)
S110 = -(A66*n^2).*invR2 - A22*(s.^2).*invR2;
S111 =  A11*s.*invR;
S112 =  A11*ones(Nx,1);

S120 = -(B12+B22+2*B66).*n.*s.*c.*invR3 - (A22 + A66).*n.*s.*invR2;
S121 =  (B12+2*B66).*n.*c.*invR2 + (A12 + A66).*n.*invR;

S130 = -(B12+B22+2*B66).*n.^2.*s.*invR3 - A22*s.*c.*invR2;
S131 =  (B12+2*B66).*n.^2.*invR2 + B22*(s.^2).*invR2 + A12*c.*invR;
S132 = -B11*s.*invR;
S133 = -B11*ones(Nx,1);

% V-equation (25b)
S210 = -(B22 - B66).*n.*s.*c.*invR3 - (A22 + A66).*n.*s.*invR2;
S211 = -(B12*n*c.*invR2 + B66*n*c.*invR2 + A12*n.*invR + A66*n.*invR);
S220 = -((D22*n.^2 - 4*D66*(s.^2)).*(c.^2)).*invR4 ...
       - ( (2*B22*n.^2 - B66*(s.^2)).*c ).*invR3 ...
       - (A22*n.^2 + A66*(s.^2)).*invR2;
S221 = -4*D66*s.*(c.^2).*invR3 - B66*s.*c.*invR2 + A66*s.*invR;
S222 =  A66 + 2*D66*(c.^2).*invR2 + 3*B66*c.*invR;

S230 = -((D22*n.^2 - 4*D66*(s.^2)).*n.*c).*invR4 ...
       - (B22*n.*(c.^2 + n.^2)).*invR3 - (A22*n.*c).*invR2;
S231 =  (D22 - 4*D66).*n.*s.*c.*invR3 + B22*n.*s.*invR2;
S232 =  (D12 + 2*D66).*n.*c.*invR2 + (B12 + 2*B66).*n.*invR;

% W-equation (25c)
S310 = ((-B22*n.^2 + 2*B66*n.^2 + B22*(s.^2)).*s).*invR3 - A22*s.*c.*invR2;
S311 = -(B12*n.^2 + 2*B66*n.^2 + B22*(s.^2)).*invR2 - A12*c.*invR;
S312 =  2*B11*s.*invR;
S313 =  B11*ones(Nx,1);

S320 = ((-D22*n.^2 + 2*D12*(s.^2) + 2*D22*(s.^2) + 8*D66*(s.^2)).*n.*c).*invR4 ...
       - A22*n.*c.*invR2 ...
       + (-B22*n.^2 - B22*(c.^2) + B22*(s.^2) + 2*B66*(s.^2)).*n.*invR3;
S321 = -(D22 + 2*D12 + 8*D66).*n.*c.*s.*invR3 - (B22 + 2*B66).*n.*s.*invR2;
S322 =  (D12 + 4*D66).*n.*c.*invR2 + (B12 + 2*B66).*n.*invR;

S330 = ((-D22*n.^2 + 2*D12*(s.^2) + 2*D22*(s.^2) + 8*D66*(s.^2)).*n.^2).*invR4 ...
       + ((-2*n.^2 + s.^2).*B22.*c).*invR3 - A22*(c.^2).*invR2;
S331 = -(D22*(s.^2) + 2*D12*n.^2 + 8*D66*n.^2).*s.*invR3;
S332 =  (2*D12*n.^2 + 4*D66*n.^2 + D22*(s.^2)).*invR2 + 2*B12*c.*invR;
S333 = -2*D11*s.*invR;
S334 = -D11*ones(Nx,1);


% Stash as vectors + sparse-diagonal helper for solve_conical_shell_modes.
S = struct( ...
    'S110',S110, 'S111',S111, 'S112',S112, ...
    'S120',S120, 'S121',S121, ...
    'S130',S130, 'S131',S131, 'S132',S132, 'S133',S133, ...
    'S210',S210, 'S211',S211, ...
    'S220',S220, 'S221',S221, 'S222',S222, ...
    'S230',S230, 'S231',S231, 'S232',S232, ...
    'S310',S310, 'S311',S311, 'S312',S312, 'S313',S313, ...
    'S320',S320, 'S321',S321, 'S322',S322, ...
    'S330',S330, 'S331',S331, 'S332',S332, 'S333',S333, 'S334',S334);
S.spdiag = @(v) spdiags(v,0,Nx,Nx);

%% ======================= Appendix 2: boundary scalars at x1 and xN =================
a1 = @(Rval)  A12*s ./ (A11*Rval);
a2 = @(Rval) (B12*n*c) ./ (A11*Rval.^2) + (A12*n) ./ (A11*Rval);
a3 = @(Rval)  A12*c ./ (A11*Rval) + (B12*n.^2) ./ (A11*Rval.^2);
a4 = @(Rval) -B12*s ./ (A11*Rval);
a5 =        - B11 / A11;

b1 = @(Rval) -B12*s ./ (D11*Rval);
b2 =        - B11 / D11;
b3 = @(Rval) -D12*n*c ./ (D11*Rval.^2) - (B12*n) ./ (D11*Rval);
b4 = @(Rval) -B12*c ./ (D11*Rval) - (D12*n.^2) ./ (D11*Rval.^2);
b5 = @(Rval)  D12*s ./ (D11*Rval);

BC.small.a = [a1(R(1)) a2(R(1)) a3(R(1)) a4(R(1)) a5];
BC.small.b = [b1(R(1)) b2           b3(R(1)) b4(R(1)) b5(R(1))];
BC.large.a = [a1(R(end)) a2(R(end)) a3(R(end)) a4(R(end)) a5];
BC.large.b = [b1(R(end)) b2          b3(R(end)) b4(R(end)) b5(R(end))];
% These feed directly into the GDQ boundary rows from Shu equations (30)-(31).

% ----------------------- package model data ----------------------------------------
% Nodal positions and radius
mdl.x = x; 
mdl.R = R; 

% Geometry helper data
mdl.geom.s = s; mdl.geom.c = c; mdl.geom.R2 = R2_big;

% GDQ derivative matrices
mdl.D1 = D1; mdl.D2 = D2; mdl.D3 = D3; mdl.D4 = D4;

% S_ijk coefficient struct
mdl.S  = S;

% Boundary-condition coefficients
mdl.BC = BC;

% Copy of original params for reference
mdl.params = params;
end
