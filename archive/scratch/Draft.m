clear;clc;
N=13;

rowU_I = (2:N-1);
disp(rowU_I)

lambda_num = zeros(10,1);
for n = 0:9
    params = struct('L',geom.L,'R1',geom.R1,'alpha',geom.alpha, ...
                    'n',n,'rho',lam.rho,'h',lam.h, ...
                    'A',lam.A,'B',lam.B,'D',lam.D,'Nx',Nx);
    mdl = step1(params);
    sol = step2(mdl,'SS','SS',1);   % only lowest mode
    lambda_num(n+1) = sol.lambda(1);
end