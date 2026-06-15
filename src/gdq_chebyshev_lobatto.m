function [x, Dx] = gdq_chebyshev_lobatto(Nx, xL, xR, maxOrder)
%GDQ_CHEBYSHEV_LOBATTO  Generalized Differential Quadrature on a 1D grid.
%   Chebyshev–Gauss–Lobatto nodes; Shu-style GDQ weights; higher orders via recurrence.
%
%   [x, Dx] = gdq_chebyshev_lobatto(Nx, xL, xR, maxOrder)
%     Nx        : number of points (>=2)
%     [xL, xR]  : interval
%     maxOrder  : highest derivative order to build (>=1)
%     x         : nodes (CGL mapped to [xL,xR])
%     Dx{k}     : k-th derivative weighting matrix (size Nx-by-Nx)

    if Nx < 2, error('Nx must be >= 2.'); end
    if ~isscalar(maxOrder) || maxOrder < 1 || maxOrder ~= floor(maxOrder)
        error('maxOrder must be a positive integer.');
    end

    % Nodes: Chebyshev–Gauss–Lobatto mapped to [xL,xR]
    x = cheb_lobatto_ascending(Nx, xL, xR);

    % Derivative matrices on x
    Dx = gdq_mats_cgl(x, maxOrder);
end

% ----------------- helpers -----------------

function x = cheb_lobatto_ascending(N, L, R)
% N Chebyshev-Lobatto nodes mapped to [L,R].
% Shu equation (37) numbers conical-shell points from x=0 to x=L.
    k  = (0:N-1).';
    s  = -cos(pi*k/(N-1));      % in [-1,1], ascending with k
    x  = (L+R)/2 + (R-L)/2*s;
end

function C = gdq_mats_cgl(x, maxOrder)
% GDQ derivative matrices using stable CGL barycentric weights.
    x = x(:); N = numel(x);
    if maxOrder > N-1
        warning('GDQ:order','maxOrder=%d exceeds N-1=%d; expect noise amplification.', ...
            maxOrder, N-1);
    end

    % First-order differentiation matrix
    w  = cgl_bary_weights(N);    % affine-invariant for CGL
    DX = x - x.';                % pairwise differences
    DX(1:N+1:end) = 1;           % protect the diagonal before division
    C1 = (w.'./w) ./ DX;         % provisional, includes bogus diagonal
    C1(1:N+1:end) = 0;           % clear diagonal exactly
    C1(1:N+1:end) = -sum(C1,2);  % row-sum closure
    C = cell(maxOrder,1); C{1} = C1;

    % Higher orders by Shu recurrence
    for k = 2:maxOrder
        Ck = zeros(N);
        for i = 1:N
            for j = 1:N
                if j ~= i
                    Ck(i,j) = k*( C1(i,j)*C{k-1}(i,i) - C{k-1}(i,j)/DX(i,j) );
                end
            end
            Ck(i,i) = -sum(Ck(i,[1:i-1,i+1:N]));
        end
        C{k} = Ck;
    end
end

function w = cgl_bary_weights(N)
% Barycentric weights for Chebyshev–Lobatto nodes (affine-invariant).
%   w_k ∝ (-1)^k * tau_k, with tau_0 = tau_{N-1} = 1/2, otherwise 1.
% Any common scaling cancels in c_ij^(1) = (w_j/w_i)/(x_i-x_j).
    k  = (0:N-1).';
    w  = (-1).^k;
    w(1)   = w(1)/2;
    w(end) = w(end)/2;
end

