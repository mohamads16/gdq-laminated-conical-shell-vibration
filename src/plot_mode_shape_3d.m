function plot_mode_shape_3d(mdl, sol, mode_index, scale)
%PLOT_MODE_SHAPE_3D  3D visualization of one conical shell mode.
%
%   plot_mode_shape_3d(mdl, sol, mode_index, scale)
%
%   Inputs:
%     mdl        - struct from build_conical_shell_gdq_model
%     sol        - struct from solve_conical_shell_modes
%     mode_index - which mode to plot (1-based index)
%     scale      - optional scale factor for displacements (default: 0.1*R2)
%
%   The plot shows the deformed cone with radial displacement:
%     w(x,theta) = W(x) * cos(n * theta)
%   where W(x) is the normalized mode shape (max |W| = 1).

    if nargin < 3
        error('plot_mode_shape_3d requires at least (mdl, sol, mode_index).');
    end

    if nargin < 4 || isempty(scale)
        scale = 0.1 * mdl.geom.R2;   % 10% of large radius by default
    end

    x = mdl.x;           % generator coordinates
    R = mdl.R;           % radius at each x
    n = mdl.params.n;    % circumferential wave number
    c = mdl.geom.c;      % cos(alpha)

    nmodes = numel(sol.lambda);
    if mode_index < 1 || mode_index > nmodes
        error('plot_mode_shape_3d:ModeOutOfRange', ...
              'Requested mode %d but only %d modes available.', ...
              mode_index, nmodes);
    end

    % W-shape for this mode (already normalized so max|W| = 1)
    W = sol.XW(:, mode_index);

    % Theta grid around circumference
    nTheta = 120;
    theta  = linspace(0, 2*pi, nTheta);

    % 2D grids: theta vs x
    [Xgrid, Tgrid] = meshgrid(x, theta);   % Xgrid: x, Tgrid: theta

    % Radius grid R(x)
    Rrow  = R(:).';                        % 1 x N
    Rgrid = repmat(Rrow, nTheta, 1);       % nTheta x N

    % W(x) grid, then modal pattern in theta: W(x)*cos(n*theta)
    Wrow   = W(:).';
    Wgrid  = repmat(Wrow, nTheta, 1);
    Wtheta = Wgrid .* cos(n * Tgrid);

    % Undeformed vs deformed radius
    r0    = Rgrid;
    r_def = r0 + scale * Wtheta;           % radial deformation

    % Axial coordinate: projection of generator along cone axis: z = x*cos(alpha)
    Zgrid = Xgrid * c;

    % Convert to Cartesian coordinates
    Xcart = r_def .* cos(Tgrid);
    Ycart = r_def .* sin(Tgrid);

    % --- plotting ---
    figName = sprintf('3D mode shape, n=%d, mode=%d', n, mode_index);
    figure('Name', figName);

    surf(Xcart, Ycart, Zgrid, Wtheta, 'EdgeColor','none');
    axis equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('Mode %d (n=%d), lambda = %.4f', ...
          mode_index, n, sol.lambda(mode_index)), ...
          'Interpreter','none');
    colorbar;
    view(135, 30);
    grid on;
end
