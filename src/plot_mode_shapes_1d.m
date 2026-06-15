function plot_mode_shapes_1d(mdl, sol, mode_indices)
%PLOT_MODE_SHAPES_1D  Plot U(x), V(x), W(x) for selected modes.
%
%   plot_mode_shapes_1d(mdl, sol, mode_indices)
%
%   Inputs:
%     mdl          - struct from build_conical_shell_gdq_model
%     sol          - struct from solve_conical_shell_modes
%     mode_indices - scalar or vector of mode numbers to plot, e.g. 1 or [1 2 3]
%
%   Notes:
%     - Uses the current figure; call figure(...) before this function.
%     - Modes are already normalized in solve_conical_shell_modes
%       (max |W| = 1 for each mode).

    % --- basic input handling ---
    if nargin < 3 || isempty(mode_indices)
        % Default: first up-to-three modes
        nmodes_avail = numel(sol.lambda);
        mode_indices = 1:min(3, nmodes_avail);
    end

    mode_indices = mode_indices(:).';  % row vector
    nm = numel(mode_indices);

    x = mdl.x;               % generator coordinate
    n = mdl.params.n;        % circumferential wave number

    % Set up a tiled layout: nm rows, 3 columns (U,V,W)
    % If sgtitle/tiledlayout don't exist (very old MATLAB), this will error;
    % for current versions it's fine.
    try
        tl = tiledlayout(nm, 3, 'TileSpacing','compact', 'Padding','compact');
    catch
        % Fallback to plain subplots if tiledlayout is not available
        tl = [];
    end

    for k = 1:nm
        m = mode_indices(k);

        if m < 1 || m > numel(sol.lambda)
            error('plot_mode_shapes_1d:ModeOutOfRange', ...
                  'Requested mode %d but only %d modes available.', ...
                  m, numel(sol.lambda));
        end

        U = sol.XU(:, m);
        V = sol.XV(:, m);
        W = sol.XW(:, m);

        % ----- U(x) -----
        if isempty(tl)
            subplot(nm,3,(k-1)*3+1);
        else
            nexttile((k-1)*3+1);
        end
        plot(x, U, '-o', 'LineWidth', 1.2);
        grid on;
        xlabel('x');
        ylabel('U');
        title(sprintf('Mode %d, n=%d, U(x)', m, n), 'Interpreter','none');

        % ----- V(x) -----
        if isempty(tl)
            subplot(nm,3,(k-1)*3+2);
        else
            nexttile((k-1)*3+2);
        end
        plot(x, V, '-o', 'LineWidth', 1.2);
        grid on;
        xlabel('x');
        ylabel('V');
        title(sprintf('Mode %d, n=%d, V(x)', m, n), 'Interpreter','none');

        % ----- W(x) -----
        if isempty(tl)
            subplot(nm,3,(k-1)*3+3);
        else
            nexttile((k-1)*3+3);
        end
        plot(x, W, '-o', 'LineWidth', 1.2);
        grid on;
        xlabel('x');
        ylabel('W');
        title(sprintf('Mode %d, n=%d, W(x)', m, n), 'Interpreter','none');
    end

    % Global title
    try
        sgtitle(sprintf('GDQ Conical Shell Mode Shapes (n=%d)', n), ...
                'Interpreter','none');
    catch
        % ignore if sgtitle not available
    end
end
