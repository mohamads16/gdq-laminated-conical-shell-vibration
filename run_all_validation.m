%RUN_ALL_VALIDATION  Run examples and tests used for repository validation.

clear; clc;

run('examples/validate_table1_isotropic.m');
run('examples/validate_table2_crossply.m');
run('examples/convergence_isotropic_ssss.m');

oldFigureVisibility = get(groot, 'DefaultFigureVisible');
cleanupObj = onCleanup(@() set(groot, 'DefaultFigureVisible', oldFigureVisibility));
set(groot, 'DefaultFigureVisible', 'off');
run('examples/plot_isotropic_mode_shapes.m');
close all;

results = runtests('tests');
disp(table(results));
assertSuccess(results);
