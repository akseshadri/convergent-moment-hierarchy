%% run_all_figures.m — Master script: reproduce every figure
%
% Companion code for:
%   Seshadri, A. K. and Pal Majumder, A. (2026). A convergent moment
%   hierarchy for the variance of weighted spatial-mean estimators under
%   missing-at-random sampling.
%
% Runs the three section scripts in order. Each opens its own MATLAB figure windows
% Nothing is written to disk. Set the working directory to src/ 
% before running (the scripts load the data from ../data/).
%
%   PFigs_Section2 — Fig 1 (one parameter rho governs convergence and
%                    accuracy), Fig 2 (good-event / rare-event split).
%   PFigs_Section3 — Fig 3 (Monte Carlo vs the four truncations), Fig 4
%                    (mixed moments), Fig 5 (source decomposition, amplified
%                    means), Fig 6 (stress regimes).
%   PFigs_SI       — Fig S1, S2 (source decomposition: raw and reduced
%                    variance).
%
% Shared functions, used consistently by every script:
%   truncations_scenarioI    uniform-weight (Scenario I) closed forms
%   truncations_general      arbitrary-weight closed forms and rho
%   ratio_var_exact_uniform  exact Var(rhat) for uniform weights (no Monte Carlo)
%   mc_truevar               Monte-Carlo moments of rhat under MAR
%   apply_fig_style          one consistent figure style
%   panel_label              bold (a)/(b)/(c) panel letters
%
% Requirements: Statistics and Machine Learning Toolbox (binornd, binopdf);
% MATLAB R2016b or later (yyaxis, histogram). Parallel Computing Toolbox is
% optional (mc_truevar uses parfor, which runs serially without a pool).

clear; close all

PFigs_Section2      % Fig 1, Fig 2
PFigs_Section3      % Fig 3, Fig 4, Fig 5, Fig 6
PFigs_SI            % Fig S1, Fig S2

fprintf('Done. All figures are open as MATLAB figure windows.\n');
