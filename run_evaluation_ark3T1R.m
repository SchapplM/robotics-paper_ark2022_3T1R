% Run the evaluation of the simulation results for the paper
% This creates the robot plots for the paper (see Readme)

% Moritz Schappler, schappler@imes.uni-hannover.de, 2021-01
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

this_path = fileparts( mfilename('fullpath') );
addpath(this_path);
run(fullfile(fileparts(which('ark2022_3T1R_dimsynth_data_dir.m')), 'dimsynth', ...
  'eval_figures_pareto.m')); close all;
run(fullfile(fileparts(which('ark2022_3T1R_dimsynth_data_dir.m')), 'dimsynth', ...
  'robot_names.m')); close all;
run(fullfile(fileparts(which('ark2022_3T1R_dimsynth_data_dir.m')), 'dimsynth', ...
  'eval_figures_pareto_groups.m')); close all;
run(fullfile(fileparts(which('ark2022_3T1R_dimsynth_data_dir.m')), 'dimsynth', ...
  'select_eval_robot_examples.m')); close all;
run(fullfile(fileparts(which('ark2022_3T1R_dimsynth_data_dir.m')), 'dimsynth', ...
  'robot_images.m')); close all;
