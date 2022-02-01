% Show properties of the kinematics model for 3T1R PKM
% Equation numbers in the comments refer to the numbers in the paper
% 
% Parts of the code of this file are taken from ParRob_nullspace_proof.m
% of https://github.com/SchapplM/robotics-paper_icinco2021
% (ICINCO-Paper "Singularity Avoidance of Task-Redundant Robots in Pointing
% Tasks: On Nullspace Projection and Cardan Angles as Orientation
% Coordinates, Schappler et al. 2021)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-01
% (C) Institut für mechatronische Systeme, Leibniz Universität Hannover

clear
clc
%% Initialization
if isempty(which('serroblib_path_init.m'))
  error('Serial robot library not in search path.');
end
if isempty(which('parroblib_path_init.m'))
  error('Parallel robot library not in search path.');
end
this_dir = fileparts(which('ParRob_3T1R_taskred_mdl_proof.m'));
addpath(fullfile(this_dir, '..'));
if isempty(which('ark2022_3T1R_dimsynth_data_dir'))
  error(['You have to create a file ark2022_3T1R_dimsynth_data_dir pointing to the ', ...
    'directory containing the results of the dimensional synthesis']);
end
resdirtotal = ark2022_3T1R_dimsynth_data_dir();
datadir = fullfile(fileparts(which('ark2022_3T1R_dimsynth_data_dir.m')),'data');
if ~exist(fullfile(datadir, 'robot_groups.mat'), 'file')
  error('Run eval_figures_pareto_groups.m first');
end
tmp = load(fullfile(datadir, 'robot_groups.mat'));
RobotGroups = tmp.RobotGroups;
% Select 4-RRRRR for evaluation. For other robots, see variable RobotGroups
GroupName = 'P4RRRRR10G';
fprintf('%d Groups available. Selected no. %d (%s)\n', size(RobotGroups,1), ...
  find(strcmp(RobotGroups(:,1), GroupName)), GroupName);

%% Initialize Robot and load data
% Load the data from the dimensional synthesis from before
detailfile = fullfile(datadir, sprintf('detail_result_group_%s.mat', GroupName));
if ~exist(detailfile, 'file')
  error(['Results of dimensional synthesis are missing. ', ...
    'Run select_eval_robot_examples.m first']);
end
res = load(detailfile);
R = res.R;
Q = res.Q;
parroblib_addtopath({R.mdlname});
q0 = Q(1,:)';
X0 = res.X(1,:)';

%% Compute IK for Initial Pose
% Set excluded Z Euler angle to zero as this must not have influence on
% the model
X0(6) = 0;
% The given value for q0 should already be correct. Compute IK using Newton
% Raphson algorithm (text below equ. 6) to check again.
[q, phi] = R.invkin_ser(X0, q0);
assert(all(abs(q-q0) < 1e-8), ['The initial value q0 from dimensional ', ...
  'synthesis should already have been correct for x0']);
% Full kinematic constraints according to equ. 5 (delta)
% and redunced task-redundant constraints psi (equ. 10)
[psi, delta] = R.constr3(q, X0);
assert(all(size(delta)==[24,1]), 'Output 2 of constr3 has to be 24x1');
assert(all(abs(psi) < 1e-10), ['residual psi corresponding to the task-', ...
  'redundant case has to be zero']);
% Obtain the reduced constraints from equ. 6. alpha_z1 ist entry 4 (first
% entry of the rotational residual of equ. 4).
I_red = [1:3, 5:6, 7:length(delta)];
delta_red = delta(I_red); % delta_t1, delta_r1red, delta_2...4
assert(all(size(delta_red)==[23,1]), 'delta_red has to be 23x1');
assert(all(abs(delta_red) < 1e-10), 'delta_red has to be zero');
% Set the last Euler angle to the actual value. This is necessary for
% some parts of the implementation
X0_corr = X0;
% Since phi_z was set to zero before, the entry delta(4) (alpha_z1 of 
% equ. 3) contains the actual phi_z. This is different from the one given
% in the input pose X0 due to the task redundancy.
X0_corr(6) = X0_corr(6) + delta(4);
% Recompute the full constraints with corrected phi_z.
[~, delta_corr] = R.constr3(q, X0_corr);
assert(abs(delta_corr(4)) < 1e-10, ['residual corresponding to the ', ...
  'Z Euler angle has to be zero after correction']);

%% Plot result
figure(10);clf;
hold on; grid on;
xlabel('x in m'); ylabel('y in m'); zlabel('z in m');
view(3);
s_plot = struct( 'ks_legs', [R.I1L_LEG; R.I2L_LEG], 'ks_platform', 1:6, 'straight', 0);
R.plot(q, X0, s_plot );
title(R.mdlname, 'interpreter', 'none');

%% Differential kinematics (Sect. III)
% TODO:
% Compute matrix psi_dq and delta_dq from first paragraph of Sect. III
[psi_dq,delta_dq] = R.constr3grad_q(q, X0);
fprintf('Direct kinematics matrix psi_dq is %dx%d and has condition %1.1f\n', ...
  size(psi_dq,1), size(psi_dq,2), cond(psi_dq));
deltared_dq = delta_dq(I_red,:);
fprintf(['Overconstraint direct kinematics matrix deltared_dq is %dx%d and ', ...
  'has condition %1.1e\n'], size(deltared_dq,1), size(deltared_dq,2), cond(deltared_dq));

% Compute corresponding matrix deltared_dy and psi_dy
[psi_dx,delta_dx] = R.constr3grad_x(q, X0);
psi_dy = psi_dx(:,1:3);
deltared_dy = delta_dx(I_red,R.I_EE_Task);

% Inverse Jacobian corresponding to all joints of the robot. See last
% sentence of first paragraph of Sect. III
Jinvtilde_x = -delta_dq\delta_dx(:,R.I_EE);
Jinv = Jinvtilde_x(R.I_qa,:); % Obtain inverse manipulator Jacobian
assert(all(size(Jinv)==[4 4]), 'inverse Jacobian Matrix has to be 4x4');
fprintf('Group %s. %s: Rank of the %dx%d Jacobian: %d/%d. condition number: %1.2e\n', ...
  GroupName, R.mdlname, size(Jinv,1), size(Jinv,2), rank(Jinv), sum(R.I_EE), cond(Jinv));
% Obtain Jacobian matrix by inversion, task-redundant Jacobian by removing
% rows. From hereon, see ICINCO-Paper (ref. see comment head)
J_x = inv(Jinv);
J_y = J_x(R.I_EE_Task,:);
% Nullspace projector from paper approach and Jacobian approach
N_y   = eye(sum(R.I_qa)) - pinv(J_y)*  J_y;
N_Psi = eye(R.NJ) -        pinv(psi_dq)*psi_dq;
% Remove numeric noise to simplify evaluation
N_Psi(abs(N_Psi)<1e-10) = 0;
N_y(abs(N_y)<1e-10) = 0;

% check the relation between both nullspace projectors.
J_q_qa = Jinvtilde_x * J_x; %#ok<MINV> 
N_Psi_from_y = J_q_qa * N_y * J_q_qa';
N_Psi_from_y(abs(N_Psi_from_y)<1e-10) = 0;
% The following assertion proves that both approaches lead to the same
% nullspace motion
N_Psi_ratio = N_Psi ./ N_Psi_from_y;
k_Psi = N_Psi_ratio(1,1);
I_notnan = ~isnan(N_Psi_ratio);
assert(all(abs(N_Psi_ratio(I_notnan(:))-k_Psi)<1e-5), ['All components of the nullspace', ...
  'projector in full coordinates have to be in constant ratio to its ', ...
  'derivation from actuator coordinates']);
