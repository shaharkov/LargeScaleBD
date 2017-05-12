%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Large-Scale Bounded Distortion Mappings".
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% init
rng(1)
clear
addpath('mex');


%% parameters
n = 100; % problem size
sigma = .2; % noise level (for initial map)
K = 1.5; % conformal distortion bound
lb = -1; % lower bound on SVs (-1 = disabled)
ub = -1; % upper bound on SVs (-1 = disabled)
iter_max = 1000; % maximal number of BD projection iterations
tol_err = 1e-10; % tolerance for stopping BD projection iterations
use_weighted_metric = false; % use a weighted metric?


%% generate problem
% generate a noisy surface
[x, y] = ndgrid(1:n,1:n);
V = [x(:), y(:)];
F = delaunay(V);
%V = V * [1 0 0.1; 0 1 0.2];
V = [V(:,1), V(:,2), 0.1*V(:,1)+2*n*(V(:,2)/n).^2.5];
V = V + sigma*randn(size(V));

% some constants
dim = size(F,2)-1;
n_vert = size(V,1);
n_tri = size(F,1);

% initial map
x0 = V(:,1:dim);

% setup linear constraints (fix centroid)
eq_lhs = kron(eye(dim),ones(1,n_vert))/n_vert;
eq_rhs = eq_lhs*colStack(x0);

% plot surface
figure;
patch('vertices',V,'faces',F,'facecolor','w');
axis equal; axis off;
cameratoolbar;
cameratoolbar('SetCoordSys','none');
title('surface');


%% solve problem
% setup BD solver
solver_bd = SolverProjectorBD(F, V, eq_lhs, eq_rhs, K, lb, ub, x0, SolverProjectorModeEnum.Tangent, use_weighted_metric);

% plot initial map
figure;
solver_bd.visualize();
title('Initial Map');
hA(1) = gca;

% run solver
solver_bd.solve(iter_max, tol_err); % solve BD projection

% plot output map
figure;
solver_bd.visualize();
title('Output Map');
hA(2) = gca;

linkaxes(hA);
