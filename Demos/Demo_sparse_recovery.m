

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Sparse_recovery_driver.m %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sparse_recovery_driver is a script for reconstructing a sparse signal   
% from noisy data using the following model
%
%                 min  ||Ax-b||_2^2 + lambda ||x||_1
%                 s.t. x in R^n
%
% with OSGA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Main body of Sparse_recovery_driver.m %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

% ==================== generating a random problem =====================
% The data is generated the same as those reported in the SpaRSA package
n                 = 2000;                    % original signal length
k                 = 1000;                    % number of observations
n_spikes          = floor(0.03*n);           % number of spikes
ff                = zeros(n,1);    
q                 = randperm(n);             % random +/- 1 signal
ff(q(1:n_spikes)) = sign(randn(n_spikes,1));
A                 = randn(k,n);              % measurement matrix
A                 = orth(A')';               % orthonormalize rows  
sigma             = 0.001;
y                 = A*ff + sigma*randn(k,1); % noisy observations
A1                = @(x) A*x;
A1t               = @(x) A'*x;
x0                = ones(n,1);

% ============= setting the parameters and executing OSGA ==============       
opt.b      = y;
opt.lambda = 0.1*max(abs(A'*y));
func       = @ (varargin) L22L1R(opt,varargin{:});
prox       = @ (varargin) Prox_Quad(varargin{:});
subprob    = @ (varargin) Subuncon(varargin{:});

options.A                 = {{A1, A1t},{}};
options.MaxNumIter        = 500;
options.MaxNumFunEval     = 1000;
options.MaxNumLinOper     = 1500;
options.MaxNumSubGradEval = 1000;
options.Stopping_Crit     = 5;
options.TimeLimit         = 10;
options.epsilon           = 1e-5;
options.flag_MSE          = 1;
options.mu                = 0;
options.xs                = ff;
options.cons              = 'unconstrained';

fprintf('Running OSGA ...\n')

tic
[ x,f,out ] = OSGA( func,prox,subprob,x0,options );         
toc

% ======== illustrating function values and MSE vs. iterations =========
F   = out.F;
MSE = out.MSE;

figure(1)
loglog(F,'k','LineWidth',3)
xlabel('iterations');
ylabel('function values');
legend('OSGA');

figure(2)
semilogx(MSE,'k','LineWidth',3)
xlabel('iterations');
ylabel('MSE');
legend('OSGA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% End of Sparse_recovery_driver.m %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



