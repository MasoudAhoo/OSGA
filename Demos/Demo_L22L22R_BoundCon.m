

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% L22L22R_BonCon_driver.m %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% L22L22R_BonCon_driver is a script for solving the following 
% minimization problem
%
%             min  ||Ax-b||_2^2 + lambda ||x||_2^2
%             s.t. x in [xl xu]
%
% with OSGA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Main body of L22L22R_BonCon_driver.m %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
format compact

% ==================== generating a random problem =====================
% The function i_laplace.m is borrowed from the package 
% Regularization Tools

n    = 2000;
xl   = 0.05*ones(n,1);
xu   = 0.95*ones(n,1);

[A,b,x] = i_laplace(n);
xt      = x; 
b       = b +0.1*rand;
x0      = max(xl,min(xu,b));

% ======== setting the parameters and executing the algorithms ========= 

lambda = 1;
maxit  = 100;

opt.b      = b;
opt.lambda = lambda;
func       = @ (varargin) L22L22R(opt,varargin{:});

% =============================== OSGA =================================
prox       = @ (varargin) Prox_Quad_Con(varargin{:});

options.A                 = {{A, A'},{}};
options.MaxNumIter        = maxit;
options.MaxNumFunEval     = 500;
options.MaxNumLinOper     = 300;
options.MaxNumSubGradEval = 500;
options.Stopping_Crit     = 1;
options.TimeLimit         = 1;
options.epsilon           = 1e-5;
options.flag_ISNR         = 0;
options.flag_x_error      = 0;
options.flag_f_error      = 0; 
options.flag_time         = 0; 
options.xl                = xl;
options.xu                = xu;
options.mu                = 0;
options.cons               = 'bound_constrained';

mex CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" Subbocon.cpp
subprob    = @ (varargin) Subbocon(varargin{:});
fprintf('Running OSGA-1 ...\n')
t = tic;
[ x,f,out ] = OSGA( func,prox,subprob,x0,options );

T1 = toc(t);
x1 = x;
f1 = f;
F1 = out.F;  

subprob    = @ (varargin) Subbocon_fzero(varargin{:});
fprintf('Running OSGA-2 ...\n')
t = tic;
[ x,f,out ] = OSGA( func,prox,subprob,x0,options );

T2       = toc(t);
x2       = x;
f2       = f;
F2       = out.F;  

% ======================== illustrating outputs ========================

% ============================= figure 1 ===============================
vx0    = feval(LOPMVM({{A, A'},{}}),x0,1);
f0     = func(x0,vx0);

fs=min([F1;F2;]);
fs=fs-0.001*fs;

figure(7)
semilogy((F1-fs)/(f0 - fs),'k','LineWidth',2) 
hold on
semilogy((F2-fs)/(f0 - fs),'--b','LineWidth',2)
xlabel('iterations');
ylabel('function values');
legend('OSGA-1','OSGA-2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% End of L22L22R_BonCon_driver.m %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



