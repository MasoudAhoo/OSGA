

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Ridge_regression_driver.m %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Deblurring_driver is a script for restoring a deblurred/noisy image 
% using the following model
%
%             min  ||Ax-b||_2^2
%             s.t. ||x|| <= xi
%
% with OSGA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Main body of Ridge_regression_driver.m %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
format compact

% ======================= setting the parameters ======================= 
maxit = 500;
xi    = 15;  

% ======================= generating the problem =======================
% The function i_laplace.m is borrowed from the package 
% Regularization Tools

n        = 1000; 
[A1,b,x] = i_laplace(n);
xt       = x; 
b        = b + 0.1*rand;
A        = @(x) A1*x;
At       = @(x) A1'*x;
x0       = b;

% =============================== OSGA =================================
opt.b   = b;
func    = @ (varargin) L22R(opt,varargin{:});
prox    = @ (varargin) Prox_Quad_Con(varargin{:});
subprob = @ (varargin) Subebcon(varargin{:});

options.A                 = {{A, At},{}};
options.xs                = xt;
options.MaxNumIter        = maxit;
options.MaxNumFunEval     = 500;
options.MaxNumLinOper     = 300;
options.MaxNumSubGradEval = 500;
options.Stopping_Crit     = 1;
options.TimeLimit         = 5;
options.epsilon           = 1e-5;
options.mu                = 0;
options.flag_ISNR         = 0;
options.flag_x_error      = 1;
options.flag_f_error      = 0; 
options.flag_time         = 0; 
options.xi                = xi;
options.cons              = 'eball_constrained';

fprintf('Running OSGA ...\n')
t = tic;
[ x,f,out ] = OSGA( func,prox,subprob,x0,options );

T1       = toc(t);
x1       = x;
f1       = f;
F1       = out.F; 
x_error1 = out.x_error;

% ======================== illustrating outputs ========================

% ============================= figure 1 ===============================

vx0    = feval(LOPMVM({{A, At},{}}),x0,1);
f0     = func(x0,vx0);

fs = min(F1);
fs = fs-0.2*fs;

figure(1)
semilogy((F1-fs)/(f0 - fs),'-r','LineWidth',3) 
xlabel('iterations');
ylabel('function values');
legend('OSGA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% End of Ridge_regression_driver.m %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



