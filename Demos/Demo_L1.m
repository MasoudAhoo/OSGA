

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L1_driver.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% L1_driver is a script for reconstructing a sparse signal from noisy  
% data using the following model
%
%                 min  ||Ax-b||_2^2 + lambda ||x||_1
%                 s.t. x in R^n
%
% with OSGA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Main body of L1_driver.m %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

% ==================== generating a random problem =====================
%flag = 'sparse';
flag = 'dense';
m    = 5000; 
n    = 10000;

if strcmp(flag , 'dense')
   A  = rand(m,n);
   b  = rand(m,1);
   x0 = rand(n,1);       
elseif strcmp(flag , 'sparse')
   A  = sprand(m,n,0.05);
   b  = sprand(m,1,0.05);
   x0 = sprand(n,1,0.05);
end

% ============= setting the parameters and executing OSGA ==============       
opt.b      = b;
opt.lambda = 1;
func       = @ (varargin) L22L1R(opt,varargin{:});
prox       = @ (varargin) Prox_Quad(varargin{:});
subprob    = @ (varargin) Subuncon(varargin{:});

options.A                 = {{A, A'},{}};
options.MaxNumIter        = 500;
options.MaxNumFunEval     = 1000;
options.MaxNumLinOper     = 1500;
options.MaxNumSubGradEval = 1000;
options.Stopping_Crit     = 5;
options.TimeLimit         = 30;
options.epsilon           = 1e-5;
options.mu                = 0;
options.cons              = 'unconstrained';

fprintf('Running OSGA ...\n')
tic
[ x,f,out ] = OSGA( func,prox,subprob,x0,options );         
toc

% ============ illustrating function values vs. iterations =============
F = out.F;

figure(1)
loglog(F,'k','LineWidth',3)
xlabel('iterations');
ylabel('function values');
legend('OSGA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% End of L1_driver.m %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



