

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% OSGA_Init.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OSGA_Init is a function for initializing the parameters of OSGA. If
% some parameters specified by the user OSGA_Init uses these parameters.
% Otherwise, the default values will be employed.
%
% INPUT: 
%
% x0                   % initial point
% options              % structure including the parameteres of OSGA
%
%   .A                 % cell including matrices or linear operators
%   .MaxNumIter        % maximum number of iterations
%   .MaxNumFunEval     % maximum number of function evaluations
%   .MaxNumSubGradEval % maximum number of subgradient evaluations
%   .MaxNumLinOper     % maximum number of linear operators
%   .TimeLimit         % maximum running time
%   .StepLenLowBound   % lower bound for step size
%   .StepDiffLowBound  % lower bound for difference of two last points
%   .epsilon           % accuracy parameter
%   .f_target          % lower bound on function values
%   .delta             % tuning parameter for updating scheme PUS
%   .alpha_max         % maximum step size
%   .kappa             % tuning parameter for updating scheme PUS
%   .kappap            % tuning parameter for updating scheme PUS
%   .mu                % strong convexity parameter
%   .D                 % preconditioner matrix for the quadratic norm
%   .q0                % value of prox-function Q in its center
%   .xs                % optimizer (e.g., clean image or signal) 
%   .cons              % type of the constraint
%   .flag_MSE          % 1 : saves MSE
%                      % 0 : do not saves MSE     (default)
%   .flag_ISNR         % 1 : saves ISNR
%                      % 0 : do not saves ISNR    (default)
%   .flag_x_error      % 1 : saves x_error
%                      % 0 : do not saves x_error (default)
%   .flag_f_error      % 1 : saves f_error
%                      % 0 : do not saves f_error (default) 
%   .flag_time         % 1 : saves time for each iteration
%                      % 0 : do not saves time (default)  
%   .Stopping_Crit     % stopping criterion
%
%                      % 1 : stop if MaxNumIter is reached (default)
%                      % 2 : stop if MaxNumFunEval is reached
%                      % 3 : stop if MaxNumSubGradEval is reached
%                      % 4 : stop if MaxNumLinOper is reached
%                      % 5 : stop if TimeLimit is reached
%                      % 6 : stop if eta <= epsilon
%                      % 7 : stop if Norm_dx/max(1,Norm_x) <= epsilon
%                      % 8 : stop if fx <= f_target     
%
% OUTPUT:
%
% MaxNumIter           % maximum number of iterations
% MaxNumFunEval        % maximum number of function evaluations
% MaxNumSubGradEval    % maximum number of subgradient evaluations
% MaxNumLinOper        % maximum number of linear operators
% TimeLimit            % maximum running time
% StepLenLowBound      % lower bound for steplength
% StepDiffLowBound     % lower bound for difference of two last points
% epsilon              % accuracy parameter
% f_target             % lower bound on function values
% delta                % tuning parameter for updating scheme PUS
% alpha_max            % maximum step size
% kappa                % tuning parameter for updating scheme PUS
% kappap               % tuning parameter for updating scheme PUS
% mu                   % strong convexity parameter
% D                    % preconditioner matrix for the quadratic norm
% q0                   % value of prox-function Q in its center
% xs                   % clean image or signal
% x_opt                % optimizer 
% f_opt                % optimum 
% cons                 % type of the constraint
% flag_MSE             % 1 : saves MSE
%                      % 0 : do not saves MSE     (default)
% flag_ISNR            % 1 : saves ISNR
%                      % 0 : do not saves ISNR    (default)
% flag_x_error         % 1 : saves x_error
%                      % 0 : do not saves x_error (default)
% flag_f_error         % 1 : saves f_error
%                      % 0 : do not saves f_error (default)  
% flag_time            % 1 : saves time for each iteration
%                      % 0 : do not saves time (default)  
% Stopping_Crit        % stopping criterion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [delta,alpha_max,kappa,kappap,mu,D,q0,xs,x_opt,f_opt,epsilon, ...
  f_target,cons,MaxNumIter,MaxNumFunEval,MaxNumLinOper, ...
  MaxNumSubGradEval,StepLenLowBound,StepDiffLowBound,TimeLimit, ... 
  Stopping_Crit,flag_MSE, flag_ISNR, flag_x_error, flag_f_error, ...
  flag_time ] = OSGA_Init(x0,options)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Main body of OSGA_Init.m %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(x0);

if isfield(options,'delta') 
    delta = options.delta;
else
    delta = 0.9;
end

if isfield(options,'alpha_max') 
    alpha_max = options.alpha_max;
else
    alpha_max = 0.7;
end    

if isfield(options,'kappa') 
    kappa = options.kappa;
else
    kappa = 0.5;
end

if isfield(options,'kappap') 
    kappap = options.kappap;
else
    kappap = kappa;
end

if isfield(options,'mu') 
    mu = options.mu;
else
    mu = 0;
end

if isfield(options,'D') 
    D = options.B;
else
    %D1 = eye(n);
    D1 = speye(n);
    D  = {{D1, D1'},{}};
end

if isfield(options,'q0') 
    q0 = options.q0;
else
    Dx01 = feval(LOPMVM(D),x0,1);
    Dx0  = Dx01{1}{1};
    q0   = 0.5 * sqrt(eps + x0'*Dx0);
end

if isfield(options,'xs') 
    xs = options.xs;
else
    xs = ones(n,1);
end

if isfield(options,'x_opt') 
    x_opt = options.x_opt;
else
    x_opt = ones(n,1);
end

if isfield(options,'f_opt') 
    f_opt = options.f_opt;
else
    f_opt = 0;
end

if isfield(options,'epsilon') 
    epsilon = options.epsilon;
else
    epsilon = 10^(-8);
end

if isfield(options,'f_target') 
    f_target = options.f_target;
else
    f_target = -inf;
end

if isfield(options,'TimeLimit') 
    TimeLimit = options.TimeLimit;
else
    TimeLimit = inf;
end

if isfield(options,'MaxNumIter') 
    MaxNumIter = options.MaxNumIter;
else
    MaxNumIter = 5000;
end

if isfield(options,'MaxNumFunEval') 
    MaxNumFunEval = options.MaxNumFunEval;
else
    MaxNumFunEval = 10000;
end

if isfield(options,'MaxNumLinOper') 
    MaxNumLinOper = options.MaxNumLinOper;
else
    MaxNumLinOper = 15000;
end

if isfield(options,'MaxNumSubGradEval') 
    MaxNumSubGradEval = options.MaxNumSubGradEval;
else
    MaxNumSubGradEval = 10000;
end

if isfield(options,'StepLenLowBound') 
    StepLenLowBound = options.StepLenLowBound;
else
    StepLenLowBound = 1e-50;
end

if isfield(options,'StepDiffLowBound') 
    StepDiffLowBound = options.StepDiffLowBound;
else
    StepDiffLowBound = 1e-50;
end

if isfield(options,'Stopping_Crit') 
    Stopping_Crit = options.Stopping_Crit;
else
    Stopping_Crit = 1;
end

if isfield(options,'flag_MSE') 
    flag_MSE = options.flag_MSE;
else
    flag_MSE = 0;
end

if isfield(options,'flag_ISNR') 
    flag_ISNR = options.flag_ISNR;
else
    flag_ISNR = 0;
end

if isfield(options,'flag_x_error') 
    flag_x_error = options.flag_x_error;
else
    flag_x_error = 0;
end

if isfield(options,'flag_f_error') 
    flag_f_error = options.flag_f_error;
else
    flag_f_error = 0;
end

if isfield(options,'flag_time') 
    flag_time = options.flag_time;
else
    flag_time = 0;
end

if isfield(options,'cons') 
    cons = options.cons;
else
    cons = 'unconstrained';
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% End of OSGA_Init.m %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


