

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OSGA.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [ x,f,out ] = OSGA(func,prox,subprob,x0,options)
% OSGA, Optimal SubGradient Algorithm, is an optimal algorithm developed  
% to approximately minimize the unconstrained convex function
%                        min   f(x)  
%                        s.t.  x in X
% where X is a finite-dimensional real vector space.
%
% INPUT:
%
% func                 % function handle for the objective function
% prox                 % prox-function
% subprob              % function handle to solve the OSGA subproblem
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
%   .xs                % clean image or signal
%   .x_opt             % optimizer 
%   .f_opt             % optimum
%   .cons              % type of the constraint
%   .flag_MSE          % 1 : saves MSE
%                      % 0 : do not saves MSE     (default)
%   .flag_ISNR         % 1 : saves ISNR
%                      % 0 : do not saves ISNR    (default)
%   .flag_x_error      % 1 : saves x_error
%                      % 0 : do not saves x_error (default)
%   .flag_f_error      % 1 : saves f_error
%                      % 0 : do not saves f_error (default)
%   .flag_time         % 1 : saves f_error
%                      % 0 : do not saves f_error (default)
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
% x                    % the best approximation of the optimizer
% f                    % the best approximation of the optimum
% out                  % structure including more output information
%
%   .T                 % running time
%   .Niter             % total number of iterations
%   .Nfunc             % total number of function evaluations
%   .Nsubgrad          % total number of subgradient evaluations
%   .Nlinop            % total number of employed linear operators
%   .F                 % array including all function values            
%   .MSE               % mean squared error
%   .ISNR              % improvement signal to the noise
%   .x_error           % relative error norm(xb(:) - xs(:))/norm(xs)
%   .f_error           % relative error (fk - fs)/(f0 - fs))    
%   .Status            % reason indicating why OSGA stopped
%
% REFERENCES: 
%
% [1] M. Ahookhosh, Optimal subgradient algorithms 
%     with application to large-scale linear inverse problems, (2014) 
%     http://arxiv.org/abs/1402.7291.
%
% [2] M. Ahookhosh, A. Neumaier, An optimal subgradient algorithm 
%     for large-scale bound-constrained convex optimization, (2015)
%     http://arxiv.org/abs/1501.01497.
%
% [3] M. Ahookhosh, A. Neumaier, An optimal subgradient algorithm 
%     for large-scale convex optimization in simple domains, (2015)
%     http://arxiv.org/abs/1501.01497.
%             
% [4] A. Neumaier, OSGA: A fast subgradient algorithm  
%     with optimal complexity, (2014)
%     http://arxiv.org/abs/1402.1125.
%             
% WRITTEN BY: 
%
% Masoud Ahookhosh
% Faculty of Mathematics, University of Vienna, Austria
%
% LAST UPDATE: 
%
% January 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ x,f,out ] = OSGA(func,prox,subprob,x0,options)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Initializing and setting the parameters %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long;

% ================ Error messages for input and output =================
if nargin > 5
    erorr('The number of input arguments is more than what is needed');
elseif nargin < 4
    erorr('The number of input arguments is not enough');
end;

if isempty(func)
    error('OSGA needs the function handle func to be defined');
elseif ~isa(func,'function_handle')
    error('func should be a function handle');
end

if isempty(prox)
    error('OSGA needs the function handle prox to be defined');
elseif ~isa(prox,'function_handle')
    error('prox should be a function handle');
end

if isempty(subprob)
    error('OSGA needs the function handle subprob to be defined');
elseif ~isa(subprob,'function_handle')
 error('subprob should be a function handle');
end

if isempty(x0)
    error('OSGA needs an starting point x0 to be defined');
elseif ~isa(x0,'numeric')
    error('x0 should be a numeric vector');
end

% ======================================================================
A     = options.A;
n     = length(x0);
xb    = x0; 
Niter = 0;

% ======================= Computing f0 and g0 ==========================
vxb             = feval(LOPMVM(A),xb,1);
Nlinop          = 1;
[fxb gxb1 gxb2] = func(xb,vxb);
gxb             = SubGradEval(A, gxb1, gxb2);
Nfunc           = 1;
Nsubgrad        = 1;
Nlinop          = Nlinop + 1;
f0              = fxb;
F               = f0;
options.f0      = f0;

% =================== Initializing the parameters ======================
[delta,alpha_max,kappa,kappap,mu,D,q0,xs,x_opt,f_opt,epsilon, ...
  f_target,cons,MaxNumIter,MaxNumFunEval,MaxNumLinOper, ...
  MaxNumSubGradEval,StepLenLowBound,StepDiffLowBound,TimeLimit, ...
  Stopping_Crit,flag_MSE,flag_ISNR,flag_x_error,flag_f_error, ...
  flag_time] = OSGA_Init(x0,options);
% ===== User has requested viewing the default values of "options" =====
if flag_MSE == 1
    MSE(1) = 1/n * sum((xb - xs).^2);
end

if flag_ISNR == 1
    Nxs_x0 = sum((x0 - xs).^2);
    ISNR(1) = 0;
end

if flag_x_error == 1
    Nxs        = sqrt(sum(x_opt(:).^2));
    x_error(1) = sqrt(sum((x0(:) - x_opt(:)).^2))/Nxs;
end

if flag_f_error == 1
    f_error(1) = 1;
end

if flag_time == 1
    Time(1) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Main body of OSGA.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0 = tic;

if fxb <= f_target 
    fprintf('fxb <= f_target');
else 
    
    pars.D         = D;
    pars.q0        = q0;
    [qxb,gqxb,NLO] = prox(x0,xb,pars);
    Nlinop         = Nlinop + NLO;
    h              = gxb - mu * gqxb;
    gamma          = fxb - mu * qxb - h' * xb;
    gammab         = gamma - fxb;
    
    % ==================== Solving the subproblem ======================
    switch cons
        case 'unconstrained'
            [u, Egh] = subprob(D,h,gammab,x0,q0);
        case 'bound_constrained'
            xl = options.xl;
            xu = options.xu;
            [u, Egh] = subprob(h,gammab,x0,q0,xl,xu);
        case 'eball_constrained'
            xi = options.xi;
            [u, Egh] = subprob(h,gammab,q0,xi);
        case 'nonnegativity_constrained'
            [u, Egh] = subprob(h,gammab,q0);
    end
    
    eta      = Egh - mu; 
    alpha    = alpha_max;
    StopFlag = 0;
   
    % ===================== Start of the main loop =====================
    while ~StopFlag
        
        x            = xb + alpha * (u - xb);
        vx           = feval(LOPMVM(A),x,1);
        Nlinop       = Nlinop + 1;
        [fx gx1 gx2] = func(x,vx);
        gx           = SubGradEval(A, gx1, gx2);
        Nsubgrad     = Nsubgrad + 1;
        Nlinop       = Nlinop + 1;
        Nfunc        = Nfunc + 1;
        
        [qx,gqx,NLO] = prox(x0,xb,pars);
        %Nlinop       = Nlinop + NLO;
        g            = gx - mu * gqx; 
        h_bar        = h + alpha * (g - h);
        gamma_bar    = gamma + alpha * (fx - mu * qx - g' * x - gamma); 
        
        if fx <= fxb
            xpb  = x;
            fxpb = fx;
        else
            xpb  = xb;
            fxpb = fxb;
        end
        
        gammapb = gamma_bar - fxpb;
        
        % =================== Solving the subproblem ===================
        parms.gamma = gammapb;
        parms.h     = h_bar;
        switch cons
            case 'unconstrained'
                [up, Egh] = subprob(D,h_bar,gammapb,x0,q0);
            case 'bound_constrained'
                [up, Egh] = subprob(h_bar,gammapb,x0,q0,xl,xu);
            case 'eball_constrained'
                [up, Egh] = subprob(h_bar,gammapb,q0,xi);
            case 'nonnegativity_constrained'
                [up, Egh] = subprob(h_bar,gammapb,q0);
        end
                
        xp     = xpb + alpha * (up - xpb);
        vxp    = feval(LOPMVM(A),xp,1);
        Nlinop = Nlinop + 1;
        fxp    = func(xp,vxp);
        Nfunc  = Nfunc + 1;
        
        if fxp <= fxpb
            xb_bar  = xp;
            fxb_bar = fxp;
        else
            xb_bar  = xpb;
            fxb_bar = fxpb;
        end
        
        gammab_bar = gamma_bar - fxb_bar;
        
        % =================== Solving the subproblem ===================       
        switch cons
            case 'unconstrained'
                [u_bar, Egh] = subprob(D,h_bar,gammab_bar,x0,q0);
            case 'bound_constrained'
                [u_bar, Egh] = subprob(h_bar,gammab_bar,x0,q0,xl,xu);
            case 'eball_constrained'
                [u_bar, Egh] = subprob(h_bar,gammab_bar,q0,xi);
            case 'nonnegativity_constrained'
                [u_bar, Egh] = subprob(h_bar,gammab_bar,q0);
        end
               
        % ================== Accepting the new point ===================
        eta_bar = Egh - mu;
        x_old   = xb;
        xb      = xb_bar;
        %xb = max(xl,min(xu,xb));
        %sum(xb<0) + sum(xb>255)
        Niter   = Niter + 1;
        fxb     = fxb_bar; 
        
        % =============== Gathering output information =================
        F(Niter+1) = fxb;
        if flag_time == 1
            Time(Niter+1) = toc(T0);
        end
                       
        if flag_MSE == 1
            MSE(Niter+1) = 1/n * sum((xb - xs).^2);
        end
        
        if flag_ISNR == 1
            ISNR(Niter+1) = 10*log10(Nxs_x0/sum((xb(:) - xs(:)).^2));          
        end
        
        if flag_x_error == 1
            Nx_opt = norm(x_opt);
            x_error(Niter+1) = sqrt(sum((xb(:)-x_opt(:)).^2))/Nx_opt;          
        end

        if flag_f_error == 1
            f_error(Niter+1) = (fxb-f_opt)/(f0-f_opt);          
        end
        
        % =================== Updating parameters ======================
        if fxb <= f_target 
            fprintf('fxb <= f_target');
            break;
        else
            [alpha,h,gamma,eta,u] = ...
                Update_Pars(delta,kappa,kappap,eta,alpha_max,alpha, ...
                            h,gamma,u,h_bar,gamma_bar,eta_bar,u_bar);
        end
        
        % ================ Checking stopping criteria ==================
        if xb == x_old
            flag = 0;
        else 
            flag = 1;
        end  
        
        Time1    = toc(T0);
        %Norm_dx = sqrt(sum((xb - x_old).^2));
        %Norm_xb = sqrt(sum(xb.^2));
        Norm_dx = 100;
        Norm_xb = 100;
                        
        [StopFlag, Status] = ...
         StopCrit(flag,fxb,Niter,Nfunc,Nlinop,Nsubgrad,eta, ...
         Norm_dx,Norm_xb,Time1,MaxNumIter,MaxNumFunEval, ...
         MaxNumLinOper,MaxNumSubGradEval,StepLenLowBound, ...
         StepDiffLowBound,TimeLimit,alpha,epsilon,f_target,Stopping_Crit);  
        
    end
    % ======================= End of main loop =========================
    
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Status
x            = xb;
f            = fxb;
T            = toc(T0);
out.T        = T;
out.F        = F';
out.Niter    = Niter;
out.Nfunc    = Nfunc;
out.Nsubgrad = Nsubgrad;
out.Nlinop   = Nlinop;
out.Status   = Status;

if flag_MSE == 1
    out.MSE = MSE;
end
        
if flag_ISNR == 1
    out.ISNR = ISNR;          
end
        
if flag_x_error == 1
    out.x_error = x_error;          
end

if flag_f_error == 1
    out.f_error = f_error;          
end

if flag_time == 1
    out.Time = Time;
    size(Time)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% End of OSGA.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


