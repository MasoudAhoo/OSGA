

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% StopCrit.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% StopCrit is a function checking that one of the stopping criteria 
% holds to terminate OSGA. It also perepare the status determining why
% the algorithm is stoped.
%
% INPUT:
% 
% flag                  % 1: if xb ~= x_old
%                       % 0: if xb == x_old
% fx                    % function value at the current point x
% Niter                 % number of iterations 
% Nfunc                 % number of function evaluations
% Nsubgrad              % number of subgradient evaluations
% Nlinop                % number of linear operators
% eta                   %
% Norm_dx               % norm of xb - x_old
% Norm_xb               % norm of xb
% Time                  % running time
% alpha                 % step size
% x0                    % initial point
% MaxNumIter            % maximum number of iterations
% MaxNumFunEval         % maximum number of function evaluations
% MaxNumSubGradEval     % maximum number of subgradient evaluations
% MaxNumLinOper         % maximum number of linear operators
% TimeLimit             % maximum running time
% StepLenLowBound       % lower bound for step size
% StepDiffLowBound      % lower bound for difference of two last points
% epsilon               % accuracy parameter
% f_target              % lower bound on function values
% Stopping_Crit         % stopping criterion
%
%                       % 1 : stop if MaxNumIter is reached
%                       % 2 : stop if MaxNumFunEval is reached
%                       % 3 : stop if MaxNumSubGradEval is reached
%                       % 4 : stop if MaxNumLinOper is reached
%                       % 5 : stop if TimeLimit is reached
%                       % 6 : stop if eta <= epsilon
%                       % 7 : stop if Norm_dx / max(1,Norm_x) <= epsilon
%                       % 8 : stop if fx <= f_target     
%
% OUTPUT:
%
% StopFlag              % 1: if one of the stopping criteria holds
%                       % 0: if none of the stopping criteria holds
% Status                % the reason indicating why OSGA stopped
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


function [StopFlag, Status] = StopCrit ...
   (flag, fx, Niter, Nfunc, Nlinop, Nsubgrad, eta,Norm_dx, Norm_x, ... 
   Time,MaxNumIter, MaxNumFunEval, MaxNumLinOper, MaxNumSubGradEval, ...
   StepLenLowBound, StepDiffLowBound, TimeLimit, alpha, epsilon, ...
   f_target, Stopping_Crit)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Main body of StopCrit.m %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch Stopping_Crit
  case 1
    if Niter >= MaxNumIter
      StopFlag = 1;
      Status   = 'Maximum number of iterations is reached';
    else
      StopFlag = 0;
      Status   = [];
    end 
  case 2
    if Nfunc >= MaxNumFunEval
      StopFlag = 1;
      Status   = 'Maximum number of function evaluations is reached';
    else
      StopFlag = 0;
      Status   = [];
    end 
  case 3
    if Nsubgrad >= MaxNumSubGradEval
      StopFlag = 1;
      Status   = 'Maximum number of subgradient evaluations is reached';
    else
      StopFlag = 0;
      Status   = [];
    end
  case 4
    if Nlinop >= MaxNumLinOper
      StopFlag = 1;
      Status   = 'Maximum number of linear operators is reached';
    else
      StopFlag = 0;
      Status   = [];
    end 
  case 5
    if Time >= TimeLimit
      StopFlag = 1;
      Status   = 'Time limit is reached';
    else
      StopFlag = 0;
      Status   = [];
    end
  case 6
    if eta <= epsilon
      StopFlag = 1;
      Status   = 'eta <= epsilon';
    else
      StopFlag = 0;
      Status   = [];
    end 
  case 7
    if flag && Norm_dx / max(1,Norm_x) <= epsilon
      StopFlag = 1;
      Status = '||x_k - x_k-1|| / max(1,norm(x_k)) <= epsilon';
    else
      StopFlag = 0;
      Status   = [];
    end
  case 8
    if fx <= f_target
      StopFlag = 1;
      Status = 'f(x) <= f_target';
    else
      StopFlag = 0;
      Status   = [];
    end
end

if alpha <= StepLenLowBound
    StopFlag = 1;
    Status = 'Lower bound of step length is reached';
elseif flag && Norm_dx <= StepDiffLowBound
    StopFlag = 1;
    Status = '||x_k - x_k-1|| <= StepDiffLowBound';
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% End of StopCrit.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
