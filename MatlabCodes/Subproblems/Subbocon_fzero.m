

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Subuncon_fzero.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subbocon is a function computing the golbal maximizer of the rational  
% bound-constrained subproblem 
%            max   -(gamma + <h,x>)/Q(x)   (1)
%            s.t.  x in [xl,xu]
% of OSGA. The problem is solved by finding a zero of a nonlinear
% equation.
% 
% INPUT:
%  
% parms
%
%      .x        % current point
%      .h        % a parameter of the subproblem
%      .gamma    % a parameter of the subproblem
%      .x0       % initial point
%      .q0       % center of prox-function
%      .xl       % lower bound on x
%      .xu       % upper bound on x
%
% OUTPUT:
%
% u        % global maximizer of the subproblem (1)
% Egh      % optimimum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ u, Egh ] = Subbocon_fzero(h,gamma,x0,q0,xl,xu )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Main body of Subuncon_fzero.m %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ~= 6
    erorr('The number of input arguments is inconsistent');
end

Nh2  = sum(h.^2);
beta = gamma + h' * x0;
    
if beta <= 0
    e = (-beta + sqrt(beta^2 + 2 * q0 * Nh2)) / (2*q0);
else
    e = Nh2 / (beta + sqrt(beta^2 + 2 * q0 * Nh2));
end
u  = x0 - h/e;

PC = @(varargin) projc(h,x0,xl,xu,varargin{:});

n  = length(x0);
a  = sum((u >= xl).*(u <= xu));
if (n-a)/n >= 0.9
    Egh = e;
else
    fun = @(varargin) subfunc(gamma,q0,h,x0,varargin{:});
    lambda2 = 1/e;
    lambda_b = fzero(fun,lambda2);
      
    if lambda_b >= 0
        Egh = 1/lambda_b;
        u   = PC(lambda_b);
    else
        u        = max(xl,min(xu,u));
        e        = -(gamma+sum(h.*u))/(0.5*sum(u.^2)+q0);
        Egh      = e;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  f = subfunc(gamma,q0,h,x0,varargin)
    lambda = varargin{1};
    p      = PC(lambda);
    p_x0   = p-x0;
    f      = (1/lambda)*(0.5*sum(p_x0.^2) + q0) + gamma + sum(h.*p);      
end

function  p = projc(h,x0,xl,xu,varargin)
    lambda = varargin{1};
    y = x0-lambda*h;
    p = max(xl,min(xu,y));     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% End of Subuncon_fzero.m %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

