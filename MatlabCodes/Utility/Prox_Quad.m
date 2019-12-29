

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Prox_Quad.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Prox_Quad is a function calculating the function and gradient values 
% of a quadratic prox-function.
%
% INPUT:
%
% x0          % initial point
% x           % current point
% pars        % structure including the parameteres
%
%     .q0     % value of prox-function Q in its center
%     .D      % preconditioner matrix for the quadratic norm
%
% OUTPUT:
%
% qx          % value of Q at the current point x
% gqx         % grdient of Q at the current point x
% NLO         % number of linear operator used in Prox_Quad
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [qx, gqx, NLO] = Prox_Quad(x0, x, pars)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Main body of Prox_Quad.m %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    erorr('The number of input arguments is not enough');
elseif nargin == 2
    q0 = 0.5 * sqrt(eps + x0' * x0);
elseif nargin == 3
    q0 = pars.q0;
    D  = pars.D;
else
    erorr('The number of input arguments is inconsistent');
end

x_x0   = (x - x0);
Dx_x01 = feval(LOPMVM(D),x_x0,1);
Dx_x0  = Dx_x01{1}{1};
NLO    =  1;
qx     = q0 + 0.5 * sqrt(x_x0' * Dx_x0);

if nargout > 1
    gqx = Dx_x0;
end    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% End of Prox_Quad.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


