

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Prox_Quad_Con.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Prox_Quad_Con is a function calculating the function and gradient 
% values of a quadratic prox-function.
%
% INPUT:
%
% x0          % initial point
% x           % current point
% pars        % structure including the parameteres
%
%     .q0     % value of prox-function Q in its center
%
% OUTPUT:
%
% qx          % value of Q at the current point x
% gqx         % grdient of Q at the current point x
% NLO         % number of linear operator used in Prox_Quad
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [qx, gqx, NLO] = Prox_Quad_Con(x0, x, pars)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Main body of Prox_Quad_Con.m %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    erorr('The number of input arguments is not enough');
elseif nargin == 2
    q0 = 0.5 * sqrt(eps + x0' * x0);
elseif nargin == 3
    q0 = pars.q0;
else
    erorr('The number of input arguments is inconsistent');
end

qx = q0 + 0.5 * sum(x.^2);

if nargout > 1
    gqx = x;
end  

NLO = 0;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% End of Prox_Quad_Con.m %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


