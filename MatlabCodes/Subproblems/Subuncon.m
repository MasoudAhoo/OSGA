

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subuncon.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subbocon is a function computing the golbal maximizer of the rational  
% unconstrained subproblem 
%            max   -(gamma + <h,x>)/Q(x)   (1)
%            s.t.  x in R^n
% of OSGA.
% 
% INPUT:
% 
% x        % current point
% h        % a parameter of the subproblem
% gamma    % a parameter of the subproblem
% x0       % initial point
% q0       % center of prox-function
%
% OUTPUT:
%
% u        % global maximizer of the subproblem (1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ u, Egh ] = Subuncon( D,h,gamma,x0,q0 )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Main body of Subuncon.m %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ~= 5
    error('The number of input arguments is inconsistent');
end

D1   = D{1}{1};
d    = D1 \ h;
Nh2  = h' * d;
beta = gamma + h' * x0;
    
if beta <= 0
    Egh = (-beta + sqrt(beta^2 + 2 * q0 * Nh2)) / (2 * q0);
else
    Egh = Nh2 / (beta + sqrt(beta^2 + 2 * q0 * Nh2));
end 
    
u = x0 - Egh^(-1) * d;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% End of Subuncon.m %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

