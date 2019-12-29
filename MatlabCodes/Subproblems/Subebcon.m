

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subebcon.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subbocon is a function computing the golbal maximizer of the rational  
% unconstrained subproblem 
%            max   -(gamma + <h,x>)/Q(x)   (1)
%            s.t.  ||x|| <= xi
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


function [ u, Egh ] = Subebcon( h,gamma,q0,xi )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Main body of Subebcon.m %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ~= 4
    error('The number of input arguments is inconsistent');
end

beta1 = q0;
beta2 = gamma;
beta3 = -0.5*sum(h.^2);
   
if beta2 <= 0
   Egh = (-beta2 + sqrt(beta2^2 - 4*beta1*beta3))/(2*beta1);
else
   Egh = -2*beta3/(beta2 + sqrt(beta2^2 - 4*beta1*beta3));
end

if norm(Egh^(-1)*h) <= xi
    u = - Egh^(-1)*h;
else
    nh  = norm(h);
    u   = - (xi/nh)*h; 
    Egh = -2*(gamma+xi*nh)/(xi^2+2*q0);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% End of Subebcon.m %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


