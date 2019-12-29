

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subnocon.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


function [ u, Egh ] = Subnocon( h,gamma,q0 )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Main body of Subnocon.m %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ~= 3
    error('The number of input arguments is inconsistent');
end

hm    = min(0,h);
beta1 = q0;
beta2 = gamma;
beta3 = 0.5*sum(hm.^2) - sum(h.*hm);
   
if beta2 <= 0
   Egh = (-beta2 + sqrt(beta2^2 - 4*beta1*beta3))/(2*beta1);
else
   Egh = -2*beta3/(beta2 + sqrt(beta2^2 - 4*beta1*beta3));
end 
   
u = - Egh^(-1)*hm;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% End of Subnocon.m %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


