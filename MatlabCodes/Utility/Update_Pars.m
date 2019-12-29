

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update_Pars.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Update_Pars is a function for updating the parameters of OSGA
%
% INPUT:
%
% delta       % global tuning parameter
% kappa       % global tuning parameter
% kappap      % global tuning parameter
% alpha_max   % maximum step size
% h_bar       % parameter determined in OSGA
% gamma_bar   % parameter determined in OSGA
% eta_bar     % parameter determined in OSGA
% u_bar       % search direction
%
% OUTPUT:  
%
% alpha       % step size
% h           % OSGA's parameter
% gamma       % OSGA's parameter
% eta         % OSGA's parameter
% u           % current search direction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [alpha,h,gamma,eta,u] = ...
         Update_Pars(delta,kappa,kappap,eta,alpha_max,alpha, ...
         h,gamma,u,h_bar,gamma_bar,eta_bar,u_bar)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Main body of Update_Pars.m %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = (eta-eta_bar)/(delta*alpha*eta);

if R < 1
    alpha_bar = alpha*exp(-kappa);
else
    alpha_bar = min(alpha*exp(kappap*(R-1)),alpha_max);
end

alpha = alpha_bar;
if eta_bar < eta
    h     = h_bar;
    gamma = gamma_bar;
    eta   = eta_bar;
    u     = u_bar;
end
   
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% End of Update_Pars.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

