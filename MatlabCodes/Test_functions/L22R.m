

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L22R.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% L22L1R is a function generating both function values and subgradient  
% evaluations of the convex test function:
% 
%                  f(x) = 1/2 ||Ax-b||_2^2.
%
% INPUT:
%
% xk         % current point;
% vxk        % vxk = A'x;
% opt        % structure includes required parameters;
%    .b      % observation;
%    .lambda % regularization parameter;
%
% OUTPUT:
%
% fk  % function value of f at xk
%
% gk1 % cell array of the following form
%     %      
%     %    gk1 = {{gs1, ..., gsn1}, {gns1, ..., gnsn2}},
%     %
%     % where   gk1{1} = {gs1, ..., gsn1} is a cell array including 
%     % gradients of smooth functions including affine terms,   and 
%     % gk1{2} = {gns1, ..., gnsn1}      is a cell array containing
%     % subgradients of nonsmooth functions involving affine terms;
%     
% gk2 % cell array of the form
%     %
%     %    gk2 = {gs, gns},
%     %
%     % where gk2{1} = gs is the gradient of the smooth part of  f which 
%     % is not including any affine term and     gk2{2} = gns     is the
%     % subgradient of the nonsmooth part of f including no affine term; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [fk gk1 gk2] = L22R(opt,varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Main body of L22R.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ~= 3
    error('The number of input arguments is not valid');
end   

if nargout == 0 || nargout == 2 || nargout > 3
    error('The number of output arguments is not valid');
end

xk  = varargin{1};
vxk = varargin{2};
   
% ======================== Function evaluation =========================
b      = opt.b;
svxk   = vxk{1}{1};
Axk_b  = svxk - b;
fk     = 0.5 * sum(Axk_b.^2);
    
% ====================== Subgradient evaluation ========================
if nargout > 1
    gk1 = {{Axk_b},{}};
    gk2 = {{},{}};
end 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of L22R.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

