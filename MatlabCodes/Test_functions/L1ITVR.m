

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L1ITVR.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% L1ITVR is a function generating both function values and subgradients  
% of the convex test function:
% 
%                  f(x) = ||Ax-b||_1 + lambda ||x||_{ITV}.
%
% INPUT:
%
% xk         % current point;
% vxk        % vxk = A'x;
% opt        % structure includes required parameters;
%
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


function [fk gk1 gk2] = L1ITVR(opt,varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Main body of L1ITVR.m %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ~= 3
    error('The number of input arguments is not valid');
end   

if nargout ==0 || nargout ==2 || nargout > 3
    error('The number of output arguments is not valid');
end

xk  = varargin{1};
vxk = varargin{2};
   
% ========================= Function Evaluation ========================
b       = opt.b;
lambda  = opt.lambda;

Lxk     = length(xk);
size_Xk = [sqrt(Lxk) sqrt(Lxk)];

mat     = @(x) reshape(x,size_Xk);
vec     = @(x) reshape(x,Lxk,1);
xk      = mat(xk);
n       = length(xk);
Axk_b   = vxk{2}{1} - b;

ind     = [1:n-1];
DV      = xk(ind + 1,ind) - xk(ind,ind);
DH      = xk(ind,ind + 1) - xk(ind,ind);
NITVxk  = sum(sum(sqrt(DV.^2 + DH.^2)));
DVB     = abs(xk(ind + 1,n) - xk(ind,n));
DHB     = abs(xk(n,ind + 1) - xk(n,ind));
NITVxk  = NITVxk + sum(DVB) + sum(DHB);

fk      = sum(abs(Axk_b)) +  lambda * NITVxk;

%NITVxk = sum(sum(sqrt(diffh(xk).^2+diffv(xk).^2)));
%fk = sum(abs(Axk_b)) + lambda *  NITVxk;

% ====================== subgradient Evaluation ========================
if nargout > 1
    
    g = zeros(n,n);

    % ============================ g(1,1) ==============================
    DV     = xk(2,1) - xk(1,1);
    DH     = xk(1,2) - xk(1,1);
    g(1,1) = -(DV + DH)/(sqrt(DV^2 + DH^2) + realmin);

    % ================= g(1,j) for j = 2, 3, ..., n-1 ==================
    
    ind      = [2:n-1];
    DV       = xk(2,ind) - xk(1,ind);
    DH       = xk(1,ind + 1) - xk(1,ind);
    g(1,ind) = -(DV+DH)./(sqrt(DV.^2 + DH.^2)+realmin);
    
    DV       = xk(2,ind-1) - xk(1,ind-1);
    DH       = xk(1,ind) - xk(1,ind-1);
    g(1,ind) = g(1,ind)+DH./(sqrt(DV.^2 + DH.^2)+realmin);
    
    % ================= g(i,1) for i = 2, 3, ..., n-1 ==================
    
    ind      = [2:n-1];
    DV       = xk(ind+1,1) - xk(ind,1);
    DH       = xk(ind,2) - xk(ind,1);
    g(ind,1) = -(DV+DH)./(sqrt(DV.^2 + DH.^2)+realmin);
    
    DV       = xk(ind,1) - xk(ind-1,1);
    DH       = xk(ind-1,2) - xk(ind-1,1);
    g(ind,1) = g(ind,1)+DV./(sqrt(DV.^2 + DH.^2)+realmin);
    
    % ====== g(i,1) for i = 2, 3, ..., n-1 and j = 2, 3, ..., n-1 ======
    
    ind        = [2:n-1];
    DV         = xk(ind,ind) - xk(ind-1,ind);
    DH         = xk(ind-1,ind + 1) - xk(ind-1,ind);
    g(ind,ind) = DV./(sqrt(DV.^2 + DH.^2)+realmin);
    
    DV         = xk(ind+1,ind) - xk(ind,ind);
    DH         = xk(ind,ind + 1) - xk(ind,ind);
    g(ind,ind) = g(ind,ind)-(DV+DH)./(sqrt(DV.^2 + DH.^2)+realmin);
    
    DV         = xk(ind+1,ind-1) - xk(ind,ind-1);
    DH         = xk(ind,ind) - xk(ind,ind-1);
    g(ind,ind) = g(ind,ind)+DH./(sqrt(DV.^2 + DH.^2)+realmin);
       
    % ================= g(n,j) for j = 1, 2, ..., n-1 ==================
        
    ind      = [1:n-1];
    DH       = abs(xk(n,ind + 1) - xk(n,ind));
    g(n,ind) = -DH./(DH+realmin);
    
    DV       = xk(n,ind) - xk(n-1,ind);
    DH       = xk(n-1,ind + 1) - xk(n-1,ind);
    g(n,ind) = g(n,ind)+DV./(sqrt(DV.^2 + DH.^2)+realmin);   

    % ================= g(i,n) for i = 1, 2, ..., n-1 ==================
    
    ind      = [1:n-1];
    DV       = abs(xk(ind+1,n) - xk(ind,n));
    g(ind,n) = -DV./(DV+realmin);
    
    DV       = xk(ind+1,n-1) - xk(ind,n-1);
    DH       = xk(ind,n) - xk(ind,n-1);
    g(ind,n) = g(ind,n)+DH./(sqrt(DV.^2 + DH.^2)+realmin);
       
    % ============================ g(n,n) ==============================
    DV     = xk(n,n) - xk(n - 1,n);
    g(n,n) =  DV / (abs(DV) + realmin);
    
    DH     = xk(n,n) - xk(n,n - 1);
    g(n,n) = g(n,n) + DH / (abs(DV) + realmin);
    % ==================================================================
    
    g   = vec(g);
    
    gk1 = {{},{sign(Axk_b)}};
    gk2 = {{},lambda*g};

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% End of L1ITVR.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




