

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Deblurring_non_driver.m %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Deblurring_non_driver is a script for restoring a deblurred/noisy  
% image using the following model
%
%             min  ||Ax-b||_2^2 + lambda ||x||_{ITV}
%             s.t. x in R^n  
% or 
%
%             min  ||Ax-b||_2^2 + lambda ||x||_{ITV}
%             s.t. x >= 0  
%
% with OSGA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Main body of Deblurring_non_driver.m %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
format compact

% ====================== reading the clean image =======================
colormap(gray);
Xe    = double(imread('MR_brain256.tif'));
Xc    = Xe/255;
[m,n] = size(Xc);
xl    = zeros(m*n,1);
xu    = inf*ones(m*n,1);

% == generating the blurred/noisy image and the convolution operator ===
% ============================= ADMM data ==============================
H      = fspecial('gaussian',[9,9],3);
sigma  = 1e-3;
center = [5,5];
mat   = @(x) reshape(x,[m,n]);
vec   = @(x) reshape(x,m*n,1);
y     = double(imfilter(Xc,H,'circular','conv')) + sigma*randn(m,n);
x0    = max(0.03,vec(y)); 
A     = @(x) vec(imfilter(mat(x),H,'circular','conv'));
At    = @(x) vec(imfilter(mat(x),H,'circular','corr'));

% ============= setting the parameters and executing OSGA ============== 

lambda = 5e-4;
maxit  = 100;

opt.b      = vec(y);
opt.lambda = lambda;
func       = @ (varargin) L22ITVR(opt,varargin{:});

% =============================== OSGA =================================
options.A                 = {{A, At},{}};
options.xs                = vec(Xc);
options.MaxNumIter        = maxit;
options.MaxNumFunEval     = 500;
options.MaxNumLinOper     = 300;
options.MaxNumSubGradEval = 500;
options.Stopping_Crit     = 1;
options.TimeLimit         = 5;
options.epsilon           = 1e-5;
options.mu                = 0;
options.flag_ISNR         = 1;
options.flag_f_error      = 1;
options.xl                = xl;
options.xu                = xu;

options.cons              = 'unconstrained';
prox                      = @ (varargin) Prox_Quad(varargin{:});
subprob                   = @ (varargin) Subuncon(varargin{:});

fprintf('Running OSGA (unconstrained) ...\n')
t = tic;
[ x,f,out ] = OSGA( func,prox,subprob,x0,options );

T1       = toc(t);
x1       = x;
f1       = f;
f_error1 = out.f_error;  
ISNR1    = out.ISNR;
PSNR1    = 20*log10(1*sqrt(numel(Xc))/norm(mat(x1)-Xc,'fro'));

options.cons              = 'nonnegativity_constrained';
prox                      = @ (varargin) Prox_Quad_Con(varargin{:});
subprob                   = @ (varargin) Subnocon(varargin{:});

fprintf('Running OSGA (nonnegativity constrained) ...\n')
t = tic;
[ x,f,out ] = OSGA( func,prox,subprob,x0,options );

T2       = toc(t);
x2       = x;
f2       = f;
f_error2 = out.f_error; 
ISNR2    = out.ISNR;
PSNR2    = 20*log10(1*sqrt(numel(Xc))/norm(mat(x2)-Xc,'fro'));

% ======================== illustrating outputs ========================

% ============================= figure 1 ===============================

figure(1);
colormap gray; 
splot = @(n) subplot(2,2,n);

splot(1);
imagesc(Xc); axis off;
title('Original image');
drawnow;

splot(2);
imagesc(y); axis off;
title(sprintf('Blurred/noisy image'))
drawnow;

splot(3);
x1 = mat(x1);
imagesc(x1); axis off;
title(sprintf('OSGA (f = %2.2f, PSNR = %2.2f)',f1,PSNR1))
drawnow;

splot(4);
x2 = mat(x2);
imagesc(x2); axis off;
title(sprintf('OSGA (f = %2.2f, PSNR = %2.2f)',f2,PSNR2))
drawnow;

% ============================= figure 2 ===============================
figure(2)
splot = @(n) subplot(1,2,n);

splot(1);
semilogy(f_error1,'k','LineWidth',2) 
hold on
semilogy(f_error2,'--b','LineWidth',2)
xlabel('iterations');
ylabel('f-error');
legend('OSGA-1', 'OSGA-2');

splot(2);
plot(ISNR1,'k','LineWidth',2) 
hold on
plot(ISNR2,'--b','LineWidth',2)
xlabel('iterations');
ylabel('ISNR');
legend('OSGA-1', 'OSGA-2');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% End of Deblurring_non_driver.m %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



