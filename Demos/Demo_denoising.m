

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Denoising_driver.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Denoising_driver is a script for restoring a noisy image using the
% following model
%
%             min  ||Ax-b||_2^2 + lambda ||x||_{ITV}
%             s.t. x in R^n
%
% with OSGA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Main body of Denoising_driver.m %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

% ====================== reading the clean image =======================
colormap(gray);
X_exact = double(imread('Pirate_1024.tif'));
%X_exact = double(imread('Cameraman_512.tif'));
%X_exact = double(imread('Goldhill_512.gif'));
X_exact = X_exact/255;
SizeX   = size(X_exact);
M       = SizeX(1);
N       = SizeX(2);

% ================ generating the blurred/noisy image ==================
if ~license('test','communication_toolbox')
    awgn = @(x,snr,ignore) x + ...
        10^( (10*log10(sum(abs(x(:)).^2)/length(x(:))) - snr)/20 ) ...
        *randn(size(x));
end

randn('state',sum(double('denoising test')) );
SNR     = 15;
X_noisy = awgn( X_exact, SNR, 'measured' );
y       = X_noisy;            

% ======================= generating operators =========================
mat  = @(X) reshape(X,SizeX);
vec  = @(X) reshape(X,M*N,1);
I    = @(x) x;
It   = @(x) x;
       
% ============= setting the parameters and executing OSGA ============== 
x0         = vec(y);      
opt.b      = vec(y);
opt.lambda = 5e-2;
func       = @ (varargin) L22ITVR(opt,varargin{:});
prox       = @ (varargin) Prox_Quad(varargin{:});
subprob    = @ (varargin) Subuncon(varargin{:});

options.A                 = {{I, It},{}};
options.xs                = vec(X_exact);
options.MaxNumIter        = 20;
options.MaxNumFunEval     = 500;
options.MaxNumLinOper     = 300;
options.MaxNumSubGradEval = 500;
options.Stopping_Crit     = 1;
options.TimeLimit         = 5;
options.epsilon           = 1e-5;
options.mu                = 0;
options.flag_MSE          = 1;
options.flag_ISNR         = 1;
options.flag_x_error      = 1;
options.flag_f_error      = 1; 
options.cons              = 'unconstrained';

fprintf('Running OSGA 1 ...\n')
[ x,f,out ] = OSGA( func,prox,subprob,x0,options );
x1          = x;
f1          = f;
F1          = out.F;  
ISNR1       = out.ISNR;
x_error1    = out.x_error;
f_error1    = out.f_error;
PSNR1        =20*log10(sqrt(numel(X_exact))/norm(mat(x)-X_exact,'fro'));

options.MaxNumIter        = 50;
fprintf('Running OSGA 2 ...\n') 
[ x,f,out ] = OSGA( func,prox,subprob,x0,options );
x2          = x;
f2          = f;
F2          = out.F;  
ISNR2       = out.ISNR;
x_error2    = out.x_error;
f_error2    = out.f_error;
PSNR2        =20*log10(sqrt(numel(X_exact))/norm(mat(x)-X_exact,'fro'));

% ======================== illustrating outputs ========================

% ============================= figure 1 ===============================
figure(1);
colormap gray; 
splot = @(n) subplot(2,2,n);

splot(1);
imagesc(X_exact); axis off;
title('Original image');
drawnow;

splot(2);
imagesc(y); axis off;
title(sprintf('Noisy image'))
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
splot = @(n) subplot(2,2,n);

splot(1);
semilogy(F1,'c','LineWidth',4) 
hold on
semilogy(F2,'--b','LineWidth',2)
xlabel('iterations');
ylabel('function values');
legend('OSGA-1', 'OSGA-2');

splot(2);
semilogy(f_error1,'c','LineWidth',4) 
hold on
semilogy(f_error2,'--b','LineWidth',2)
xlabel('iterations');
ylabel('f-error');
legend('OSGA-1', 'OSGA-2');

splot(3);
plot(x_error1,'c','LineWidth',4) 
hold on
plot(x_error2,'--b','LineWidth',2)
xlabel('iterations');
ylabel('x-error');
legend('OSGA-1', 'OSGA-2');

splot(4);
semilogy(ISNR1,'c','LineWidth',4) 
hold on
semilogy(ISNR2,'--b','LineWidth',2)
xlabel('iterations');
ylabel('ISNR');
legend('OSGA-1', 'OSGA-2');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% End of Denoising_driver.m %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



