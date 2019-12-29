

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Deblurring_driver.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Deblurring_driver is a script for restoring a deblurred/noisy image 
% using the following model
%
%             min  ||Ax-b||_2^2 + lambda ||x||_{ITV}
%             s.t. x in R^n
%
% with OSGA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Main body of Deblurring_driver.m %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

% ====================== reading the clean image =======================
colormap(gray);
%X_exact = double(imread('Cameraman_512.tif'));
X_exact = double(imread('Goldhill_512.gif'));
%X_exact = double(imread('Elaine_512.tif'));
SizeX   = size(X_exact);
M       = SizeX(1);
N       = SizeX(2);

% ================ generating the blurred/noisy image ==================
% The data is generated the same as those reported in the TwIST package
middle = N/2+1;                                % N must be even;
B      = zeros(N);                             % Uniform blurr matrix;
lx     = 4;                                    % blur x-size;
ly     = 4;                                    % blurr y-size;
B((middle - ly):(middle + ly), (middle - lx):(middle + lx)) = 1;
B      = fftshift(B);                          % circularly center;
B      = B/sum(sum(B));                        % normalize;
y      = real(ifft2(fft2(B).*fft2(X_exact)));  % convolve;
BSNR   = 40;                                   % set BSNR;
Py     = var(y(:));
sigma  = sqrt((Py/10^(BSNR/10)));

y      = y + sigma*randn(N);                   % add noise;                     

% =============== generating the convolution operator ==================
K    = fft2(B);
KC   = conj(K);

% Function handles for OSGA convolution operators (with matrix input)
%A   = @(x) real(ifft2(K.*fft2(x)));
%At  = @(x) real(ifft2(KC.*fft2(x)));

% Function handles for OSGA convolution operators (with vector input)
mat  = @(X) reshape(X,SizeX);
vec  = @(X) reshape(X,M*N,1);
A    = @(x) vec(real(ifft2(K.*fft2(mat(x)))));
At   = @(x) vec(real(ifft2(KC.*fft2(mat(x)))));
       
% ============= setting the parameters and executing OSGA ============== 
x0         = vec(y);      
opt.b      = vec(y);
opt.lambda = 5e-2;
func       = @ (varargin) L22ITVR(opt,varargin{:});
prox       = @ (varargin) Prox_Quad(varargin{:});
subprob    = @ (varargin) Subuncon(varargin{:});

options.A                 = {{A, At},{}};
options.xs                = vec(X_exact);
options.MaxNumIter        = 20;
options.MaxNumFunEval     = 500;
options.MaxNumLinOper     = 300;
options.MaxNumSubGradEval = 500;
options.Stopping_Crit     = 5;
options.TimeLimit         = 5;
options.epsilon           = 1e-5;
options.mu                = 0;
options.flag_ISNR         = 1;
options.flag_x_error      = 1;
options.flag_f_error      = 1; 
options.cons              = 'unconstrained';

fprintf('Running OSGA 1 ...\n')
[ x,f,out ] = OSGA( func,prox,subprob,x0,options );
x1       = x;
f1       = f;
F1       = out.F;  
ISNR1    = out.ISNR;
x_error1 = out.x_error;
f_error1 = out.f_error;
PSNR1    =20*log10(255*sqrt(numel(X_exact))/norm(mat(x)-X_exact,'fro'));

options.TimeLimit         = 15;
fprintf('Running OSGA 2 ...\n') 
[ x,f,out ] = OSGA( func,prox,subprob,x0,options );
x2          = x;
f2          = f;
F2          = out.F;  
ISNR2       = out.ISNR;
x_error2    = out.x_error;
f_error2    = out.f_error;
PSNR2    =20*log10(255*sqrt(numel(X_exact))/norm(mat(x)-X_exact,'fro'));

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
%%%%%%%%%%%%%%%%%%%%% End of Deblurring_driver.m %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



