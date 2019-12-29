

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OSGA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OSGA - Optimal SubGradient Algorithm for unconstrained  or  
%        constrained convex optimization problems. It approximately 
%        solves the smooth or nonsmooth convex optimization problem 
%                        min   f(x),  
%                        s.t.  x in C,
%        where C is a simple constraint. It is assumed that the  
%        first-order oracle is available that rerturns for a given 
%        x in C a function value f(x) and a subgradient g(x).
%
% Software written in 2014 by 
%
%    Masoud Ahookhosh, 
%    Faculty of Mathematics, University of Vienna, Austria,
%
% based on ideas related to optimal first-order methods and traditional
% subgradient methods.
% 
% Free for academic use, provided you cite the accompanying paper
% in all publications using results obtained with this package.
%
% The driver programs L1_driver.m and Deblurring_driver.m demonstrate 
% the usage of the solver.
%
% The interface to each subprogram is fully documented in the 
% corresponding file. Please also see the OSGA user manual.
%
% LAST UPDATE:
% 
% January 2015
%
% REFERENCES:
% 
% [1] M. Ahookhosh, Optimal subgradient algorithms 
%     with application to large-scale linear inverse problems, (2014) 
%     http://arxiv.org/abs/1402.7291.
%
% [2] M. Ahookhosh, A. Neumaier, An optimal subgradient algorithm 
%     for large-scale bound-constrained convex optimization, (2015)
%     http://arxiv.org/abs/1501.01497.
%
% [3] M. Ahookhosh, A. Neumaier, An optimal subgradient algorithm 
%     for large-scale convex optimization in simple domains, (2015)
%     http://arxiv.org/abs/1501.01497.
%             
% [4] A. Neumaier, OSGA: A fast subgradient algorithm  
%     with optimal complexity, (2014)
%     http://arxiv.org/abs/1402.1125.
%
% ACKNOWLEDGMENTS: 
%
% I am very grateful to Arnold Neumaier for encouraging me and 
% instructing me to write this software package.  
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% contents %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Programs in this toolbox:
%
%
% MAIN PROGRAM:
%
% OSGA.m                    % the local solver
%
% DRIVERS:
%                      
% Deblurring_driver.m       % shows how to use OSGA for solving 
%                           % deblurring problem with the isotropic 
%                           % total variation
% Denoising_driver.m        % shows how to use OSGA for solving 
%                           % denoising problem with the isotropic 
%                           % total variation
% Inpainting_driver.m       % shows how to use OSGA for solving 
%                           % inpainting problem with the isotropic 
%                           % total variation
% L1_driver.m               % shows how to use OSGA for solving 
%                           % L1-regularized problem
% Sparse_recovery_driver.m  % shows how to use OSGA for solving bound-            
%                           % constrained problems
% Deblurring_non_driver.m   % shows how to use OSGA for solving 
%                           % nonnegativity constrained problems
% L22L22R_BounCon_driver.m  % shows how to use OSGA for sparse recovery
%                           % problems
% Ridge_regression_driver.m % shows how to use OSGA for solving Ridge
%                           % regression problems
%       
% AUXILIARY PROGRAMS: 
%
% OSGA_Init.m        % initializes the OSGA parameters
% Update_Pars.m      % updates parameters of OSGA
% Prox_Quad.m        % calculates functions and gradients of a quadratic 
%                    % prox-function
% Subuncon.m         % solves unconstrained OSGA's subproblem
% Subnocon.m         % solves nonnegativity constrained OSGA's 
%                    % subproblem
% Subbocon.cpp       % solves bound-constrained OSGA's subproblem 
%                    % exactly
% Subbocon_fzero.m   % solves bound-constrained OSGA's subproblem 
%                    % inexactly
% Subebcon.m         % solves $l_2$-ball constrained OSGA's subproblem               
% LOPMVM.m           % implements matrix-vector multiplications or 
%                    % forward and adjoint operators
% SubGradEval.m      % calculates the subgradients
% StopCrit.m         % checks stopping criteria
% L22ITVR.m          % computes function values and subgradients of the
%                    % isotropic total variation regularized problem
% L22L1R.m           % computes function values and subgradients of the
%                    % L1-regularized problem
% L22L22R.m          % computes function values and subgradients of the
%                    % L22-regularized problem
% L22R.m             % computes function values and subgradients of the
%                    % L22 problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

