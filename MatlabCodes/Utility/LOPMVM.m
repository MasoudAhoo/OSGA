

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOPMVM.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LOPMVM is a function to implement matrix-vector multiplications or
% forward and adjoint operations. 
% 
% INPUT:
%
% A % cell array of the following form
%   %        
%   %     A = {{A1,...,An1,A1',...,An1'},{W1,...,Wn2,W1',...,Wn2'}}
%   %
%   %  where V = {A1,..., An1,A1', ..., An1'} is a cell array including 
%   %  numeric matrices or linear operators of smooth part of objective 
%   %  function, and W = {W1, ..., Wn2, W1', ..., Wn2'} is a cell array 
%   %  including numeric matrices or linear operators of nonsmooth part  
%   %  of the objective function
%       
% OUTPUT: 
%
% op % function handle to implement matrix-vector multiplication or
%    % forward and adjoint operations
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function op = LOPMVM( A )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Main body of LOPMVM.m %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part checks if the cell array A is appropriately constructed or
% not and also compute length of V and W.

cell_length = length(A);

if cell_length ~= 2
    error('A should be a cell array with the dimension 2');
else
    cell_length1 = length(A{1}); 
    if mod(cell_length1,2) ~= 0
        error('Length of the cell should be an even number');
    else
        clength1 = cell_length1/2;
    end
    
    cell_length2 = length(A{2});
    if mod(cell_length2,2) ~= 0
        error('Length of the cell should be an even number');
    else
        clength2 = cell_length2/2;
    end
end

% ===== Check if elements of cell arrays A{1} and A{2} are numeric =====

flag1 = 1;
for i = 1 : cell_length1
    flag1 = flag1 * isa(A{1}{i}, 'numeric');
end

flag2 = 1;
for i = 1 : cell_length2
    flag2 = flag2 * isa(A{2}{i}, 'numeric');
end

if flag1*flag2 == 1
    op = @(x,mode) numeric_matrix_vec( A, x, mode );    
else 
    % Check if elements of cell arrays A{1} and A{2} are function handle
    flag3 = 1;
    for i = 1 : cell_length1
        flag3 = flag3 * isa(A{1}{i}, 'function_handle');
    end

    flag4 = 1;
    for i = 1 : cell_length2
        flag4 = flag4 * isa(A{2}{i}, 'function_handle');
    end
    
    if flag3*flag4 == 1
        op = @(x,mode) forward_adjoint_linop( A, x, mode );
    else
    error( ...
       'elments of A{1} and A{2} must be matrices or function handles');
    end    
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% numeric_matrix_vec.m %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   
% numeric_matrix_vec implements matrix-vector multiplications when all 
% elements of V and W are matrixes.
%
% INPUT:
%
% A        % cell array involving matrices
% x        % current point for matrix-vector multiplications
% mode     % 1 : forward multiplication Ax
%          % 2 : adjoint multiplication A'x
%
% OUTPUT:
%
% y        % cell array contains y = Ax or y = A'x
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


function y = numeric_matrix_vec( A, x, mode )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Main body of numeric_matrix_vec.m %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch mode

    case 1
        if isa(x, 'numeric')
            if ~isempty(A{1})
                y1 = cell(1,clength1);
                for j = 1 : clength1
                    y1{j} = A{1}{j} * x;
                end
            else
                y1 = {};
            end
            
            if ~isempty(A{2})
                y2 = cell(1,clength2);
                for j = 1 : clength2
                    y2{j} = A{2}{j} * x;
                end
            else
                y2 = {};
            end
            
        else
            error('x has to be a numeric vector');
        end
        y = {y1 y2};

    case 2 
        if (isa(x, 'cell')) && (length(x) == 2)
           % y1 is a vector containing sum of gradients for smoorh part
           % of the objective function
           if ~isempty(A{1})
               y1 = cell(1,clength1);
               %y1 = A{1}{clength1 + 1} * x{1}{1};
               for j = 1 : clength1
                   y1{j} = A{1}{clength1 + j} * x{1}{j};
               end
           else
               y1 = {};
           end
            
           % y1 is a vector containing sum of subgradients for nonsmooth
           % part of the objective function
           if ~isempty(A{2})
               y2 = cell(1,clength2);
               %y2 = A{2}{clength2 + 1} * x{2}{1};
               for j = 1 : clength2
                   y2{j} = A{2}{clength2 + j} * x{2}{j};
                   %y2 = y2 + A{2}{clength2 + j} * x{2}{j};
               end
           else
               y2 = {};
           end
        else
            error('x has to be a cell including two cell arrays');
        end
        y = {y1 y2};
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% End of numeric_matrix_vec.m %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% forward_adjoint_linop.m %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    
% forward_adjoint_linop implements forward or adjoint operatoions when  
% all elements of V and W are linear operators.
%
% INPUT:
%
% A        % cell array involving linear operators
% x        % current point for forward or adjoint operatoions
% mode     % 1 : forward operation A(x)
%          % 2 : adjoint operation A*(x)
%
% OUTPUT:
%
% y        % cell array contains y = A(x) or y = A*(x)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


function y = forward_adjoint_linop( A, x, mode )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main body of forward_adjoint_linop.m %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
switch mode
    case 1
        if isa(x, 'numeric')
            if ~isempty(A{1})
                y1 = cell(1,clength1);
                for j = 1 : clength1
                    y1{j} = A{1}{j}(x);
                end
            else
                y1 = {};
            end
            
            if ~isempty(A{2})
                y2 = cell(1,clength2);
                for j = 1 : clength2
                    y2{j} = A{2}{j}(x);
                end
            else
                y2 = {};
            end
        else
            error('x has to be a numeric vector');
        end
        y = {y1 y2};
    case 2 
        if (isa(x, 'cell')) && (length(x) == 2)
           % y1 is a vector containing sum of gradients for smoorh part
           % of the objective function
           if ~isempty(A{1})
               y1 = cell(1,clength1);
               %y1 =  A{1}{clength1 + 1}(x{1}{1});
               for j = 1 : clength1
                   y1{j} = A{1}{clength1 + j}(x{1}{j});
                   %y1 = y1 + A{1}{clength1 + j}(x{1}{j});
               end
           else
               y1 = {};
           end
            
           % y1 is a vector containing sum of subgradients for nonsmooth 
           % part of the objective function
           if ~isempty(A{2})
               y2 = cell(1,clength2);
               %y2 = A{2}{clength2 + 1}(x{2}{1});
               for j = 1 : clength2
                   y2{j} = A{2}{clength2 + j}(x{2}{j});
                   %y2 = y2 + A{2}{clength2 + j}(x{2}{j});
               end
           else
               y2 = {};
           end
        else
            error('x has to be a cell including two cell arrays');
        end
        y = {y1 y2};
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% End of forward_adjoint_linop.m %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% End of Subfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% End of LOPMVM.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
