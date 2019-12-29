

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SubGradEval.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SubGradEval is a function to construct gradient or subgradient of a 
% convex optimization problem. 
% 
% INPUT:
%
% A   % cell array of the following form   
%     %    
%     %    A = {{A1,...,An1,A1',...,An1'},{W1,...,Wn2,W1',...,Wn2'}}
%     %
%     % where V = {A1,...,An1,A1',...,An1'} is a cell array including
%     % numeric matrixes or linear operators of smooth part 
%     % of objective function, and W = {W1, ..., Wn2, W1', ..., Wn2'} 
%     % is a cell array including numeric matrixes or linear operators 
%     % of smooth part of the objective function
%     
% gx1 % cell array of the following form
%     %      
%     %    gx1 = {{gs1, ..., gsn1}, {gns1, ..., gnsn2}},
%     %
%     % where  gx1{1} = {gs1, ..., gsn1} is a cell array including 
%     % gradients of smooth functions including affine terms,  and 
%     % gx1{2} = {gns1, ..., gnsn1}     is a cell array containing
%     % subgradients of nonsmooth functions involving affine terms
%     
% gx2 % cell array of the form
%     %
%     %    gx2 = {gs, gns},
%     %
%     % where gx2{1} = gs is the gradient of the smooth part of f which 
%     % is not including any affine term and     gx2{2} = gns    is the
%     % subgradient of the nonsmooth part of f including no affine term
%       
% OUTPUT: 
%
% g   % gradient or subgradient of the objective function in x
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function g = SubGradEval( A, gx1, gx2 )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Main body of SubGradEval.m %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IEgx1_1 = isempty(gx1{1});
IEgx1_2 = isempty(gx1{2});
IEgx2_1 = isempty(gx2{1});
IEgx2_2 = isempty(gx2{2});

Lgx1    = feval(LOPMVM(A),gx1,2);

if     ~IEgx1_1 && IEgx1_2 && IEgx2_1 && IEgx2_2

    g_length    = length(Lgx1{1}{1});
    g           = zeros(g_length,1);

    Lgx1_length1 = length(Lgx1{1});
    for i = 1 : Lgx1_length1
        g = g + Lgx1{1}{i};
    end
    
elseif IEgx1_1 && ~IEgx1_2 && IEgx2_1 && IEgx2_2

    g_length     = length(Lgx1{2}{1});
    g            = zeros(g_length,1);

    Lgx1_length2 = length(Lgx1{2});
    for i = 1 : Lgx1_length2
        g = g + Lgx1{2}{i};
    end
       
elseif IEgx1_1 && IEgx1_2 && ~IEgx2_1 && IEgx2_2
    
    g = gx2{1};
    
elseif IEgx1_1 && IEgx1_2 && IEgx2_1 && ~IEgx2_2
    
    g = gx2{2};
     
elseif ~IEgx1_1 && ~IEgx1_2 && IEgx2_1 && IEgx2_2
    
    g_length    = length(Lgx1{1}{1});
    g           = zeros(g_length,1);

    Lgx1_length1 = length(Lgx1{1});
    for i = 1 : Lgx1_length1
        g = g + Lgx1{1}{i};
    end

    Lgx1_length2 = length(Lgx1{2});
    for i = 1 : Lgx1_length2
        g = g + Lgx1{2}{i};
    end
    
elseif ~IEgx1_1 && IEgx1_2 && ~IEgx2_1 && IEgx2_2
    
    g_length    = length(Lgx1{1}{1});
    g           = zeros(g_length,1);

    Lgx1_length1 = length(Lgx1{1});
    for i = 1 : Lgx1_length1
        g = g + Lgx1{1}{i};
    end

    g = g + gx2{1};
    
elseif ~IEgx1_1 && IEgx1_2 && IEgx2_1 && ~IEgx2_2

    g_length    = length(Lgx1{1}{1});
    g           = zeros(g_length,1);

    Lgx1_length1 = length(Lgx1{1});
    for i = 1 : Lgx1_length1
        g = g + Lgx1{1}{i};
    end

    g = g + gx2{2};
    
elseif IEgx1_1 && ~IEgx1_2 && ~IEgx2_1 && IEgx2_2
    
    g_length     = length(Lgx1{2}{1});
    g            = zeros(g_length,1);

    Lgx1_length2 = length(Lgx1{2});
    for i = 1 : Lgx1_length2
        g = g + Lgx1{2}{i};
    end

    g = g + gx2{1};
    
elseif IEgx1_1 && ~IEgx1_2 && IEgx2_1 && ~IEgx2_2
    
    g_length     = length(Lgx1{2}{1});
    g            = zeros(g_length,1);

    Lgx1_length2 = length(Lgx1{2});
    for i = 1 : Lgx1_length2
        g = g + Lgx1{2}{i};
    end

    g = g + gx2{2};
    
elseif IEgx1_1 && IEgx1_2 && ~IEgx2_1 && ~IEgx2_2 
    
    g = gx2{1} + gx2{2};
    
elseif IEgx1_1 && ~IEgx1_2 && ~IEgx2_1 && ~IEgx2_2
    
    g_length     = length(Lgx1{2}{1});
    g            = zeros(g_length,1);

    Lgx1_length2 = length(Lgx1{2});
    for i = 1 : Lgx1_length2
        g = g + Lgx1{2}{i};
    end

    g = g + gx2{1} + gx2{2};
    
elseif ~IEgx1_1 && IEgx1_2 && ~IEgx2_1 && ~IEgx2_2
    
    g_length    = length(Lgx1{1}{1});
    g           = zeros(g_length,1);

    Lgx1_length1 = length(Lgx1{1});
    for i = 1 : Lgx1_length1
        g = g + Lgx1{1}{i};
    end

    g = g + gx2{1} + gx2{2};
    
elseif ~IEgx1_1 && ~IEgx1_2 && IEgx2_1 && ~IEgx2_2
    
    g_length    = length(Lgx1{1}{1});
    g           = zeros(g_length,1);

    Lgx1_length1 = length(Lgx1{1});
    for i = 1 : Lgx1_length1
        g = g + Lgx1{1}{i};
    end

    Lgx1_length2 = length(Lgx1{2});
    for i = 1 : Lgx1_length2
        g = g + Lgx1{2}{i};
    end

    g = g + gx2{2};
    
elseif ~IEgx1_1 && ~IEgx1_2 && ~IEgx2_1 && IEgx2_2
    
    g_length    = length(Lgx1{1}{1});
    g           = zeros(g_length,1);

    Lgx1_length1 = length(Lgx1{1});
    for i = 1 : Lgx1_length1
        g = g + Lgx1{1}{i};
    end

    Lgx1_length2 = length(Lgx1{2});
    for i = 1 : Lgx1_length2
        g = g + Lgx1{2}{i};
    end

    g = g + gx2{1};
    
elseif ~IEgx1_1 && ~IEgx1_2 && ~IEgx2_1 && ~IEgx2_2
    
    g_length    = length(Lgx1{1}{1});
    g           = zeros(g_length,1);

    Lgx1_length1 = length(Lgx1{1});
    for i = 1 : Lgx1_length1
        g = g + Lgx1{1}{i};
    end

    Lgx1_length2 = length(Lgx1{2});
    for i = 1 : Lgx1_length2
        g = g + Lgx1{2}{i};
    end

    g = g + gx2{1} + gx2{2};
    
else
    
    error('The objective function is empty')
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% End of SubGradEval.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



