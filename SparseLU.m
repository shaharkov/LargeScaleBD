%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Large-Scale Bounded Distortion Mappings".
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef SparseLU < handle
    properties
        LHS;
        L;
        U;
        p;
        q;
    end
    methods
        function obj = SparseLU(LHS)
            if ~issparse(LHS)
                error('Argument must be a sparse matrix')
            end              
            obj.LHS = LHS;
            [obj.L,obj.U, obj.p, obj.q] = lu(LHS, 'vector');
        end
        function x = solve(obj,RHS)
            x(obj.q,:) = obj.U\(obj.L\(RHS(obj.p,:)));
        end
    end
    
end

