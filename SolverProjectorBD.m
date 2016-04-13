%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Large-Scale Bounded Distortion Mappings".
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef SolverProjectorBD < SolverProjector
    
    properties
        F;
        K;
        lb;
        ub;
        dim;
        distortions;
        flips;
        minsv;
        maxsv;
        use_weighted_metric;
    end
    
    methods
        function obj = SolverProjectorBD(F, V, eqLHS, eqRHS, K, lb, ub, x0, mode, use_weighted_metric)
            % setup BD solver
            obj.report(1,'-----------------------------------------------------------------\n');
            obj.report(1,'Constructing SolverProjectorBD (vertices %d,  elements %d)\n', size(V,1), size(F,1));
            t_start = tic;
            obj.F = F;
            obj.eqLHS = eqLHS;
            obj.eqRHS = eqRHS;
            obj.x0 = x0;
            obj.mode = mode;
            obj.use_weighted_metric = use_weighted_metric;
            obj.K = K;
            obj.lb = lb;
            obj.ub = ub;
            obj.dim = size(F,2)-1;
            [obj.T,areas] = computeMeshTranformationCoeffsMex(F, V); % compute mapping of vertices to differentials (T)
            weights = kron(sqrt(areas),ones(obj.dim^2,1));
            if use_weighted_metric,
                obj.W = sparse(1:length(weights),1:length(weights),weights);
            else
                obj.W = 1;
            end
            obj.initSolver(); % initialize
            obj.report(1,'SolverProjectorBD is ready (%.3g secs)\n', toc(t_start));
            obj.report(1,'-----------------------------------------------------------------\n\n');
        end
        
        function projectD_(obj)
            % project onto BD
            [obj.pTx, obj.distortions, obj.flips, obj.minsv, obj.maxsv] = projectBDMex(obj.Tx, obj.dim, obj.K, obj.lb, obj.ub);
        end
        
        function solve(obj, iter_max, tol_err)
            obj.report(1,'-----------------------------------------------------------------\n');
            obj.report(1,'BD PROJECTION (K=%g):\n', obj.K);
            obj.report(1,'-----------------------------------------------------------------\n');
            obj.report(1,'initial max dist %g,  flips %d  (infeasible %d)\n', max(obj.distortions), nnz(obj.flips), nnz((obj.distortions>obj.K)|obj.flips));
            obj.report(1,'(min sv %g,  max sv %d)\n', min(abs(obj.minsv)), max(obj.maxsv));
            obj.report(1,'-----------------------------------------------------------------\n');
            for iter = 1:iter_max
                obj.iterate();
                err = mse(obj.tanNormal);
                obj.report(1,'iter %d -- err: %g    dist: %g    flips: %g   time: %g sec (pD:%d%%, pL:%d%%)\n',iter,err,max(obj.distortions),nnz(obj.flips),obj.t_iter,round(100*obj.t_projectD/obj.t_iter), round(100*obj.t_projectLinear/obj.t_iter));
                if err<tol_err
                    obj.report(1,'err<tol_err --> stopping...\n');
                    break;
                end
            end
            obj.report(1,'-----------------------------------------------------------------\n');
            obj.report(1,'initial max dist %g,  flips %d  (infeasible %d)\n', max(obj.distortions), nnz(obj.flips), nnz((obj.distortions>obj.K)|obj.flips));
            obj.report(1,'-----------------------------------------------------------------\n\n');
        end
        
        function visualize(obj)
            [dist_colors,flip_color] = getColors;
            switch size(obj.F,2)
                case 3 % triangular mesh
                    sides = [1 2 3];
                case 4 % tetraherdral mesh
                    sides = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
                otherwise
                    error('Invalid problem type')
            end
            for kk = 1:size(sides,1)
                patch('vertices',obj.y,'faces',obj.F(:,sides(kk,:)),'FaceVertexCData',obj.distortions,'facecolor','flat','EdgeAlpha',0.1,'facealpha',0.6);
                patch('vertices',obj.y,'faces',obj.F(logical(obj.flips),sides(kk,:)),'FaceVertexCData',obj.distortions,'facecolor','none','edgecolor',flip_color);
            end
            axis equal;
            axis off;
            colormap(dist_colors);
            colorbar;
            caxis([1 1.5*obj.K]);
        end
    end
    
end

