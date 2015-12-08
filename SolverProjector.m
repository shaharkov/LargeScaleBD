classdef (Abstract) SolverProjector < handle
    
    properties
        % problem
        T; % lift operator
        W; % norm weights
        eqLHS;
        eqRHS;
        x0;
        % solver
        mode;
        nVars;
        nEq;
        x;
        Tx;
        pTx;
        tanNormal;
        tanLHS;
        tanRHS;
        preFactorization;
        TWW; % T'*W'*W
        TWWT; % T'*W'*W*T
        % display / log
        verbose = 1; % (0 = silent)
        t_iter;
        t_projectD;
        t_projectLinear;
        t_factorization;
    end
    
    properties (Dependent)
        y;
    end
    
    methods (Abstract)
        projectD_(obj)
    end
    
    methods
        function obj = SolverProjector
            % empty constructor
        end
        
        function value = get.y(obj)
            value = reshape(obj.x,size(obj.x0));
        end
        
        function initSolver(obj)
            obj.nVars = numel(obj.x0);
            obj.nEq = size(obj.eqLHS,1);
            obj.x = colStack(obj.x0);
            obj.Tx  = obj.T*obj.x;
            obj.projectD();
            obj.tanNormal = zeros(obj.nVars,1);
            obj.tanLHS = zeros(obj.nVars,1);
            obj.tanRHS = 0;
            obj.updateProblem();
        end
        
        function updateProblem(obj)
            t_start = tic;
            obj.TWW = obj.T'*obj.W'*obj.W;
            obj.TWWT = obj.TWW*obj.T;
            obj.factorize();
            obj.t_factorization = toc(t_start);
            obj.report(2,'Prectorization took (%.3g secs)\n', obj.t_factorization);
        end
        
        function factorize(obj)
            % construct KKT matrix
            LHS = [obj.TWWT, obj.eqLHS'; obj.eqLHS, sparse(obj.nEq,obj.nEq)];
            % factorize
            obj.preFactorization = SparseLU(LHS);
        end
        
        function projectD(obj)
            t_start = tic;
            obj.projectD_();
            obj.t_projectD = toc(t_start);
        end
        
        function projectLinear(obj)
            t_start = tic;
            switch obj.mode
                case SolverProjectorModeEnum.AltProj
                    temp_RHS = [obj.TWW*obj.pTx; obj.eqRHS];
                    temp_x_lambda = obj.preFactorization.solve(temp_RHS);
                    obj.x = temp_x_lambda(1:obj.nVars);
                case SolverProjectorModeEnum.Tangent
                    if any(obj.tanNormal) % avoid solving linear system if projecting on D didn't do anything
                        obj.tanLHS = obj.tanNormal'*obj.T;
                        obj.tanRHS = obj.tanNormal'*obj.pTx;
                        temp_rhs = [obj.TWW*obj.pTx; obj.eqRHS];
                        temp_Au = [obj.tanLHS'; zeros(obj.nEq,1)];
                        Fm_c = obj.preFactorization.solve(temp_rhs);
                        Fm_n = obj.preFactorization.solve(temp_Au);
                        temp_x_lambda = Fm_c - (temp_Au'*Fm_c - obj.tanRHS)/(temp_Au'*Fm_n)*Fm_n;
                    else % if tanNormal=0
                        temp_RHS = [obj.TWW*obj.pTx; obj.eqRHS];
                        temp_x_lambda = obj.preFactorization.solve(temp_RHS);
                    end
                    obj.x = temp_x_lambda(1:obj.nVars);
                otherwise
                    error('invalid mode');
            end
            obj.t_projectLinear = toc(t_start);
        end
        
        function iterate(obj)
            t_start = tic;
            obj.Tx = obj.T*obj.x; % lift
            obj.projectD(); % project onto D
            obj.tanNormal = obj.Tx - obj.pTx; % compute normal (=error)
            obj.projectLinear(); % project onto linear constraints
            obj.t_iter = toc(t_start);
        end
        
        function report(obj,verbosity,varargin)
            if verbosity<=obj.verbose
                fprintf(varargin{:});
            end
        end
    end
end

