classdef pointSolver < handle
    % point Solver determines point sources from Cauchy data of Helmholtz
    % equation.
    
    properties (Access = public)
        refcIdx
        points
        strengths
        
        model
        cache
        
        measure
        
        result
    end
   
    
    methods
        function obj = pointSolver(opt)
            % The domain is unit square with fine mesh.
            if nargin < 1 
                opt = struct('deg', 2, 'qdeg', 4, 'min_area', 5e-5,...
                      'edge', [0 1 1 0; 0 0 1 1], 'hull',...
                      [0 0 1 0 1 1 0 1]',...
                'idx', [0 1 1 2 2 3 3 0]');
            
                opt.pointStyle   = 'random';
                opt.pointDist    = 0.2;
                opt.minStrength  = 0.5;
 
                opt.waveNum      = 5 ; %4 * sqrt(-1);
                opt.refcIdxStyle = 'random';
            end

            % generate mesh
            cprintf('cyan',                 '1. Generating triangulation mesh   ...\t');
            tic;   
            obj.model = femm(opt);  
            obj.measure = {};
            obj.result  = {};
            obj.result.convergence = [];
            % generate boundary information.
            obj.cache.n = size(obj.model.space.nodes, 2);
            obj.cache.ndof = unique(obj.model.space.edges);
            obj.cache.dof = setdiff(1:obj.cache.n, obj.cache.ndof);
            obj.cache.k   = opt.waveNum;
            t = toc; 
            cprintf('cyan',   sprintf('%f seconds  \n', t));
            
            % generate point source
            tic;
            if strcmp(opt.pointStyle, 'random')
                nPoints = 4;
                obj.points = 0.5 + 0.4 * (2 * rand(2, nPoints) - 1);
                searchCnt = 0;
                while obj.minDist(obj.points) < opt.pointDist && searchCnt < 100
%                    nPoints = randi(10);
                   obj.points = 0.5 + 0.25 * (2 * rand(2, nPoints) - 1);
                   searchCnt = searchCnt + 1;
                end
                if searchCnt == 100
                    cprintf('cyan',         '2. Cannot find a configuration     ...\t');
                else
                    cprintf('cyan', sprintf('2. Find a configuration at %d step  ...\t', searchCnt));
                end
                obj.strengths = opt.minStrength + rand(1, size(obj.points, 2));
            else
                obj.points = opt.points;
                obj.strengths = opt.strengths;
                cprintf('cyan',             '2. Points loaded from option       ...\t');
            end
            t = toc;
            cprintf('cyan',   sprintf('%f seconds  \n', t));   
            
            % generate load vector.
            gaussian_source = zeros(size(obj.model.quad2d, 2), 1);
            for i = 1:size(obj.points, 2)
                gaussian_source = gaussian_source +...
                    obj.gaussian(obj.model.quad2d, obj.points(:, i), obj.strengths(i));
            end
           
            
            if strcmp(opt.refcIdxStyle, 'random')
                obj.refcIdx = 1.0 + 1.0 * rand(size(obj.model.space.nodes, 2), 1);
            else
                if ismatrix(opt.refcIdx)
                    obj.refcIdx = opt.refcIdx;
                else
                    obj.refcIdx = opt.refcIdxFunc(obj.model.space.nodes);
            
                end
            end
            [obj.cache.u_ground, obj.measure.dirichlet] = ...
                obj.forwardSolve(obj.refcIdx, gaussian_source);
            figure(1);
            title('solution');
            obj.visualize(obj.cache.u_ground);
        end
        
        function [u, f] = forwardSolve(obj, Idx, gaussian_source)
            % given the refc Idx and zero Neumann, return the Neumann data on full boundary.
            % partial boundary is possible for later.
            
            % projection by interpolation, find quadrature point values.
            
            qIdx = obj.mapping(Idx, obj.model.space.elems, obj.model.facet.ref');
            M    = obj.model.build('m', qIdx); % mass matrix
            S    = obj.model.build('s', 1);    % stiffness matrix
            
            obj.cache.A    = obj.cache.k^2 * M - S;
            % solve Neumann problem.
            load  = obj.model.build('l', gaussian_source);
            u = obj.cache.A \load; 
            
            % return Dirichlet data
            f = u(obj.cache.ndof);
            
        end
        
        function visualize(obj, u)
            trisurf(obj.model.space.elems(1:3,:)', ...
                obj.model.space.nodes(1,:), obj.model.space.nodes(2,:), u,...
                'EdgeColor', 'None');
            colormap jet;view(2);colorbar;shading interp;
            hold on;
            s = scatter3(obj.points(1,:), obj.points(2,:), 100 * ones(size(obj.strengths)), 120, 'o');
            s.LineWidth = 2;
            s.MarkerEdgeColor='k';
            s.MarkerFaceColor=[0 .75 .75];
            view(2);
            hold off;
        end

        function [loss, grad] = backwardSolve(obj, rps)
            % generate loss and gradient
            rps_stack = reshape(rps, 3, numel(rps)/3);

            gaussian_source = zeros(size(obj.model.quad2d, 2), 1);
            for i = 1:size(obj.points, 2)
                gaussian_source = gaussian_source +...
                    obj.gaussian(obj.model.quad2d, rps_stack(1:2, i), rps_stack(3, i));
            end
                       

            
            figure(2);
            title('reconstruction');
            s1 = scatter3(rps_stack(1, :), rps_stack(2,:), rps_stack(3,:), 120, 'o');
            s1.LineWidth = 2;
            xlim([-0.2 1.2]);
            ylim([-0.2 1.2]);
            hold on;
            s2 = scatter3(obj.points(1,:), obj.points(2,:), obj.strengths, 120, 'x', 'LineWidth', 1);
            s2.LineWidth = 2;
            hold off;
            colormap jet;view(2);
            drawnow();
            
            [cur_u, curDirichlet] = forwardSolve(obj, obj.refcIdx, gaussian_source);
            mismatch = (curDirichlet - obj.measure.dirichlet);
            
                        
            figure(3);
            title('difference')
            obj.visualize(log10(abs(cur_u - obj.cache.u_ground)));
            
            loss = 0.5 * (mismatch' * mismatch);
            % regularization to penalty the gradient if point is outside?
 
            obj.result.convergence = [ obj.result.convergence , loss];
            
            
           
            
            % gradient is in group of 3: x, y, s.
            
            % augment
            augment_mismatch = zeros(size(obj.model.space.nodes, 2), 1);
            augment_mismatch(obj.cache.ndof) = mismatch;
            feedback = obj.cache.A \ augment_mismatch;
            
            p_grad = zeros(length(rps), size(obj.model.space.nodes, 2));
            for p_id = 1:length(obj.strengths)
                % x 
                px = obj.px_gaussian(obj.model.quad2d, ...
                    rps((3 * p_id-2):(3*p_id-1)), rps(3 * p_id));
                p_grad(3 * p_id - 2, :) = obj.model.build('l', px);
                
                p_grad(3 * p_id - 2, :) = p_grad(3 * p_id - 2, :) ;
                
                
                py = obj.py_gaussian(obj.model.quad2d, ...
                    rps((3 * p_id-2):(3*p_id-1)), rps(3 * p_id));
                p_grad(3 * p_id - 1, :) =  obj.model.build('l', py) ;
                
                ps = obj.ps_gaussian(obj.model.quad2d, ...
                    rps((3 * p_id-2):(3*p_id-1)));
                p_grad(3 * p_id - 0, :) = obj.model.build('l', ps);
            end
            
            grad = p_grad * feedback;
%             
        end
        
        
        function [rps, hist] = reconstruction_lbfgs(obj, rps)

            if nargin < 2
                pre_rps = reshape([obj.points; obj.strengths], 3 * length(obj.strengths), 1) ;
                rps = pre_rps + 0.1 * (2 * rand(size(pre_rps)) - 1);
%                   rps = rand(30, 1);
            end
            
            options    = struct( 'factr', 1e0, 'pgtol', 1e-10, 'm', 20, ...
                'x0', rps, 'maxIts', 200, 'maxTotalIts', 1e5);
            
            options.printEvery     = 1;            
            
            [rps, ~, hist] =...
                lbfgsb_c(@obj.backwardSolve, -inf * ones(size(rps)), inf * ones(size(rps)), options);  
        end
        
        function [rps,fval,exitflag,output] = reconstruction(obj, rps)
            % reconstruction for point sources.
            % use BFGS for quasi-Newton's method?
            options = optimset('Diagnostics','on','DerivativeCheck','off',...
                'FinDiffType','central','LargeScale','off',...
                'GradObj','on','Display','iter-detailed',...
                'TolFun',1e-12,'TolX',1e-9,'MaxFunEvals',4000,....
                'MaxIter',500,'HessUpdate','bfgs');
            if nargin < 2
                pre_rps = reshape([obj.points; obj.strengths], 3 * length(obj.strengths), 1) ;
                rps = pre_rps + 0.1 * (2 * rand(size(pre_rps)) - 1);
%                   rps = rand(30, 1);
            end
            [rps,fval,exitflag,output] = fminunc(@obj.backwardSolve, rps,options);
            disp(reshape(rps, 3, numel(rps)/3));
            disp([obj.points; obj.strengths]);
        end
        
    end
    
    methods (Static)
        function d = minDist(points)
            % Brute-force search for minimum distance between points.
            d = inf;
            n = size(points, 2);
            for i = 1:(n-1)
                for j = (i+1):n
                    if norm(points(:, i) - points(:, j)) < d
                        d = norm(points(:, i) - points(:, j));
                    end
                end
            end
        end
        
        function [interpolate] = mapping(func, elems, trans_ref)
            numberofqnodes = size(trans_ref, 1);
            interpolate = zeros(numberofqnodes, size(elems, 2));
            for i = 1: size(elems, 2)
                interpolate(:, i) = trans_ref * func(elems(:, i));
            end
        end        
        
        function [source] = gaussian(nodes, point, strength)
            diff   = bsxfun(@minus, nodes, point);
            dist   = sum((diff').^2, 2);
            sigma  = 0.01;
            source = strength * exp(-dist / sigma^2) / (sigma^2);
        end
        
        function [source] = px_gaussian(nodes, point, strength)
            diff   = bsxfun(@minus, nodes, point);
            dist   = sum((diff').^2, 2);
            sigma  = 0.01;
            source = strength * exp(-dist / sigma^2) / (sigma^2) .*...
                (2 * diff(1,:)' / sigma^2);
        end
        
                
        function [source] = py_gaussian(nodes, point, strength)
            diff   = bsxfun(@minus, nodes, point);
            dist   = sum((diff').^2, 2);
            sigma  = 0.01;
            source = strength * exp(-dist / sigma^2) / (sigma^2) .*...
                (2 * diff(2,:)' / sigma^2);
        end
        
               
        function [source] = ps_gaussian(nodes, point)
            diff   = bsxfun(@minus, nodes, point);
            dist   = sum((diff').^2, 2);
            sigma  = 0.01;
            source = exp(-dist / sigma^2) / (sigma^2);
        end
        
        
    end
    
end

