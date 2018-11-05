classdef pointSolver < handle
    % point Solver determines point sources from Cauchy data of Helmholtz
    % equation.
    
    properties (Access = public)
        refcIdx
        points
        strengths
        
        model
        cache
    end
   
    
    methods
        function obj = pointSolver(opt)
            % The domain is unit square with fine mesh.
            if nargin < 1 
                opt = struct('deg', 4, 'qdeg', 8, 'min_area', 1e-4,...
                      'edge', [0 1 1 0; 0 0 1 1], 'hull',...
                      [0 0 1 0 1 1 0 1]',...
                'idx', [0 1 1 2 2 3 3 0]');
            
                opt.pointStyle  = 'random';
                opt.pointDist   = 0.2;
                opt.minStrength = 0.5;
                opt.refcStyle   = 'random';
                opt.waveNum     = 8;
            end

            % generate mesh
            cprintf('cyan',                 '1. Generating triangulation mesh   ...\t');
            tic;   
            obj.model = femm(opt);  
            
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
                nPoints = randi(5);
                obj.points = rand(2, nPoints);
                searchCnt = 0;
                while obj.minDist(obj.points) < opt.pointDist && searchCnt < 100
                   obj.points = rand(2, nPoints);
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
            
            obj.cache.load  = obj.model.build('l', gaussian_source); % init
            
        end
        
        function [f] = forwardSolve(obj, Idx)
            % given the refc Idx and zero Neumann, return the Neumann data on full boundary.
            % partial boundary is possible for later.
            
            % projection by interpolation, find quadrature point values.
            
            qIdx = obj.mapping(Idx, obj.model.space.elems, obj.model.facet.ref');
            M    = obj.model.build('m', qIdx); % mass matrix
            S    = obj.model.build('s', 1);    % stiffness matrix
            
            A    = obj.cache.k^2 * M - S;
            % solve Neumann problem.
            
            u = A \obj.cache.load; 
            
            % return Dirichlet data
            f = u(obj.cache.ndof);
            
        end
        
        function visualize(obj, u)
            trisurf(obj.model.space.elems(1:3,:)', ...
                obj.model.space.nodes(1,:), obj.model.space.nodes(2,:), u,...
                'EdgeColor', 'None');
            colormap jet;view(2);colorbar;shading interp;
            hold on;
            scatter(obj.points(1,:), obj.points(2,:));
            hold off;
        end

        function [p, s] = backwardSolve(obj)
            % todo
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
        
    end
    
end

