function [ Z ] = landscape(obj, rps)
%LANDSCAPE Summary of this function goes here
%   Detailed explanation goes here

% We pick the x-coordinates of first two points in rps and formulate  the
% landscape.
N = 42;
V = linspace(0,1,N);
[X,Y]= meshgrid(V(2:end-1));
X_ =  X(:) ;
Y_ =  Y(:) ;

Z = zeros(length(X_), 1);
tic;
for i = 1:length(Z)
    rps_ = rps;
    if mod(i, 20) == 0
        toc;
        disp(i);
        tic;
    end
    rps_(1:2) = [X_(i);Y_(i)];

    Z(i) = obj.backwardSolve(rps_);
end
toc;

Z = reshape(Z, N-2, N-2);

surf(X,Y,imgaussian(Z,1) / 80, 'EdgeColor', 'None');shading interp; colormap(flipud(jet));view(2);colorbar;

end

