function [h ] = dxx_l( x,u,params )
%hessian of loss wrt x for cart pole problem for a single state/control
%x is state
%u is control
%params is structure of parameters

h = zeros(4,4);
h(2,2) = 2*params.xcost;


end

