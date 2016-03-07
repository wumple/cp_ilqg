function [ d ] = du_l( x,u,params )
%derivative of loss for cart pole problem wrt to u evaluated
%at x,u with parameters params

d(1) = 2*params.Fcost*u(1);%2*u*0.001;
d(2) = 2*params.Tcost*u(2);
d=d';

end

