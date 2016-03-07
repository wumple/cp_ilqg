function [ h ] = duu_l( x,u,params)
%hessian of loss wrt inputs u for cp problem at state x with parameters
%params

h = zeros(2,2);
h(1,1) = params.Fcost*2;
h(2,2)= params.Tcost*2;


end

