function [ d ] = dx_l( x,u,params )
%derivative of cart pole loss wrt x
%returns a 1x4 vector
%u is applied controls
%params is structure of system parameters
%x is current state

d = zeros(1,4);
d(1) = 0;
d(2) = 2*params.xcost*(x(2)-pi);
d(3) = 0;
d(4) = 0;
d = d';


end

