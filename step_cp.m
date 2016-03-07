function [ xn ] = step_cp(x,u,params)
%returns next state vector starting from state x,applying control u, and
%doing euler integration with step size params.dt
%x= [x,theta,xdot,thetadot]
dt = params.dt;
m1 = params.m1;
m2 = params.m2;
M = m1+m2;
l = params.l;
g = params.g;
mu = params.mu; %viscous friction coefficient

x4dot = (u(1)+m2*l*x(4)*x(4)*sin(x(2))+M*g*tan(x(2)) - M*u(2)/cos(x(2)))/(m2*l*cos(x(2))-(M*l+M*mu)/cos(x(2)));
x3dot = (-g*sin(x(2))-(l+mu)*x4dot+u(2))/cos(x(2));
xn = zeros(4,1);
xn(1) = x(1)+dt*x(3);
xn(2) = x(2) + dt*x(4);
xn(3) = x(3) + dt*x3dot;
xn(4) = x(4) + dt*x4dot;

xn = xn';


end

