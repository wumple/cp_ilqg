function [ jac ] = jac_ilqg_f( x,u,params )
%analytic jacobian of f as defined in ilqg paper (full title in main script) at point
%x,u with system parameters 'params'

dt = params.dt;
m1 = params.m1;
m2 = params.m2;
M = m1+m2;
l = params.l;
g = params.g;
mu = params.mu;

x4dot = (u(1)+m2*l*x(4)*x(4)*sin(x(2))+M*g*tan(x(2)) - M*u(2)/cos(x(2)))/(m2*l*cos(x(2))-(M*l+M*mu)/cos(x(2)));

jac = zeros(4,4);

%derivs wrt x(3) (first element of f)
jac(1,1) = 0;
jac(2,1) = 0;
jac(3,1) = 1;
jac(4,1) = 0;

%derivs wrt x(4) (second element of f)
jac(1,2) = 0;
jac(2,2) = 0;
jac(3,2) = 0;
jac(4,2) = 1;

%derivs wrt fourth output of f
jac(1,4) = 0;
den = m2*l*cos(x(2)) - M*(l+mu)/cos(x(2));
t1 = -(-M*(l+mu)*(1/(cos(x(2))*cos(x(2))))*sin(x(2))-m2*l*sin(x(2)))*(u(1)+m2*l*x(4)*x(4)*sin(x(2))+M*g*tan(x(2))-M*u(2)/cos(x(2)))/(den*den);
t2 = (m2*l*x(4)*x(4)*cos(x(2))+M*g*sec(x(2))*sec(x(2))-M*u(2)*(1/(cos(x(2))*cos(x(2))))*sin(x(2)))/den;
jac(2,4) = (t1+t2);
jac(3,4) = 0;
jac(4,4) = (2*m2*l*x(4)*sin(x(2)))/den;

%derivs wrt third output of f
jac(1,3) = 0;
jac(2,3) = (-g*sec(x(2))*sec(x(2)) - ((tan(x(2))/cos(x(2)))*(l+mu)*x4dot + (l+mu)*(jac(2,4))*(1/cos(x(2))))+u(2)*tan(x(2))/cos(x(2)));
jac(3,3) = 0;
jac(4,3) = (-(l+mu)*((jac(4,4))))/cos(x(2));


end

