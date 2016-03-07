function [ l ] = loss_cp(x,u,params)

%cart pole loss for driving system to unstable equilibrium evaluated
%at state x applying control u under parameters 'params'
%u should be a 2x1 vector, with u(1) = force applied to cart
%u(2) = torque applied to pole joint

l = params.xcost*(x(2)-pi)*(x(2)-pi) + params.Tcost*u(2)*u(2) + params.Fcost*u(1)*u(1);

end

