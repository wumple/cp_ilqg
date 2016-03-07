function [ dus,dxs ] = compute_dus(As,Bs,Is,Ls,params)
%computes the updated controls for ilqg
%As are linearized dynamics matrices,
%Bs are linearized control matrices,
%Is are open loop controls
%Ls are closed loop gains
%params is set of system parameters

    dxs = zeros(4,params.nsteps);
    dus = zeros(2,params.T);
    for i=1:params.nsteps-1
        dus(:,i) = Is{i} + Ls{i}*dxs(:,i);
        dxs(:,i+1) = As{i}*dxs(:,i) + Bs{i}*dus(:,i);
    end


end

