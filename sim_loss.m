function [ loss,xs] = sim_loss( x0,us,params)
%loss of entire trajectory of cart pole system started
%at x0 using open loop controls 'us' with system parameters 'params'

    nsteps = params.nsteps;
    loss = 0;
    xs = zeros(4,nsteps);
    xs(:,1) = x0;
    for i=1:nsteps-1
        loss = loss + loss_cp(xs(:,i),us(:,i),params);
        xs(:,i+1) = step_cp(xs(:,i),us(:,i),params);
    end

    loss = loss + loss_cp(xs(:,nsteps),zeros(2,1),params);

end

