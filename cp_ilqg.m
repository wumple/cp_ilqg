clear all;

params.dt = 0.005; %euler time step length (s)
params.m1 = 0.5; %mass of cart (kg)
params.m2 = 20.5; %mass on top of pole (kg)
params.l = 1; %length of pole (m)
params.g = 9.81; %gravity, (m/s/s)
params.mu = 0.6; %viscous damping coefficient
params.nsteps = 400; %number of states
params.T = params.nsteps-1; %number of state transitions
params.Tcost = 10000; %control cost weight on the pure torque at the pole joint
params.Fcost = 1e-4; %control cost weight on the force applied to the cart
params.alpha = 0.9; %multiply step size by this amount on each iteration of
                    %the line search
params.xcost = 100;%cost weight on the distance from upright

%see paper for definitions, "A generalized iterative LQG method for locally
%-optimal feedback control of constrained nonlinear stochastic systems" 
%Todorov/Li
A = cell(params.nsteps,1);
B = cell(params.nsteps,1);
q = cell(params.nsteps,1);
q_twiddle = cell(params.nsteps,1);
Q = cell(params.nsteps,1);
P = cell(params.nsteps,1);
r = cell(params.nsteps,1);
R = cell(params.nsteps,1);
g = cell(params.nsteps,1);
G = cell(params.nsteps,1);
H = cell(params.nsteps,1);
I = cell(params.nsteps,1);
L = cell(params.nsteps,1);
s = cell(params.nsteps,1);
s_twiddle = cell(params.nsteps,1);
S = cell(params.nsteps,1);

x0 = [0;0;0;0]; %[x (pos of cart);theta;x dot; theta dot]
us = (rand(2,(params.nsteps-1))-0.5)*50; %[force on cart;torque on pole joint]
niter = 10; %number of optimization iterations
losses = zeros(niter,1);
u_bar = us; %linearization point
backtrack_lim = 50; %maximum number of backtracking line search attempts

for j = 1:niter
    fprintf('%d/%d\n',j,niter);
    [~,x_bar] = sim_loss( x0,u_bar,params); %get the linearization points
    %x_bar and u_bar
    
    
    %first pass is a forward pass
    %here we compute all of the required matrices for the 
    %linearized dynamics and quadraticized costs

    for i = 1:params.nsteps
        if i == params.nsteps
            q{i} = params.dt*loss_cp(x_bar(:,i),zeros(2,1),params);
            q_twiddle{i} = dx_l(x_bar(:,i),zeros(2,1),params);
            Q{i} = dxx_l(x_bar(:,i),zeros(2,1),params);
            A{i} = zeros(4,4);
            B{i} = zeros(4,4);
            P{i} = zeros(2,4);
            r{i} = zeros(2,1);
            R{i} = zeros(2,2);
        else
            A{i} = eye(4,4) + params.dt*jac_ilqg_f(x_bar(:,i),u_bar(:,i),params)';
            B{i} = params.dt*ilqg_du_f(x_bar(:,i),u_bar(:,i),params)';
            P{i} = zeros(2,4); 
            r{i} = params.dt*du_l(x_bar(:,i),u_bar(:,i),params);
            R{i} = params.dt*duu_l(x_bar(:,i),u_bar(:,i),params);
            q{i} = params.dt*loss_cp(x_bar(:,i),u_bar(:,i),params);
            q_twiddle{i} = params.dt*dx_l(x_bar(:,i),u_bar(:,i),params);
            Q{i} = params.dt*dxx_l(x_bar(:,i),u_bar(:,i),params);
        end
    end

    %second pass is backwards, and we do the actual optimization
    %of the controls here
   
    for i = 1:params.nsteps
        k = (params.nsteps-(i-1));
        if k == params.nsteps
            S{k} = Q{k};
            s_twiddle{k} = q_twiddle{k};
            s{k} = q{k};
        else
            g{k} = r{k} + B{k}'*s_twiddle{k+1};
            G{k} = P{k} + B{k}'*S{k+1}*A{k};
            H{k} = R{k} + B{k}'*S{k+1}*B{k};
            [V,D] = eig(H{k});
            
            if sum(D < 0) > 0
                disp('got negative eigenvalues!');
            end
            %make the eigenvalues of the hessian positive...
            %they might not be due to our approximations
            min_eig = min(min(diag(D)),0);
            eps = 1e-6;
            I{k} = -V*(inv(D+eye(2,2)*(eps-min_eig)))*V'*g{k};
            L{k} = -V*(inv(D+eye(2,2)*(eps-min_eig)))*V'*G{k};
            s{k} = q{k} + s{k+1} + 0.5*I{k}'*H{k}*I{k}+I{k}'*g{k};
            s_twiddle{k} = q_twiddle{k} + A{k}'*s_twiddle{k+1}+L{k}'*H{k}*I{k}+L{k}'*g{k}+G{k}'*I{k};
            S{k} = Q{k} + A{k}'*S{k+1}*A{k} + L{k}'*H{k}*L{k}+L{k}'*G{k}+G{k}'*L{k};
        end
    end
    
    [dus,dxs] = compute_dus(A,B,I,L,params); %make a third pass to calculate 
    %the open loop and the feedback terms 
    
    us_twiddle = u_bar + dus; %the updated control term, but do 
    % a backtracking linesearch before actually committing to this set of
    % controls...that happens next
    
    %backtracking line search:
    if j > 1
        count = 1;
        while(sim_loss(x0,us_twiddle,params) > losses(j-1))
            
            for p = 1:params.nsteps
                I{p} = params.alpha*I{p};
            end
            
            [dus,dxs] = compute_dus(A,B,I,L,params);
            us_twiddle = u_bar + dus;
            count = count+1;
            
            if count > backtrack_lim
                fprintf('Took more than %d backtracking steps\n',backtrack_lim);
                for p = 1:params.nsteps
                    I{p} = 0;
                end
                break;
            end
            
        end
    end
    
    %all done!  update controls and repeat...
    u_bar = us_twiddle;
    losses(j) = sim_loss(x0,u_bar,params);
end


plot(losses);
xlabel('Iteration');
ylabel('Trajectory cost');
u_bar(2,:) = 0;
[loss,xs] = sim_loss(x0,u_bar,params);
times = (1:1:params.nsteps)*params.dt;
figure;
plot(times,xs(2,:));
xlabel('Time (s)');
ylabel('Theta (rad)');
    
