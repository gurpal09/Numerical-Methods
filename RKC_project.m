% Gurpal Singh
% Math 567
% 12/13/2017

% About:
%   This project uses 3 different methods to approximate the solution
%   to diffusion in a one-dimensional rod. The methods used are:
%       i.   Explicit Method:           Forward Euler + Central Difference
%       ii.  Implicit/Direct Method:    Backward Euler + Central Difference
%       iii. Runge-Kutta + Chebychev:   Central Difference + RKC
method = 'rkc';

% RKC Settings
order = 2;

plot_soln = false;
T = .4; % Total Time
% m = 1;
x_start = 0; % Starting x
x_end = 1;   % Ending x


% Eigenmode solution
eigmode = @(x,t,m) exp(-m^2*pi^2*t).*sin(m*pi*x);

m = 1;
% Exact Solution
uexact = @(x,t) eigmode(x,t,m);


uxx = @(x,t) -(pi^2)*(m^2).*sin(m.*pi.*x)*exp(-(m^2)*(pi^2)*t);

% Set Spatial Scale and Convergence Study values
N0 = 8;
I = 0:5;
Nvec = N0*2.^I;

% Print to screen
l = double('-')*ones(1,40);
fprintf('%s\n',l);
fprintf('%6s %6s %8s %10s %6s\n','N','M','k/h','err','rate');
fprintf('%s\n',l);

% Plot Title Settings
if plot_soln
    close all
    switch method
        case 'fe'
            tstr = 'Forward Euler';
        case 'be'
            tstr = 'Backward Euler';
        case 'rkc'
            tstr = 'RKC Method';
    end
end


for i = 1:length(I) % Convergence Study Loop
    N = Nvec(i);    % Set Nvalue
%     N = 10;
    h = (x_end-x_start)/N; % Define h based off of N

    xe = linspace(x_start,x_end,N+1)'; % x coordinates
    x = xe(2:end-1);                   % Interior Coordinates
        
    % Build A matrix with Dirichlet conditions
    z = ones(N-1,1);
    A = spdiags([z -2*z z],-1:1,N-1,N-1);
    
    switch method % Pick which Method
        case 'fe'
            k_est = 0.45*h^2;
            M = round(T/k_est) + 1;
            k = T/M;
            B = eye(N-1) + (k/h^2)*A;
        case 'be'
            k_est = h;
            M = round(T/k_est) + 1;
            k = T/M;
            B_lhs = eye(N-1) - (k/h^2)*A;
        case 'rkc'
            
            % Set k to 0.5h
            k_est = 0.5*h;
            M = round(T/k_est) + 1;
            k = T/M;
           
            
            % Calculate Spectral Radius of A
            A = (1/(h^2)) .*A;
            eigvals = eig(A);
            spec_rad = max(abs(eigvals));
            
            % Calculate RKC Parameters
            params = rkc_params(order, k, spec_rad);
    end
    
    t = zeros(1,M+1); % Time Vector
   
    un = uexact(x,0);    % Interior values only
    
    % Prepare Figure for Plot
    if (plot_soln)
        clf
        clear a b
        ue = uexact(xe,0);
        b = plot(xe,ue,'black','linewidth',2);
        hold on;
        a = plot(xe,[0;un;0],'ro','linewidth',2,'markersize',10);
        lh = legend([a,b],{'Approximate solution','Exact solution'});
        title(tstr,'fontsize',18);
        axis([0 1 -1 1]);
        set(gca,'fontsize',16);
    end 
    
   
    for n = 1:M % Time Loop
        t(n+1) = n*k; % Set Current Time
        tn = t(n+1);
        switch method % Select Method
            case 'fe'
                unp1 = B*un;
            case 'be'
                unp1 = B_lhs\un;
            case 'rkc'
                
                % Set wn to initial condition for 1st timestep
                if (n == 1)
                wn = uexact(x,0);
                end
                
                % Calculate W0 and W1
                W0 = wn;
                MuT_1 = params.MuT(2);
                F0 = A*W0;
                W1 = W0 + MuT_1 * k * F0;

                W_jm1 = W1;
                W_jm2 = W0;
           
                for j = 2:params.s
                    % Set Coefficients
                    Mu_j = params.Mu(j+1);
                    Nu_j = params.Nu(j+1);
                    MuT_j = params.MuT(j+1);
                    GammaT_j = params.GammaT(j+1);
                    c_j = params.c(j+1);
                
                    % Calculate F_jm1
                    F_jm1 = A*W_jm1;
                    
                    % Calculate Current Stage
                    Wj = (1 - Mu_j - Nu_j)*W0 + Mu_j*W_jm1 + Nu_j*W_jm2 ...
                        + MuT_j*k*F_jm1 + GammaT_j*k*F0;
                    
                    % Update
                    W_jm2 = W_jm1;
                    W_jm1 = Wj;
                end
                
                
                % Save W_s
                unp1 = Wj;
                wn = Wj;
        end
        
        ue = uexact(xe,t(n+1)); % Compute Exact Solution at that timestep
        
        % Plotting
        U = [0; unp1; 0]; % Add on Boundary Conditions
            if (plot_soln)
                set(a,'ydata',U);
                set(b,'ydata',ue);
            end
        pause(0.6)
        un = unp1;
    end
    
    % Compute error at end of time step
    err(i) = norm(ue(2:end-1)-unp1,Inf);
    if (i == 1)
        rate = '---';
    else
        r = log(err(i-1)/err(i))/log(2);
        rate = sprintf('%6.2f',r);
    end
    fprintf('%6d %6d %7.4f %10.2e %6s\n',N,M,k/h,err(i),rate); 
end
    