 function rkc_params = RKC_parameters(order,dt,sprad,s)

if (nargin < 4)
    if (order == 1)
        s = 1 + ceil(sqrt(dt*sprad/0.19)); %s = 1+ceil(sqrt(dt*4/(dx^2*0.653)));
    else
        s = 1 + ceil(sqrt(dt*sprad/0.653)); %s = 1+ceil(sqrt(dt*4/(dx^2*0.653)));
    end
end

rkc_params.s = s;
rkc_params.order = order;

% --------------------------------------------
% Construct Chebyshev polynomials recursively
% --------------------------------------------

if order == 1
    ep = 0.05;
else
    ep = 2/13;
end

% --------------------------------------------
% Construct Chebyshev polynomials recursively
% --------------------------------------------

s_idx = s+1; % for indexing matlab arrays
idx_0 = 1;
idx_1 = 2;
idx_2 = 3;
w0 = 1 + ep/s^2;


T(idx_0) = 1; 
T(idx_1) = w0;
T(idx_2) = 2*w0^2-1;

dT(idx_0) = 0;
dT(idx_1) = 1;
dT(idx_2) = 4*w0;

dT2(idx_0) = 0;
dT2(idx_1) = 0;
dT2(idx_2) = 4;

for j = 3:s
    j_idx = j+1;
    T(j_idx) = 2*w0*T(j_idx-1)-T(j_idx-2);
    dT(j_idx) = 2*w0*dT(j_idx-1)-dT(j_idx-2) + 2*T(j_idx-1);
    dT2(j_idx) = 2*w0*dT2(j_idx-1)-dT2(j_idx-2) + 4*dT(j_idx-1); 
end

% --------------------------------------
% Get parameters for each order
% --------------------------------------
if order == 1
    w1 = T(s_idx)/dT(s_idx);
    
    b = NaN(1,s+1);
    for j = 0:s,
        j_idx = j+1;
        b(j_idx) = 1./T(j_idx);
    end
    Mu(idx_0) = NaN;    
    MuT(idx_0) = NaN;   

    Mu(idx_1) = b(idx_1)*w1;
    MuT(idx_1) = w1/w0;
    
    for j = 2:s
        j_idx = j+1;
        Mu(j_idx) = 2*b(j_idx)*w0/b(j_idx-1);
        Nu(j_idx) = -b(j_idx)/b(j_idx-2);
        MuT(j_idx)= 2*b(j_idx)*w1/b(j_idx-1);       
        GammaT(j_idx) = 0;
    end

    c = NaN(1,s+1);
    for j = 1:s
        j_idx = j+1;
        c(j_idx) = T(s_idx)/dT(s_idx)*dT(j_idx)/T(j_idx);
    end
    c(idx_0) = 0;
    c(s_idx) = 1;   % already 1, but still might be a good idea to set to 1
else

    arg = s*log(w0+sqrt(w0^2-1));
    w1 = sinh(arg)*(w0^2-1)/(cosh(arg)*s*sqrt(w0^2-1)-w0*sinh(arg));
    % w1 = dT(s_idx)/dT2(s_idx);

    b = NaN(1,s+1);
    for j = 2:s
        j_idx = j+1;
        b(j_idx) = dT2(j_idx)/dT(j_idx)^2;
    end
    b(idx_0) = b(idx_2);
    b(idx_1) = b(idx_2);

    Mu(idx_0) = NaN;   
    MuT(idx_0) = NaN;   
    MuT(idx_1) = b(idx_1)*w1;
    for j = 2:s
        j_idx = j+1;
        Mu(j_idx) = 2*b(j_idx)*w0/b(j_idx-1);
        Nu(j_idx) = -b(j_idx)/b(j_idx-2);
        MuT(j_idx)= 2*b(j_idx)*w1/b(j_idx-1);       
        GammaT(j_idx) = -(1-b(j_idx-1)*T(j_idx-1))*MuT(j_idx);
    end

    c = NaN(1,s+1);
    for j = 2:s
        j_idx = j+1;
        c(j_idx) = dT(s_idx)/dT2(s_idx)*dT2(j_idx)/dT(j_idx);
    end
    c(idx_1) = c(idx_2)/dT(idx_2);   % needed here, but not for RKC1
    c(idx_0) = 0;
    c(s_idx) = 1;
end


% Store parameters
rkc_params.Mu = Mu;
rkc_params.Nu = Nu;
rkc_params.MuT = MuT;
rkc_params.GammaT = GammaT;
rkc_params.c = c;


end






















