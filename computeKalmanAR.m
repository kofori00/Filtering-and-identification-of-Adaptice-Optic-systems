function [A,Cw,K]=computeKalmanAR(C_phi, C_phi1, G, sigma_e)
% Calculation of Kalman gain , A matrix and process noise covariance Cw
% INPUT:
%   C_phi - covariance matrix of phi(k)*phi(k)^T
%   C_phi1 - covariance matrix of phi(k)*phi(k+1)^T
%   sigma_e - variance of measurement noise
% OUTPUT:
%   A - system matrix in model phi(k+1) = Aphi*(k) + w(k)
%   Cw - process noise covariance
%   K - Kalman gain
% 


A=C_phi1/C_phi; 
C=G;
Cw=round(C_phi-A*C_phi*A', 8);
S=0;
R=sigma_e^2 * eye(size(G,1));


% Solve the Ricatti equation
[P,L,K]=dare(A',C',Cw, R);
K = K';
end

\end{lstlisting}

\begin{lstlisting}[frame=single]
function [var_eps] = AOloopAR(G, H, C_phi, C_phi1, sigma_e, phi)
% Closed loop simulation of AO system based on Kalman gain
% INPUT:
%   G - matrix relation between phi and slopes measurement s
%   H - matrix relation between input u(k-1) and deformable mirror induced
%       wavefront phi_dm(k)
%   C_phi, C_phi1 - covariace matrices of wavefront phi
%   sigma_e - measurement noise variance
%   phi - incoming turbulent wavefront used for simulation only
% OUTPUT:
%   var_eps - residual wavefron variance

%% Define necessary matrices
    
    W_e = 1/sigma_e^2;
    W = inv(sqrtm(C_phi)*sqrtm(C_phi)');
%     A = C_phi1\C_phi;
    
    %% Closing the loop
    
    % Allocating memory
    u_opt = zeros(size(phi,1),size(phi,2));
    eps = zeros(size(phi,1),size(phi,2));
    eps_hat = zeros(size(phi,1),size(phi,2)+1);
    e = zeros (size(G,2), size(phi,2));
    
    % Define matrix required for least squares solution 
    M = inv(H);
    % Compute Kalman gain
    [A, Cw, K] = computeKalmanAR(C_phi, C_phi1, G, sigma_e);
    
    % Simulate the initial conditions of closed loop
    e = random('norm', zeros(size(G,1),1), sigma_e * ones(size(G,1),1)); 
    eps(:,1) = phi(:,1);
    s = G*eps(:,1) + e;
    u_opt(:,1) = M*((A-K*G)*eps_hat(:,1) + K*s);
    eps_hat(:,2) = (A-K*G)*eps_hat(:,1) - H*u_opt(:,1) + K*s;
    
    % Close the loop and simulate for the remaining time steps
    for i = 2:size(phi,2)
        eps(:,i) = phi(:,i) - H*u_opt(:,i-1);
        e = random('norm', zeros(size(G,1),1), sigma_e * ones(size(G,1),1)); 
        s = G*eps(:,i) + e;
        u_opt(:,i) = M*((A-K*G)*eps_hat(:,i) + A*H*u_opt(:,i-1) + K*s);
        eps_hat(:,i+1) = (A-K*G)*eps_hat(:,i) - H*u_opt(:,i) + A*H*u_opt(:,i-1) + K*s;
    end

    % Remove the mean of residual
    eps = detrend(eps, 'constant');
    
    % Compute the variance
    var_eps = mean(var(eps));
end
