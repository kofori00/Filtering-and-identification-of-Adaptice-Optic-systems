function [var_eps] = AOloopMVM(G, H, C_phi_0, sigma_e, phi, lambda);
% Function for closing the loop based on random walk wavefront model
% INPUT:
%   G - matrix relation between phi and slopes measurement s
%   H - matrix relation between input u(k-1) and deformable mirror induced
%       wavefront phi_dm(k)
%   C_phi, C_phi1 - covariace matrices of wavefront phi
%   sigma_e - measurement noise variance
%   phi - incoming turbulent wavefront used for simulation only
%   lambda - weighting factor of regularization problem
% OUTPUT:
%   var_eps - residual wavefron variance 
%
% SEE THE REPORT FOR THEORETICAL DESCRIPTION



    W_e = 1/sigma_e^2;
    
    C_phi = C_phi_0;
    W_phi = inv(sqrtm(C_phi)*sqrtm(C_phi)');
    
    u_opt = zeros(size(phi,1),size(phi,2));
    eps = zeros(size(phi,1),size(phi,2));
    e = zeros (size(G,2), size(phi,2));
    
    e = random('norm', zeros(size(G,1),1), sigma_e * ones(size(G,1),1)); 
    eps(:,1) = phi(:,1);
    s = G*eps(:,1) + e;
    u_opt(:,1) = inv(H'*W_phi*H)*H'*W_phi*inv(G'*W_e*G + lambda*W_phi)*G'*W_e*s;
    
    for i = 2:size(phi,2)
        eps(:,i) = phi(:,i) - H*u_opt(:,i-1);
        e = random('norm', zeros(size(G,1),1), sigma_e * ones(size(G,1),1)); 
        s = G*eps(:,i) + e;
        u_opt(:,i) = inv(H'*W_phi*H)*H'*W_phi*inv(G'*W_e*G + lambda*W_phi)*G'*W_e*s + u_opt(:,i-1);
    end

    eps = detrend(eps, 'constant');
    
    var_eps = mean(var(eps));
end
