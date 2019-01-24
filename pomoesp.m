function [Aest, Cest, K, vaf] = pomoesp(sk, Nid, Nval, s, n)
% POMOESP
% Not really po-moesp, but a combination with N4SID. 
% WARNING: WRITTEN SPECIALLY FOR PRACTICAL, DO NOT USE FOR GENERAL SYSTEMS
% Estimates A, C and Kalman gain matrices of an autonomous system from data.
%
% INPUT:
%   sk - data
%   Nid - amount of data used for identification (rule of thumb: 2/3 of
%   total data)
%   Nval - data used for validation
%   s - number of block rows of Hankel matrix
%   n - desired output system order
% OUTPUT:
%   Aest - estimated n-by-n A matrix
%   Cest - estimated (number-of-outputs)-by-n matrix
%   K - Kalman gain matrix, n-by-(number-of-outputs)
%   vaf - variance-accounted-for - the more, the better
%
% For theoretical description, see book "Filtering and System
% Identification: A Least Squares Approach (M. Verhaegen, V. Verdult),
% pages 332-334.


%% 1. Detrend data

skdecd=sk;
for k=1:length(sk(:,1))
    skdecd(k,:)=detrend(skdecd(k,:));
end  

sk = skdecd;



%% 2. Hankel matrix and instrumental variable
length_out = size(sk,2);

% Sanity check - don't use more data then you have
if (Nid + Nval > length_out)
    error('Number of data points for identification PLUS number of data points for validation must be less than supplied data length');
end

% Number of outputs of the identified system - denoted by l in the book
nout = size(sk,1);

% Construct a big Hankel matrix of outputs. Hankel matrix will have 2*s
% block-rows, and Nid - 2*s +1 columns
y_hankel_outputs = zeros(nout*2*s, (Nid-2*s+1));
    for i = 1:Nid-2*s+1
       y_hankel_outputs(:,i) = reshape(sk(:,i:i+2*s-1), 2*s*nout,1); 
    end
    
 % Get dimensions of output Hankel matrix - just in case
[y1, y2] = size(y_hankel_outputs);

% Instrumental variable - used to eliminate inovation signal from data
% equation
instr_var = y_hankel_outputs(1:s*nout,:);


%% 3. RQ-decomposition 

% Perform RQ decomposition of the big Hankel matrix. You don't really use
% Q, so don't save it for the sake of saving memory.
R=triu(qr(y_hankel_outputs'))';

% Extract the matrices you need for further calculations
R11 = R(1:nout*s, 1:nout*s);
R21 = R(nout*s+1:end, 1:nout*s);


%% 4. SVD decomposition

% Perform SVD decomposition - will be used to calculate state sequence and
% system matrices. You can also approximate system order by looking at how
% many nonzero singular values there are in S matrix. U is not used in this
% case, so don't save it.

[~, S, V] = svd(R21*(R11\eye(length(R11)))*instr_var);


%% 5. Estimation of system matrices

% Estimated state sequence
Xest = sqrtm(S(1:n, 1:n))*V(:,1:n)';

% Extract its dimensions
[x1, x2] = size(Xest);

% Construct two chopped matrices of states - used further for least squares
Xest1 = Xest(:, 2:x2);
Xest0 = Xest(:,1:x2-1);

% Cut a part of output Hankel matrix
Ycut = y_hankel_outputs(s*nout+1:(s+1)*nout,1:y2-1);

% Estimated state matrices from least squares problem. Used solution for
% Frobenius norm - see page 302 of the book, top of the page.
sysest = [Xest1;Ycut]*pinv(Xest0);

% Extract system matrices
Aest = sysest(1:n, 1:n);
Cest = sysest(n+1:n+nout,1:n);


%% 6. Estimation of covariances and Kalman gain

resid_est = [Xest1;Ycut] - sysest*Xest0;
covest = 1/(Nid-1) * (resid_est*resid_est');
% Covariance of process noise
Qest = covest(1:n, 1:n);
% Covariance of measurement noise
Rest = covest(n+1:n+nout, n+1: n+nout);
% Cross-covariance
Sest = covest(1:n, n+1: n+nout);

% Kalman gain calculation - see documentation of dare() and compare it to
% equation on page 334 of the book to see why exaxtly A' and C' are sent to
% the function
[~, ~, K] = dare(Aest', Cest', Qest, Rest, Sest);

K = K';


%% 7. Simulation and verification

phi_sim = zeros(n,Nval);
s_sim = zeros(nout, Nval);
for i=1:Nval
    phi_sim(:,i+1) = (Aest - K*Cest)*phi_sim(:,i) + K*sk(:,Nid + i);
    s_sim(:,i) = Cest*phi_sim(:,i);
end

% Compute variance-accounted-for
vaf = computevaf(s_sim, sk(:,Nid+1:Nid+Nval));


end
