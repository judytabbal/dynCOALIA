function [ newCov,iW ] = TruncateAndRegularizeCovariance( Cov )

% needed method for wMNE computing

% This code was originally developped by Sahar Yassine and Ahmad Mheich based on the brainstorm toolbox codes 
% (Tadel, F., Baillet, S., Mosher, J.C., Pantazis, D., Leahy, R.M., 2011.
% Brainstorm: A user-friendly application for MEG/EEG analysis. Computational Intelligence and Neuroscience 2011. https://doi.org/10.1155/2011/879716)

% contact: saharyassine94@gmail.com
%          mheich.ahmad@gmail.com

newCov = (Cov + Cov')/2;
[Un,Sn2] = svd(newCov,'econ');
Sn = sqrt(diag(Sn2)); % singular values
tol = length(Sn) * eps(single(Sn(1))); % single precision tolerance
Rank_Noise = sum(Sn > tol);
Un = Un(:,1:Rank_Noise);
Sn = Sn(1:Rank_Noise);
newCov = Un*diag(Sn.^2)*Un';

%diag method
% newCov = diag(diag(newCov)); % strip to diagonal
% iW = diag(1./sqrt(diag(newCov)));

%reg Method
NoiseReg=0.1;
RidgeFactor = mean(diag(Sn2)) * NoiseReg;
newCov = newCov + RidgeFactor * eye(size(newCov,1));
iW = Un*diag(1./sqrt(Sn.^2 + RidgeFactor))*Un'; 

end

