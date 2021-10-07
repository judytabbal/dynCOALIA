function [F,S]=P_SAUD(X,tau,P,Nit,alpha_min,alpha_max)
%computes the ICA of X using the P-SAUD algorithm
%
%[F,S]=P_SAUD(X,tau,P,Nit,alpha_min,alpha_max)
%
%computes the ICA of X where X is not centered using the P-SAUD algorithm, 
%which is based on the constrast function
%psi=cum4(s_i)^2+cum_4(s_j)^2+lambda*cov(s_i(t)*s_i(t-tau))^2.
%The penalization permits to ensure that sources with high autocorrelation
%are extracted first. Starting with an initial, high penalization parameter
%alpha_max, the influence of the penalization term is gradually decreased
%over the number of iterations. The final penalization parameter alpha_min
%is used for the last iterations.
%
%INPUT: X - data matrix
%       optional:
%       tau - delay between signals that is considered for the penalization
%             term (default: 1)
%       P - dimension of the signal subspace, permits to reduce the
%           dimension of the problem by truncating the SVD at the 
%           prewhitening step (default: min(size(X)))
%       Nit - number of iterations (default:20)
%       alpha_min - final relative penalization parameter (default:0)
%       alpha_max - inital relative penalization parameter (default:5)
%       
%
%OUTPUT: F - estimated mixing matrix
%        S - estimated signal matrix
%
% Hanna Becker, January 2013

%set default parameters
if nargin<6
    %initial relative penalization parameter
    alpha_max=5;
    if nargin<5
        %final relative penalization parameter
        alpha_min=0;
        if nargin<4
            %number of iterations
            Nit=20;
            if nargin<3
                %number of sources to extract
                P=size(X,1);
                if nargin<2
                    %delay
                    tau=1;
                end
            end
        end
    end       
end

%get size of the data matrix
[N,K]=size(X);

%center data
Xm=X-mean(X,2)*ones(1,K);

%prewhitening
[V,Sigma,U]=svd(Xm,'econ');
F=V(:,1:P)*(Sigma(1:P,1:P));
Z=pinv(F)*X;
Zm=pinv(F)*Xm;

%initialization of the matrix G_tilde
G1=eye(P);

%set penalization strategy (relative penalization parameter is gradually
%decreased)
if Nit>3
    lambda1=fliplr(linspace(alpha_min,alpha_max,Nit-3));
    lambda1=[lambda1 alpha_min*ones(1,3)];
%     lambda1=alpha_max*ones(1,Nit);
else
    lambda1=[alpha_max, alpha_min*ones(1,Nit-1)];
end

%loop over number of sources
for n=1:min(P,N-1)
        %compute matrix G for the n-th source
        Z1=G1*Zm;
        
        N1=size(Z1,1);

        %initialization of the rotation matrix
        G=eye(N1);

        %loop over number of iterations that are used to update the rotation matrix
        for nit=1:length(lambda1)
            %relative penalization parameter
            lambda=lambda1(nit);

            %modified data matrix that is considered for the nit-th iteration
            S1=G*Z1;
            
            %determine kurtosis of current source estimate
            y2=S1(1,:).^2;
            Cyyyy=(1/K*sum(y2.*y2)-3*(1/K*sum(y2))^2)^2;
    
            %determine penalization term for current source estimate
            for k=1:length(tau)
                alpha1(k)=1/(K-tau(k))*(S1(1,1:end-tau(k))*S1(1,tau(k)+1:end).');
            end
            pen=sum(alpha1.^2);
            
            %adjust penalization parameter based on kurtosis and
            %autocorrelation of the current source estimate
            lambda=Cyyyy/pen*lambda1(nit);

            %loop over reference sources
            for k=2:N1
                %compute optimal value for the rotation angle that
                %maximizes the contrast
                beta=get_optimal_rotation_angle_mult_cov(S1(1,:),S1(k,:),lambda,tau);
                beta=[cos(beta) sin(beta); -sin(beta) cos(beta)];

                %update data matrix and rotation matrix
                Gg=eye(N1);
                Gg(1,1)=beta(1,1);
                Gg(1,k)=beta(1,2);
                Gg(k,1)=beta(2,1);
                Gg(k,k)=beta(2,2);

                %update data matrix and rotation matrix
                G=Gg*G;
                S1=Gg*S1;
            end
        end

        %compute n-th vector of the orthogonal mixing matrix
        Gs(:,n)=G1.'*G(1,:).';

        %determine n-th signal estimate
        S(n,:)=Gs(:,n).'*Z;

        %update G_tilde
        G1=G(2:end,:)*G1;

end
if P==N
    %compute last vector of the orthogonal mixing matrix
    Gs(:,N)=G1.';
    S(N,:)=Gs(:,N).'*Z;
end

%compute complete mixing matrix
F=F*Gs;