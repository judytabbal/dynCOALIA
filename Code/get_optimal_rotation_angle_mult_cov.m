function theta=get_optimal_rotation_angle_mult_cov(y,z,lambda,tau)
%comp. optimal rotation angle for indep. of sources y and z with constraint
%
%theta=get_optimal_rotation_angle_mult_cov(y,z,lambda,tau)
%
%computes the optimal rotation angle for the Givens rotation that makes the
%vectors y and z maximally independent using the penalized constrast
%psi = cum_4(s1(t))^2 + cum_4(s2(t))^2 + lambda*cov(s1(t)*s1(t-tau))^2. 
%
%INPUT: y - data vector
%       z - data vector
%       lambda - regularization parameter
%       optional:
%       tau - considered delay (default: 1)
%       
%OUTPUT: theta - optimal rotation angle
%
% Hanna Becker, January 2013

if nargin <4
    %default delay
    tau = 1;
end

%determine number of time samples of the data vectors
K=length(y);

%compute required 4-th order cumulants of the data vectors
y2=y.^2;
z2=z.^2;
yz=y.*z;
Cyyyy=1/K*sum(y2.*y2)-3*(1/K*sum(y2))^2;
Czzzz=1/K*sum(z2.*z2)-3*(1/K*sum(z2))^2;
Cyyyz=1/K*sum(y2.*yz)-3*1/K*sum(y2)*1/K*sum(yz);
Cyyzz=1/K*sum(y2.*z2)-1/K*sum(y2)*1/K*sum(z2)-2*(1/K*sum(yz))^2;
Cyzzz=1/K*sum(z2.*yz)-3*1/K*sum(z2)*1/K*sum(yz);

%dertermine 3 sets of coefficients that are used in [1] to compute the
%contrast
a0=Cyyyy; a1=4*Cyyyz; a2=6*Cyyzz; a3=4*Cyzzz; a4=Czzzz;

%determine covariance coefficients for several penalization terms
for k=1:length(tau)
    alpha1(k)=1/(K-tau(k))*(y(1:end-tau(k))*y(tau(k)+1:end).');
    beta1(k)=1/(K-tau(k))*(y(1:end-tau(k))*z(tau(k)+1:end).'+z(1:end-tau(k))*y(tau(k)+1:end).');
    gamma1(k)=1/(K-tau(k))*(z(1:end-tau(k))*z(tau(k)+1:end).');
end
alpha=sum(lambda.*alpha1.^2);
beta=sum(lambda.*beta1.^2);
gamma=sum(lambda.*gamma1.^2);
delta=sum(lambda.*alpha1.*beta1);
epsilon=sum(lambda.*alpha1.*gamma1);
phi=sum(lambda.*gamma1.*beta1);

%compute polynomial coefficients of the contrast function
d8=a0^2+a4^2+gamma;
d7=2*(a3*a4-a0*a1+phi);
d6=a1^2+a3^2+2*a0*a2+2*a2*a4+2*epsilon+2*gamma+beta;
d5=2*(a1*a4+a2*a3-a0*a3-a1*a2+delta+2*phi);
d4=2*(a2^2+2*a0*a4+2*a1*a3)+alpha+2*beta+4*epsilon+gamma;
d3=-2*(a1*a4+a2*a3-a0*a3-a1*a2)+4*delta+2*phi;
d2=a1^2+a3^2+2*a0*a2+2*a2*a4+2*alpha+2*epsilon+beta;
d1=-2*(a3*a4-a0*a1)+2*delta;
d0=a0^2+a4^2+alpha;

%compute coefficients of the derivative of the contrast function
coeff(1)=-d7;
coeff(2)=-2*d6+8*d8;
coeff(3)=-3*d5+7*d7;
coeff(4)=-4*d4+6*d6;
coeff(5)=-5*d3+5*d5;
coeff(6)=-6*d2+4*d4;
coeff(7)=-7*d1+3*d3;
coeff(8)=-8*d0+2*d2;
coeff(9)=d1;

%determine roots of the derivative polynomial and find optimal rotation
%angle theta
theta1=roots(coeff).';
theta1=theta1(find(abs(imag(theta1))<1e-14));
psi2=polyval([d8 d7 d6 d5 d4 d3 d2 d1 d0],theta1)./(1+theta1.^2).^4;
[psi2max,idt]=max(psi2);
idt=find(psi2>=psi2max-psi2max*1e-4);
theta1=theta1(idt);
[~,idt]=min(abs(theta1));
theta=atan(theta1(idt));