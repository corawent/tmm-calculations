%Edited from transfer matrix function 12/18/18 Cora Went
%Returns R, T, A and reflection amplitude r_amp (R=r_amp^2).

%inputs:
%angle of incidence: thetai (degrees) 
%wavelength of incident light: lambda (nm)
%thicknesses of the layers: h
%refractive index of the layers: n (may be absorbing and dispersive, but this may require a subfunction)
%polarization: s or p (need one calculation for each for unpolarized light)

% thetai=0:3.0:90; %angle of incidence (degrees)
% lambdai=1200:10:1440; %vacuum wavelength (nm)
% h=[NaN,50,NaN,NaN]; %film thicknesses in nm, equal in length to n, start and end with NaN
% pol=1; %polarization, 1 for p and 0 for s

function [R,T,A,r_amp]=transfermatrix(lambda,thetai,h,n,pol)

%Initialize matrices
theta=zeros(1,length(n)-1);
Fr=zeros(1,length(n)-1);
Ft=zeros(1,length(n)-1);
delta=zeros(1,length(n)-2);

%Snell's law:
theta(1)=thetai*pi/180;
for a=1:length(n)-1
    theta(a+1)=real(asin(n(a)/n(a+1)*sin(theta(a))))-1i*abs(imag(asin(n(a)/n(a+1)*sin(theta(a)))));
end

%Fresnel coefficients:
if pol==0 %formulas for s polarization
    for a=1:length(n)-1
        Fr(a)=(n(a+1)*cos(theta(a+1))-n(a)*cos(theta(a)))/(n(a)*cos(theta(a))+n(a+1)*cos(theta(a+1)));
        Ft(a)=2*n(a+1)*cos(theta(a+1))/(n(a)*cos(theta(a))+n(a+1)*cos(theta(a+1)));
    end
elseif pol==1 %formulas for p polarization
    for a=1:length(n)-1
        Fr(a)=(n(a+1)*cos(theta(a))-n(a)*cos(theta(a+1)))/(n(a)*cos(theta(a+1))+n(a+1)*cos(theta(a)));
        Ft(a)=2*n(a+1)*cos(theta(a+1))/(n(a)*cos(theta(a+1))+n(a+1)*cos(theta(a)));
    end
end

%Phase shift factors:
for a=2:length(n)-1
    delta(a)=2*pi*h(a)/lambda*n(a)*cos(theta(a));
end

%Build up transfer matrix
M=[1,0;0,1]; %start with unity matrix
for a=1:length(n)-1
    
    if a==1
        M=(1/Ft(1)*[1,Fr(1);Fr(1),1])*M; %no propagation matrix on first one
    else
        M=(1/Ft(a)*[1,Fr(a);Fr(a),1]*[exp(1i*delta(a)),0;0,exp(-1i*delta(a))])*M; 
    end

end

%Total Fresnel coefficients:
Frtot=-M(2,1)/M(2,2);
Fttot=M(1,1)-M(1,2)*M(2,1)/M(2,2);

%special case of single interface:
if length(n)==2
    Frtot=Fr(1);
    Fttot=Ft(1);
end

%total Fresnel coefficients in intensity:
R=(abs(Frtot))^2;
T=(abs(Fttot))^2*real(n(length(n))*cos(theta(length(n))))/real(n(1)*cos(theta(1)));
A=1-R-T;

%reflection coefficient amplitude
r_amp=Frtot;

end
