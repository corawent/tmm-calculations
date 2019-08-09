%Edited from absorption function 1/15/18 Cora Went
%Returns generation in units of 1/nm. Electrons created/photons incident per z step.
%Absorption is photons absorbed/photons incident, here no z integration.

%Note: need to add angle-dependent absorption, this is only normal
%incidence
%Also need to add single interface case
%Add generation for substrate case

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
% subs=0; %0 means no infinite substrate, 1 means infinite substrate 
% r_amp=0.5; %reflection amplitude calculated from transfer matrix (r_amp^2=R)

function [R,T,A,absi,gen]=generation(lambda,thetai,h,n,pol,r_amp,subs)

%% Initialize matrices
theta=zeros(1,length(n)-1);
Fr=zeros(1,length(n)-1);
Ft=zeros(1,length(n)-1);
delta=zeros(1,length(n)-2);
absi=zeros(1,length(n)-2);
gen = cell(1,length(n)-2);
gen{1,length(n)-2} = []; %use cell array for generation b/c matrix of g(x) will be different length for each layer

%% Snell's law:
theta(1)=thetai*pi/180;
for a=1:length(n)-1
    theta(a+1)=real(asin(n(a)/n(a+1)*sin(theta(a))))-1i*abs(imag(asin(n(a)/n(a+1)*sin(theta(a)))));
end

%% Fresnel coefficients:
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

%% Phase shift factors:
for a=2:length(n)-1
    delta(a)=2*pi*h(a)/lambda*n(a)*cos(theta(a));
end

%% Build up transfer matrix, for each M calculate the coefficients and absorption
M=[1,0;0,1]; %start with unity matrix
E0=[1;r_amp]; %electric fields in layer 1 (air)

%Calculate absorption coefficient
alpha=4*pi*n/lambda;

%Calculate absorption. a=1 calculates absorption in layer 2 (first non-air
%layer)
for a=1:length(n)-1
    
    if a==1
        M=(1/Ft(1)*[1,Fr(1);Fr(1),1])*M; %no propagation matrix on first one
    else
        M=(1/Ft(a)*[1,Fr(a);Fr(a),1]*[exp(1i*delta(a)),0;0,exp(-1i*delta(a))])*M; %a=2 and above, propagation included
    end
    
    %calculate absorption in all middle layers
    if a<length(n)-1
        E=M*E0;
        
        Ai=E(1,1);
        Bi=E(2,1); %so here we're calculating A2, B2
        
        alphai=imag(alpha(a+1)); %use alphai(2)
        z=0:1:h(a+1); %integrate over height of layer, z=0 at the start of the layer
        gen{a}=zeros(1,length(z));
        
        if h(a+1)==0 %in case one layer height is set to zero
            absi(a)=0;
            %generation already set to zero
        else
            absi(a)=real(n(a+1))*alphai*trapz(z,abs(Ai)^2*exp(-alphai*z)...
                +abs(Bi)^2*exp(alphai*z)+2*real(Ai*conj(Bi)...
                *exp((1j*4*pi*real(n(a+1))/lambda)*z))); %here we're calculating Abs in layer 2
            gen{a}=real(n(a+1))*alphai*(abs(Ai)^2*exp(-alphai*z)...
                +abs(Bi)^2*exp(alphai*z)+2*real(Ai*conj(Bi)...
                *exp((1j*4*pi*real(n(a+1))/lambda)*z))); %and Gen in layer 2
        end
    end
    
    %replace with correct absorption if there's a substrate 
    %exponent will blow up for very large z
    if subs==1
        absi(length(n)-2)=real(n(a))*abs(Ai)^2;
    end
    
end

%% Total Fresnel coefficients:
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

end