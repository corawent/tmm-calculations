function [rgb,xyzval]=spectrumtocolor(wl,r,em) %wl in nm, spectrum, em=0 for reflective and em=1 for emissive

%get x, y, z values - CIE color matching spectra
[~, ~, ~, ~, xyz] = get_xyz();

%constants
c=2.998e8;
h=6.626e-34;
k=1.381e-23;
t=6500;

%emissive or reflective?
if em==0
    %halogen illuminator, ish
    %ill=100*(560./wl).^5.*(exp(1.435e7/(2848*560))-1)./(exp(1.435e7./(2848*wl))-1);
    
    %6500K blackbody illuminator
    ill=2*h*c^2./(wl*10^-9).^5*1./(exp(h*c./(k*t*wl*10^-9))-1)*10^-9;
    
elseif em==1
    ill=ones(size(wl));
    
end

%initialize arrays
xyznew=zeros(length(wl),4);
xyznew(:,1)=wl;
xyzval=zeros(1,3);

for i=1:3
    grid=griddedInterpolant(xyz(:,1),xyz(:,i+1),'linear','none');
    xyznew(:,i+1)=grid(wl);
    xyznew(isnan(xyznew))=0;
    xyzval(1,i)=trapz(wl,xyznew(:,i+1).*ill.*r);
end

%emissive or reflective?
if em==0
    xyzval=xyzval/trapz(wl,xyznew(:,2).*ill); %for reflective thing, normalize by integrated ybar  
elseif em==1
    xyzval=xyzval/xyzval(:,2); %for blackbody, luminosity=1  
end

%xyzval=xyzval/sum(xyzval);

%convert to rgb
rgb=xyz2rgb(xyzval);
rgb(rgb>1)=1;
rgb(rgb<0)=0;

end