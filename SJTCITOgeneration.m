clc;
clear;
close all;

%% inputs
%thicknesses
tito=100;
tmet1=5;
ttmd=32;
tmet2=100;

%which TMD/metal?
tmd_ind=3; %1,2,3,4 are mos2, ws2, mose2, wse2
met1_ind=1; %1,2,3,4,5,6 are au, ag, pt, in, pd, al; met1=top
met2_ind=2; %met2=bottom

%thetai=0:5:90; %angle of incidence (degrees)
thetai=0;
lambdai=transpose(400:1:800); %vacuum wavelength (nm)
h=[NaN,tito,tmet1,ttmd,tmet2,NaN]; %film thicknesses in nm, equal in length to n, start and end with NaN (met1/tmd/met2)
pol=0; %polarization, 1 for p and 0 for s
subs=0; %signifies whether final layer should be treated as a substrate (infinite thickness)

%% load and initialize
%which index is which
tmds={'MoS_2','WS_2','MoSe_2','WSe_2'};
metals={'Au','Ag','Pt','In','Pd','Al'};

%initialize r, t, a, g
r=zeros(length(lambdai),1);
t=zeros(length(lambdai),1);
a=zeros(length(lambdai),1);
a_tmd=zeros(length(lambdai),length(h)-2);
g_tmd=zeros(length(lambdai),h(4)+1);

%load data
load('nkdata'); %pd and pt from werner;
load('am1p5_raw');

%% select materials
%select correct n, k data for tmd
if tmd_ind==1
    tmd_raw=mos2_raw;
elseif tmd_ind==2
    tmd_raw=ws2_raw;
elseif tmd_ind==3
    tmd_raw=mose2_raw;
elseif tmd_ind==4
    tmd_raw=wse2_raw;
end

%select correct n, k data for metal 1
if met1_ind==1
    met1_raw=au_raw;
elseif met1_ind==2
    met1_raw=ag_raw;
elseif met1_ind==3
    met1_raw=pt_raw;
elseif met1_ind==4
    met1_raw=in_raw;
elseif met1_ind==5
    met1_raw=pd_raw;
elseif met1_ind==6
    met1_raw=al_raw;
end

%select correct n, k data for metal 2
if met2_ind==1
    met2_raw=au_raw;
elseif met2_ind==2
    met2_raw=ag_raw;
elseif met2_ind==3
    met2_raw=pt_raw;
elseif met2_ind==4
    met2_raw=in_raw;
elseif met2_ind==5
    met2_raw=pd_raw;
elseif met2_ind==6
    met2_raw=al_raw;
end

%% interpolate raw nk data,am1p5
ito=nk_interp(ito_raw,lambdai);
met1=nk_interp(met1_raw,lambdai);
tmd=nk_interp(flip(tmd_raw),lambdai);
met2=nk_interp(met2_raw,lambdai);

am1p5(:,1)=lambdai;
am1p5(:,2)=interp1(am1p5_raw(:,1),am1p5_raw(:,2),am1p5(:,1),'linear','extrap');

clear('au_raw','ag_raw','pt_raw','pd_raw','in_raw','met1_raw','met2_raw',...
    'tmd_raw','mos2_raw','ws2_raw','mose2_raw','wse2_raw','am1p5_raw',...
    'si_raw','sio2_raw','al_raw','ito_raw')

%% loop over wavelengths and run transfer matrix, then generation function
for i=1:length(lambdai)
    lambda=lambdai(i);
    
    %initialize n for a given lambda
    n=ones(1,length(h));
    n(1,2)=ito(i,2)+1i*ito(i,3);
    n(1,3)=met1(i,2)+1i*met1(i,3);
    n(1,4)=tmd(i,2)+1i*tmd(i,3);
    n(1,5)=met2(i,2)+1i*met2(i,3);
    
    %for j=1:length(thetai)
    
    [R,T,A,r_amp]=transfermatrix(lambda,thetai,h,n,pol);
    
    r(i,1)=R;
    t(i,1)=T;
    a(i,1)=A;
    
    [R,T,A,absi,gen]=generation(lambda,thetai,h,n,pol,r_amp,subs); %note: thetai must=0 for now
    
    a_tmd(i,1:4)=absi(1:4);
    g_tmd(i,:)=gen{3};
    
    % end
    
end

%% calculate generation in active layer & save
const_c=3e8;
const_h=6.63e-34;
const_q=1.6e-19;
integrand=g_tmd.*am1p5(:,2).*lambdai/(const_h*const_c);
activegen=trapz(lambdai,integrand*10^-6); %nm*(1/nm)*(W/(m^2*nm))/(W*s*m/nm)=1/m^3*1/s*1m^3/(10^2)^3cm^3=1/cm^3*1/s

%save generation in format for lumerical device
x=(0:0.01:0.1)*10^-6;
y=(0:0.01:0.1)*10^-6;
z=(0:1:h(4))*10^-9;
G=permute(repmat(repmat(flip(activegen),length(x),1),[1,1,length(y)]),[1,3,2])*10^6; %turn 1D G data into 3D, multiply by 10^6 for units of 1/(m^3*s)
save('activegen.mat','G','x','y','z')

%save active layer absorption to use for IQE calculation
activeabs=zeros(length(lambdai),2);
activeabs(:,1)=lambdai;
activeabs(:,2)=a_tmd(:,3);
%save('activeabs.mat','activeabs');

%check Jsc match
int=am1p5(:,2).*activeabs(:,2)*const_q/(const_h*const_c).*lambdai*10^(-9);
jsc1=trapz(lambdai,int)/10;
jsc2=1000*trapz(z*10^9,activegen*10^-7*const_q); %mA/A*carriers/(cm^3*s)*C/carrier*(nm/cm)*nm=mA/cm^2

%% calculate color of reflected light, plot
em=0;
%[rgb,xyzval]=spectrumtocolor(lambdai,r,em);

%% plots
%
figure(1)
plot(lambdai,[a_tmd,sum(a_tmd,2),am1p5(:,2)./max(am1p5(:,2))])
title('Simulated Absorption')
xlabel('Wavelength [nm]')
ylabel('Absorption')
legend('ITO',metals{met1_ind},tmds{tmd_ind},metals{met2_ind},'Total','Solar')
ylim([0,1])
formatpresplot(8,8)
%export_fig 20nmAu35nmmose2_absorption -transparent -m4

% figure(2)
% plot(lambdai,r)
% title('TMD Color')
% xlabel('Wavelength [nm]')
% ylabel('R')
% formatpresplot
% legend off
%
% lines = findobj(gcf, 'type', 'line');
% for i=1:length(lines)
%     set(lines(i),'Color',[0,0,0])
% end
% set(gca,'Color',rgb)

% figure(3)
% plot(lambdai,a,lambdai,sum(a_tmd,2))
% title('test')
% xlabel('Wavelength [nm]')
% ylabel('A')
% legend('A','Sum')
% formatpresplot
% export_fig test -transparent -m4
%
% figure(4)
% plot(lambdai,[r,t,a])
% title('20nm Au')
% xlabel('Wavelength [nm]')
% ylabel('R, T, or A')
% legend('R','T','A')
% formatpresplot
% export_fig 20nmAu -transparent -m4

% figure(5)
% plot(met1(:,1),[met1(:,2),met1(:,3)])
% title(metals{met1_ind})
% xlabel('Wavelength [nm]')
% ylabel('Optical Constant')
% legend('n','k')
% formatpresplot
% export_fig(strcat(metals{met1_ind},'.png'),'-transparent','-r300')

% figure(6)
% plot(z*10^9,activegen)
% title('Generation in WS_2')
% xlabel('x [nm]')
% ylabel('Generation [1/(cm^3s)')
% formatpresplot
% legend off
% %export_fig generation -transparent -m4

%% display current
disp(strcat('Jsc:',num2str(jsc1),' mA/cm^2'))


%%

%% function for nk interpolation
function nk=nk_interp(nk_raw,lambdai)

G1=griddedInterpolant(nk_raw(:,1),nk_raw(:,2));
G2=griddedInterpolant(nk_raw(:,1),nk_raw(:,3));

nk(:,1)=lambdai;
nk(:,2)=G1(nk(:,1)/1000);
nk(:,3)=G2(nk(:,1)/1000);

end

