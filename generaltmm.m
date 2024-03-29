clc;
clear;
close all;

%% inputs

%materials in stack (must correspond to name of variable in nkdata.mat)
mat_list=["SiO_2","Pt","WS_2","Ag"];

%thicknesses, starting and ending with NaN
h=[NaN,1700,10,80,100,NaN];

%angle of incidence
thetai=0;

%wavelengths
lambdai=transpose(400:1:800);

%polarization (1 for p and 0 for s)
pol=0;

%treat last layer as infinite substrate?
subs=0;

%% load and initialize

%initialize r, t, a, g
r=zeros(length(lambdai),1);
t=zeros(length(lambdai),1);
a=zeros(length(lambdai),1);
a_tmd=zeros(length(lambdai),length(h)-2);
g_tmd=zeros(length(lambdai),h(4)+1);

%load data
nk_raw=load('nkdata'); %pd and pt from werner;
load('am1p5_raw');

%% interpolate raw nk data,am1p5

for i=1:length(mat_list)
    mat=mat_list(i);
    nk.(mat)=nk_interp(nk_raw.(strcat(mat,'_raw')),lambdai);
end

am1p5(:,1)=lambdai;
am1p5(:,2)=interp1(am1p5_raw(:,1),am1p5_raw(:,2),am1p5(:,1),'linear','extrap');

clear('nk_raw')

%% loop over wavelengths and run transfer matrix, then generation function
for i=1:length(lambdai)
    lambda=lambdai(i);
    
    %initialize n for a given lambda
    n=ones(1,length(h));
    for j=1:length(mat_list)
        mat=mat_list(j);
        nk_mat=nk.(mat);
        n(1,j+1)=nk_mat(i,2)+1j*nk_mat(i,3);
    end
   
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

%% calculate color of reflected light
em=0;
[rgb,xyzval]=spectrumtocolor(lambdai,r,em);

%% plots
%filename
name='';
for j=1:length(mat_list)-1
    name=strcat(name,num2str(h(j+1)),'nm',char(mat_list(j)),'-');
end

figure(1)
plot(lambdai,[a_tmd,sum(a_tmd,2)])
title('Simulated Absorption')
xlabel('Wavelength [nm]')
ylabel('Absorption')
legend([mat_list,'Total'])
ylim([0,1])
formatpresplot(8,8)
%export_fig(strcat(name,'Abs.png'),'-transparent','-r300')
%
% figure(2)
% plot(lambdai,r)
% title('TMD Color')
% xlabel('Wavelength [nm]')
% ylabel('R')
% formatpresplot(8,8)
% % legend off
%
% lines = findobj(gcf, 'type', 'line');
% for i=1:length(lines)
%     set(lines(i),'Color',[0,0,0])
% end
% set(gca,'Color',rgb)

% figure(3)
% plot(lambdai,[r,t,a])
% title('20nm Au')
% xlabel('Wavelength [nm]')
% ylabel('R, T, or A')
% legend('R','T','A')
% formatpresplot
% export_fig 20nmAu -transparent -m4

% figure(4)
% mat="WS_2";
% nk_mat=nk.(mat);
% plot(nk_mat(:,1),[nk_mat(:,2),nk_mat(:,3)])
% title(mat)
% xlabel('Wavelength [nm]')
% ylabel('Optical Constant')
% legend('n','k')
% formatpresplot
% %export_fig(strcat(char(mat),'.png'),'-transparent','-r300')

% figure(6)
% plot(z*10^9,activegen)
% title('Generation in WS_2')
% xlabel('x [nm]')
% ylabel('Generation [1/(cm^3s)')
% formatpresplot
% legend off
%export_fig generation -transparent -m4

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

