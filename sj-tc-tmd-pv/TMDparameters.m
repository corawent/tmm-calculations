clc;
clear;
close all;

%% TMDs
ntmds=4;
tmds={'mos2','ws2','mose2','wse2'};

%% constants
k=1.381e-23; %J/K
t=300; %K
h=6.626e-34; %J/s
m=9.11e-31; %kg
q=1.602e-19; %coulombs
c=2.998e8; %m/s

%kT=0.0257; %eV

%% initialize structure
inputs={'eg','mel','met','mh','doping','ei','eps0','mu','qy'};
outputs={'me','ni','wf','brad','tnr'};

tmd(ntmds).input(length(inputs))=0;
tmd(ntmds).output(length(outputs))=0;
tmd(4).table=zeros(1,9);

%% enter inputs

eg=[1.23,1.35,1.09,1.20]; %eV, Kam
mel=[0.590,0.569,0.521,0.489]; %me, Wickramarante
met=[0.845,0.665,0.776,0.643]; %me, Wickramarante
mh=[0.838,0.832,0.973,0.997]; %me, Wickramarante
doping=[1e15,1e14,1e15,1e15]; %1/cm^3, HQ Graphene
ei=[5.44,5.26,5.08,5.08]; %eV, Keyshar CHANGE for WSe2
eps0=[7.3,6.7,8.5,7.9]; %Laturia
mu=0.01*ones(1,4); %cm^2/s, Massicotte
qy=1e-6*ones(1,4); %CHANGE

all=[eg;mel;met;mh;doping;ei;eps0;mu;qy];

%% populate structure with inputs

for i=1:ntmds
    tmd(i).inputs=inputs;
    tmd(i).outputs=outputs;
    
    for j=1:length(inputs)
        tmd(i).input(j)=all(j,i);
    end
end

%% optical constants

ev=transpose(1:0.01:3);

load('nkdata_old'); %wl in units of microns, n, k
load('ftpsdata'); %eV, something prop. to alpha

tmd(1).nk=mos2_raw;
tmd(2).nk=ws2_raw;
tmd(3).nk=mose2_raw;
tmd(4).nk=wse2_raw;

tmd(1).ftps=mos2_ftps;
tmd(2).ftps=ws2_ftps;
tmd(3).ftps=mose2_ftps;
tmd(4).ftps=wse2_ftps;

for i=1:4
    
    %set matching point and bounds for quadfit for each tmd
    if i==1
        matching=1.8;
        bound_low=1.4;
        bound_high=1.7;
        
    elseif i==2
        matching=1.8;
        bound_low=1.5;
        bound_high=1.6;
        
    elseif i==3
        matching=1.45;
        bound_low=1.25;
        bound_high=1.3;
        
    elseif i==4
        matching=1.5;
        bound_low=1.27;
        bound_high=1.32;
        
    end
    
    %interpolate low energy (ftps) and high energy (nk) data
    ev_raw_high=1.240./tmd(i).nk(:,1);
    alpha_raw_high=4*pi*tmd(i).nk(:,3)./(tmd(i).nk(:,1)*10^-4); %cm^-1
    
    ev_raw_low=tmd(i).ftps(:,1);
    alpha_raw_low=tmd(i).ftps(:,2);
    
    highe=griddedInterpolant(ev_raw_high,alpha_raw_high);
    lowe=griddedInterpolant(ev_raw_low,alpha_raw_low);
    
    alpha_highe=highe(ev);
    alpha_lowe=lowe(ev);
    
    %match low and high energy data at energy specified by "matching"
    factor=alpha_highe(ev==matching)/alpha_lowe(ev==matching);
    alpha_lowe=factor*alpha_lowe;
    alpha=zeros(length(ev),1);
    alpha(ev<=matching,1)=alpha_lowe(ev<=matching);
    alpha(ev>matching,1)=alpha_highe(ev>matching);
    
    eg=tmd(i).input(contains(inputs,'eg'));
    
    %quadfit that is set to zero at bandgap between given bounds
    quad=fittype(strcat('a*(x-',num2str(eg),')^2'),'independent','x');
    quadtail=fit(ev(ev>bound_low&ev<bound_high),alpha(ev>bound_low&ev<bound_high),...
        quad,'Startpoint',0.01);
    alpha(ev<bound_low)=quadtail(ev(ev<bound_low));
    
    %set alpha to zero below bandgap (ideal material...)
    alpha(ev<eg)=0;
    tmd(i).alpha=alpha; %gives alpha in cm^-1

end

%% calculate outputs
for i=1:4
    
    %me
    mel=tmd(i).input(contains(inputs,'mel'));
    met=tmd(i).input(contains(inputs,'met'));
    me=(mel*met^2)^(1/3);
    tmd(i).output(1)=me;
    
    %ni
    eg=tmd(i).input(contains(inputs,'eg'));
    mh=tmd(i).input(contains(inputs,'mh'));
    nc=2*(2*pi*me*m*k*t/(h^2))^(3/2)*1/(100)^3; %1/cm^3
    nv=2*(2*pi*mh*m*k*t/(h^2))^(3/2)*1/(100)^3; %1/cm^3
    ni=sqrt(nc*nv*exp(-eg*q/(k*t))); %1/cm^3
    tmd(i).output(2)=ni;
    
    %wf = fermi level of undoped semiconductor
    ei=tmd(i).input(contains(inputs,'ei'));
    ef=0.5*(2*ei-eg)+0.75*k*t/q*log(mh/me);
    tmd(i).output(3)=ef;
    
    %brad
    alpha=tmd(i).alpha;
    n=tmd(i).nk(1,2);
    bb=2*(q*ev).^2./(h^3*c^2).*1./(exp(q*ev./(k*t))-1);
    integrand=4*pi*(alpha*100).*n.^2.*bb/100^3; %multiply alpha by 100 to get m^-1, divide by 100^3 to go back to cm^-3
    rrad=trapz(q*ev,integrand);
    brad=rrad/ni^2;
    tmd(i).output(4)=brad;
    
    %tnr
    doping=tmd(i).input(contains(inputs,'doping'));
    qy=tmd(i).input(contains(inputs,'qy'));
    tau_r=1/(brad*doping);
    tau_nr=qy*tau_r/(1-qy);
    tau=1/(1/tau_nr+1/tau_r);
    tmd(i).output(5)=tau_nr;

    %fill in table values
    tmd(i).table(1)=eg;
    tmd(i).table(2)=ef;
    tmd(i).table(3)=me;
    tmd(i).table(4)=mh;
    tmd(i).table(5)=doping;
    tmd(i).table(6)=tmd(i).input(contains(inputs,'mu'));
    tmd(i).table(7)=tmd(i).input(contains(inputs,'eps0')); 
    tmd(i).table(8)=brad;
    tmd(i).table(9)=tau_nr;
    
    tmd(i).table=transpose(tmd(i).table);

end

%% create table 

%column names
cols={'MoS2', 'WS2', 'MoSe2', 'WSe2'};

%row names
rows={'Bandgap (eV)','Work Function (eV)','Electron Eff Mass (*me)', ...
    'Hole Eff Mass (*me)','Doping (cm^-3)','OOP Mobility (cm^2/Vs)',...
    'DC permittivity','Brad (cm^3/s)','SRH Liftime (s)'};

table(tmd(1).table,tmd(2).table,tmd(3).table,tmd(4).table,...
    'RowNames',rows,'VariableNames',cols)


%% clear unnecessary stuff 

clear('au_raw','ag_raw','mos2_raw','ws2_raw','mose2_raw','wse2_raw',...
    'ev_raw_high','ev_raw_low','alpha_raw_high','alpha_raw_low')