%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% last update 29 November 2021, lne %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Have a look in the book "Multi-Band Efective Mass Approximations"
% p145
% 4.5.2 Ellipticity Criteria
% => There are explanations how to deal with the spurious solutions...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program solves the Schrodinger equation with m(E,z) using different algorithms:
% -> FDM (Finite Difference Method) 
%          and / or
% -> Scaning/Shooting method
%
% The non-parabolicity is implemented via the kp model using the Kane or Luttinger approach
% The HH band remains ALL THE TIME parabolic in ZB-001 growth and even with strain
% A strain model is included. It basically shifts the conduction and valence band edge
% but also influences the coupling bw LH/SO in the kp 6bands.
% The strain is mainly interesting for InGaAs/GaAs heterostructures

% to do:
% include the strain in the kp 8band
% check the kp 6bands with DKK notation
% check that HH is still parabolic if the confinemtent is along x-axis or y-axis

% -> II-VI and cubic nitride material parameters are available but should
% be grabt in the "Library.m" file
% -> Wurtzite parameter are also availables but the code isn't optimized for it.
% In the "Library.m" file, the WZ table must be open and the meaningfull parameters
% must be taken. Also, the electric field has to be handle...
% -> Additionnal material can be added in the "materialDB_ZB.csv" file

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant J.s
hbar=h/(2*pi);
e=1.602176487E-19;              %% charge de l electron Coulomb
m0=9.10938188E-31;              %% electron mass kg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 for turn off
% 1 for turn on

% Diagonalization of the Hamiltonian with Finite Difference Method (FDM)
FDM_kp1bandCB           = 0;    % mass=m(z) => parabolic band
FDM_kp1bandHH           = 0;    % Good to solve HH that is always parabolic even with strain
FDM_kp1bandLH           = 0;    % Good to solve LH if Dso is large otherwise not valid
FDM_kp2bands_Kane       = 1;    % Good to solve Ec only
FDM_kp3bands_Kane       = 0;    % Good to solve Ec only
FDM_kp6bands_Luttinger  = 1;    % Good to solve HH, LH and SO in 1+2 band because HH is not coupled to the other bands in ZB-001 (even with strain)
FDM_kp8bands_Luttinger  = 0;    % NOT WELL WORKING!!! "spurious" solution but electron solutions are too high as well. Increasing the resolution dz helps

% Shooting/Scanning in Energy method (Euler)
ShootingCB_1band_Kane   = 0;    % mass=m(z) => parabolic band
ShootingHH_1band        = 0;    % Good to solve HH that is always parabolic in ZB-001, even with strain
ShootingCB_2bands_Kane  = 0;    % Good to solve Ec
ShootingCB_3bands_Kane  = 0;    % Good to solve Ec
ShootingCB_Luttinger    = 0;    % Good to solve Ec but gives different results compare to Kane s model because the mass is slightly different

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

StrainModel     = 0;            % Activate Strain model
PrintResultsISB = 0;            % Switch to print or not the ISB dipoles on the shell
PlotMass        = 1;            % PLot meff in the Kane and Luttinger model
PlotVB          = 1;            % plot the valence band

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=3;                     % number of solution asked per model
ScF=0.1;                  % scaling factor to plot the wave function [Without Dimension]
dz=1E-10;                 % resolution of the grid [m]
F0=0;%-6e7;               % Electric field [Volt/meter]
T=300;                    % Temperature [Kelvin], react on the band gap Eg only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Library;                  % load material parameter DB from "materialDB_ZB.csv"
ExtractParameters;        % extract parameter from the Library
TernaryAlloy;             % compute the ternary alloy    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% import the layer structure file %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first column is the material used from the "library"
% second column is the length of the layer in nm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_file;

% substrate=GaAs;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
% M=[
% AlGaAs40      5
% GaAs     20
% AlGaAs40      5
% ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE !!! %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Grabbing the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zt  = M(:,end)*1E-9;    % conversion of the length from Angstrom to meter

Egt = M(:,idx_Eg6c) - (M(:,idx_alphaG)*T^2) ./ (T+M(:,idx_betaG));   %Eg = Eg0 - (a*T.^2)./(T + b);
VBOt= M(:,idx_VBO);
CBOt= Egt+VBOt;         % CBO form band gap difference and temperature
Dsot= M(:,idx_Dso);     % Spin-Orbit shift band parameter
Ft  = M(:,idx_F);       % Gammac Luttinger parameter for the electron
g1t = M(:,idx_g1);      % Gamma1 Luttinger parameter
g2t = M(:,idx_g2);      % Gamma2 Luttinger parameter
g3t = M(:,idx_g3);      % Gamma3 Luttinger parameter

EPt_K=M(:,idx_EP_K);    % EP Kane
EPt_L=M(:,idx_EP_L);    % EP Luttinger

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Strain Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

at  = M(:,idx_a);           % lattice parameter
act = M(:,idx_ac);          % Conduction band strain offset parameter
avt = M(:,idx_av);          % Valence band strain offset parameter
bvt = M(:,idx_bv);          % Valence band strain offset parameter
c11t = M(:,idx_c11);        % strain parameter
c12t = M(:,idx_c12);        % strain parameter

a0   = substrate(idx_a);

if StrainModel == 1
  exxt =  (a0-at)/a0; % eyyt =  exxt;
  ezzt = -2*c12t./c11t.*exxt;
else
  exxt =  (a0-at)/a0 * 0; % eyyt =  exxt;
  ezzt = -2*c12t./c11t.*exxt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Discretisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here, I descretize the grid z, the potential V0 and the mass me

z(1)=0;  V0(1)=CBOt(1); Eg(1)=Egt(1);Dso(1)=Dsot(1);F(1)=Ft(1);
g1(1)=g1t(1);g2(1)=g2t(1);g3(1)=g3t(1); EP_K(1)=EPt_K(1); EP_L(1)=EPt_L(1); %me(1)=met(1);% mhh(1)=mhh_t(1);
ac=act(1); av=avt(1); bv=bvt(1); exx=exxt(1); ezz=ezzt(1);

for i=1:length(zt)
    zv  =  (z(end)+dz) : dz : (z(end) + zt(i)) ;
    z   = [ z  zv ];
    V0  = [ V0     ones(size(zv)) * CBOt(i)  ];
    Eg  = [ Eg     ones(size(zv)) * Egt(i)   ];
    EP_K= [ EP_K   ones(size(zv)) * EPt_K(i) ];
    EP_L= [ EP_L   ones(size(zv)) * EPt_L(i) ];
    Dso = [ Dso    ones(size(zv)) * Dsot(i)  ];
    F   = [ F      ones(size(zv)) * Ft(i)    ];
    g1  = [ g1     ones(size(zv)) * g1t(i)   ];
    g2  = [ g2     ones(size(zv)) * g2t(i)   ];
    g3  = [ g3     ones(size(zv)) * g3t(i)   ];
    ac  = [ ac     ones(size(zv)) * act(i)   ];
    av  = [ av     ones(size(zv)) * avt(i)   ];
    bv  = [ bv     ones(size(zv)) * bvt(i)   ];
    exx = [ exx    ones(size(zv)) * exxt(i)  ];
    ezz = [ ezz    ones(size(zv)) * ezzt(i)  ];
end

V0=V0-min(V0);             % Shift the band in order to get the bottom of the well at zero
V0=(F0*z)+V0;              % adding the electric field to the potential

eyy = exx;
DCBO   = -abs(ac).*(exx+eyy+ezz) ;                      % shift of the CB due to strain
DVBOHH = +abs(av).*(exx+eyy+ezz) - abs(bv).*(exx-ezz) ; % shift of the VB-HH due to strain
DVBOLH = +abs(av).*(exx+eyy+ezz) + abs(bv).*(exx-ezz) ; % shift of the VB-LH due to strain
DVBOSO = +abs(av).*(exx+eyy+ezz) ;                      % shift of the VB-SO due to strain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Discretisation of the Masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ne=10000;
EE=linspace(min(V0),max(V0)+0.2,Ne);

EEmat   = repmat(EE'   ,[1 length(z)]);

V0mat   = repmat(V0 , [length(EE) 1]);
Egmat   = repmat(Eg , [length(EE) 1]);
Dsomat  = repmat(Dso, [length(EE) 1]);
Fmat    = repmat(F  , [length(EE) 1]);
  
EP=EP_K;
EPmat   = repmat(EP , [length(EE) 1]);
meK2  = 1 ./ ( 1 +  EPmat ./ (EEmat + Egmat-V0mat+Dsomat/3)  ); % m(z,E) in the 2 bands Kane s model
meK3  = 1 ./ ( 1 + 2/3 * EPmat ./ (EEmat+Egmat-V0mat) + 1/3*EPmat./(EEmat+Egmat+Dsomat-V0mat) ); % m(z,E) in the 3 bands Kane s model

EP=EP_L;
EPmat   = repmat(EP , [length(EE) 1]);
%meL = 1./(1 + 2*Fmat +  EPmat.*(EEmat+Egmat-V0mat + 2*Dsomat/3)./((EEmat+Egmat-V0mat).*(EEmat+Egmat-V0mat + Dsomat))  );  % m(z,E) in the 8 bands Luttinger s model
meL = 1./(1 + 2*Fmat +  2/3 * EPmat ./ (EEmat+Egmat-V0mat) + 1/3*EPmat./(EEmat+Egmat+Dsomat-V0mat) );  % m(z,E) in the 8 bands Luttinger s model

mhhL= 1 ./ (g1-2*g2);
mlhL= 1 ./ (g1+2*g2);


if PlotMass == 1
  figure('position',[10 100 800 500],'color','w');
  subplot(1,1,1,'fontsize',15)
  hold on;grid on;box on;
  idx=round(length(z)/2);
  %idx=1;
  plot(EE,meK2(:,idx),'b-')
  plot(EE,meK3(:,idx),'g-')
  plot(EE,meL(:,idx) ,'r-')
  set(gca,'ytick',0:0.005:0.2)
  legend('Kane 2bands','Kane 3bands','Luttinger')
  xlabel('Energy (eV)')
  ylabel('meff')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Selection of the model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E=[];
j=0;FS=15;FSt=10;
s{1}=strcat('\fontsize{',num2str(FSt),'}\color{blue}Potential');
s{2}=strcat('\fontsize{',num2str(FSt),'}\color{blue}------------');

col=[
1 0 0
0 1 0
%0 0 1
0 1 1
1 0 1
1 1 0
0 0 0
0.5 0.5 0.5
1 0.5 1
0 0.8 1
0.8 0 1
0 1 0.8
1 0.5 0.3
0.8 0 0.3
];

if FDM_kp1bandCB==1
    tic
    j=j+1;
    Mass=meK2(1,:);  % m=m(z)
    %Mass=meK3(1,:); % m=m(z)
    %Mass=meL(1,:);  % m=m(z)
    [E{j},psi{j}] = Schrod_1band_f(z,V0+DCBO,Mass,n);  % m = m(z)
    s{j+2}=strcat('\fontsize{',num2str(FSt),'}\color[rgb]{',num2str(col(j,:)),'}-> FDM kp1bandCB: m=m(z)');
    display(strcat('-> FDM_kp1bandCB =',num2str(toc),'sec'))
end

if FDM_kp1bandHH==1
    tic
    j=j+1;
    Mass=mhhL(1,:);  % m=m(z)  % the HH are all the time parabolic in ZB-001, even with strain!
    [E{j},psi{j}] = Schrod_1band_f(z,-(V0-Eg+DVBOHH),Mass,n);  % m = m(z)
    E{j}=-E{j};
    s{j+2}=strcat('\fontsize{',num2str(FSt),'}\color[rgb]{',num2str(col(j,:)),'}-> FDM kp1bandHH: m=m(z)');
    display(strcat('-> FDM_kp1bandHH =',num2str(toc),'sec'))
end

if FDM_kp1bandLH==1
    tic
    j=j+1;
    Mass=mlhL(1,:);  % m=m(z)% the LH can be almost parabolic if Dso is large and without strain!
    [E{j},psi{j}] = Schrod_1band_f(z,-(V0-Eg+DVBOLH),Mass,n);  % m = m(z)
    E{j}=-E{j};
    s{j+2}=strcat('\fontsize{',num2str(FSt),'}\color[rgb]{',num2str(col(j,:)),'}-> FDM kp1bandLH: m=m(z)');
    display(strcat('-> FDM_kp1bandLH =',num2str(toc),'sec'))
end

if FDM_kp2bands_Kane==1
    tic
    j=j+1;
    [E{j},psi{j}] = Schrod_2bands_Kane_f(z,V0,Eg,EP_K,Dso,n,ac,av,bv,exx,ezz);
    s{j+2}=strcat('\fontsize{',num2str(FSt),'}\color[rgb]{',num2str(col(j,:)),'}-> FDM kp2bands Kane: m=m(z,E)');
    display(strcat('-> FDM_kp2bands_Kane =',num2str(toc),'sec'))
end

if FDM_kp3bands_Kane==1
    tic
    j=j+1;
    [E{j},psi{j}] = Schrod_3bands_Kane_f(z,V0,Eg,EP_K,Dso,n,ac,av,bv,exx,ezz);
    s{j+2}=strcat('\fontsize{',num2str(FSt),'}\color[rgb]{',num2str(col(j,:)),'}-> FDM kp3bands Kane: m=m(z,E)');
    display(strcat('-> FDM_kp3bands_Kane =',num2str(toc),'sec'))
end

if FDM_kp6bands_Luttinger==1
    tic
    j=j+1;
    Mass=mhhL(1,:);  % m=m(z)
    [E{j},psi{j}] = Schrod_1band_f(z,-(V0-Eg+DVBOHH),Mass,n);  % m = m(z) the HH are all the time parabolic in ZB-001, even with strain!
    E{j}=-E{j};
    s{j+2}=strcat('\fontsize{',num2str(FSt),'}\color[rgb]{',num2str(col(j,:)),'}-> FDM HH 1band: m=m(z)');
    j=j+1;
    [E{j},psi{j}] = Schrod_2bands_Luttinger_Kohn_f(z,V0,Eg,Dso,g1,g2,g3,n,av,bv,exx,ezz);
    s{j+2}=strcat('\fontsize{',num2str(FSt),'}\color[rgb]{',num2str(col(j,:)),'}-> FDM LH/SO kp6bands Luttinger: m=m(z,E)');
    display(strcat('-> FDM_kp6bands_Luttinger =',num2str(toc),'sec'))
end

if FDM_kp8bands_Luttinger==1
    tic
    j=j+1;
    [E{j},psi{j},Ev,psi_v] = Schrod_3bands_Luttinger_Kohn_f(z,V0,Eg,EP_L,Dso,F,g1,g2,g3,n);
    s{j+2}=strcat('\fontsize{',num2str(FSt),'}\color[rgb]{',num2str(col(j,:)),'}-> FDM kp8bands Luttinger: m=m(z,E)');
    j=j+1;
    E{j}=Ev;psi{j}=psi_v;
    s{j+2}=strcat('\fontsize{',num2str(FSt),'}\color[rgb]{',num2str(col(j,:)),'}-> FDM kp8bands Luttinger: m=m(z,E)');
    display(strcat('-> FDM kp8bands_Luttinger =',num2str(toc),'sec'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ShootingCB_1band_Kane==1
    dE=0.02;
    precision=1e-6;
    tic
    j=j+1;
    Mass=meK2(1,:);  % m=m(z)
    s{j+2}=strcat('\fontsize{',num2str(FSt),'}\color[rgb]{',num2str(col(j,:)),'}-> ShootingCB 1band Kane: m=m(z)');
    [E{j},psi{j}] = Schrod_1band_shoot_f(z,V0+DCBO,Mass,n,dE,precision);
    display(strcat('-> ShootingCB_1band_Kane =',num2str(toc),'sec'))      
end

if ShootingHH_1band==1
    dE=0.02;
    precision=1e-5;
    tic
    j=j+1;
    Mass=mhhL; % m=m(z)   % the HH are all the time parabolic in ZB-001, even with strain!
    s{j+2}=strcat('\fontsize{',num2str(FSt),'}\color[rgb]{',num2str(col(j,:)),'}-> ShootingHH 1band: m=m(z)');
    [E{j},psi{j}] = Schrod_1band_shoot_f(z,-(V0-Eg+DVBOHH),Mass,n,dE,precision);
    E{j}=-E{j};
    display(strcat('-> SShootingHH_1band =',num2str(toc),'sec'))      
end

if ShootingCB_2bands_Kane==1
    dE=0.02;
    precision=1e-5;
    tic
    j=j+1;
    Mass=meK2; % m=m(z,E)
    s{j+2}=strcat('\fontsize{',num2str(FSt),'}\color[rgb]{',num2str(col(j,:)),'}-> ShootingCB 2bands Kane: m=m(z,E)');
    [E{j},psi{j}] = Schrod_Nbands_shoot_f(z,V0+DCBO,Mass,n,EE,dE,precision);
    display(strcat('-> ShootingCB_2bands_Kane =',num2str(toc),'sec'))      
end

if ShootingCB_3bands_Kane==1
    dE=0.02;
    precision=1e-5;
    tic
    j=j+1;
    Mass=meK3; % m=m(z,E)
    s{j+2}=strcat('\fontsize{',num2str(FSt),'}\color[rgb]{',num2str(col(j,:)),'}-> ShootingCB 3bands Kane: m=m(z,E)');
    [E{j},psi{j}] = Schrod_Nbands_shoot_f(z,V0+DCBO,Mass,n,EE,dE,precision);
    display(strcat('-> ShootingCB_3bands_Kane =',num2str(toc),'sec'))      
end

if ShootingCB_Luttinger==1
    dE=0.02;
    precision=1e-5;
    tic
    j=j+1;
    Mass=meL; % m=m(z,E)
    s{j+2}=strcat('\fontsize{',num2str(FSt),'}\color[rgb]{',num2str(col(j,:)),'}-> ShootingCB Luttinger: m=m(z,E)');
    [E{j},psi{j}] = Schrod_Nbands_shoot_f(z,V0+DCBO,Mass,n,EE,dE,precision);
    display(strcat('-> ShootingCB_Luttinger =',num2str(toc),'sec'))      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:length(E)
  for i=1:length(E{k})
      PSI{k}(:,i)=abs(psi{k}(:,i)).^2/max(abs(psi{k}(:,i)).^2)*ScF + E{k}(i); % normalisation for the plotting
  end
end

if PrintResultsISB == 1

for k=1:length(E)
  for i=1:length(E{k})
    for j=1:length(E{k})
      if j>i
        z_dipole{k}(i,j) = abs(  trapz( z , psi{k}(:,i).*z'.*psi{k}(:,j) )  );
        f_dipole{k}(i,j) = 2*m0/hbar^2 * ( E{k}(j)-E{k}(i) )* e * z_dipole{k}(i,j)^2 ;
        % Take care! Some people use meff inside the oscillator strenght f
        % Actually, meff has sens in an infinite QW because there is a single mass value
        % but not in multi-QW structure with various materials
        % https://www.nextnano.com/nextnano3/tutorial/1Dtutorial_IntrabandTransitions.htm
        % https://www.nextnano.com/nextnano3/tutorial/1Dtutorial_InGaAs_MQWs.htm
      end
    end
  end
end

display('')
display('===================================================')
display('Intersubband Results:')
display('===================================================')

for k=1:length(E)

display('')
display(s{k+2}(34:end))

if k>1
  %display('===================================================')
end

for i=1:length(E{k})
  for j=1:length(E{k})
    if j>i
      
      if E{k}(i)>0
        display(strcat(...
        'e',num2str(i),'-e',num2str(j),' = ',num2str( E{k}(j)-E{k}(i),'%.3f' ),'eV;   z'...
        ,num2str(i),'-',num2str(j),' = ',num2str( z_dipole{k}(i,j)*1e9,'%.3f' ),'nm;   f'...
        ,num2str(i),'-',num2str(j),' = ',num2str( f_dipole{k}(i,j),'%.3f' ) ...
        )  )
      else
        display(strcat(...
        'h',num2str(i),'-h',num2str(j),' = ',num2str( E{k}(j)-E{k}(i),'%.3f' ),'eV;   z'...
        ,num2str(i),'-',num2str(j),' = ',num2str( z_dipole{k}(i,j)*1e9,'%.3f' ),'nm;   f'...
        ,num2str(i),'-',num2str(j),' = ',num2str( f_dipole{k}(i,j),'%.3f' ) ...
        )  )
      end
    
    end
  end
end
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure('position',[-3500 100 1000 700],'color','w');
figure('position',[10 100 1000 700],'color','w');
subplot(1,1,1,'fontsize',15)
hold on;grid on;box on;

xscale=[z(1) z(end)]*1e9;
if PlotVB == 0
  yscale=[min(V0)-0.1 max(V0)+0.1];
else
  yscale=[min(V0-Eg)-0.1 max(V0)+0.1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shift=0;%Eg(round(length(z)/2));

plot(z*1e9,V0       +shift, 'b-','linewidth',2)
plot(z*1e9,V0-Eg    +shift, 'b-' ,'linewidth',2)
plot(z*1e9,V0-Eg-Dso+shift, 'b-','linewidth',2)
  
if StrainModel == 1
  plot(z*1e9,V0   +DCBO      +shift, 'b--','linewidth',2)
  plot(z*1e9,V0-Eg+DVBOHH    +shift, 'b--' ,'linewidth',2)
  plot(z*1e9,V0-Eg+DVBOLH    +shift, 'k--' ,'linewidth',2)
  plot(z*1e9,V0-Eg-Dso+DVBOSO+shift, 'm--' ,'linewidth',2)
end

for j=1:length(E)
  for i=1:length(E{j})
    
      plot(z*1e9,PSI{j}(:,i)+shift,'color',col(j,:),'linewidth',1)
      
  end
end

xlabel('z (nm)');
ylabel('Energy (eV)');

if StrainModel == 1
  title(strcat('T=',num2str(T),'K ; with STRAIN'))
else
  title(strcat('T=',num2str(T),'K ; without STRAIN'))
end

xlim(xscale)
ylim(yscale)

text(-2,(yscale(2)-yscale(1))*1.0+yscale(1),s,'background','w','edgecolor','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%