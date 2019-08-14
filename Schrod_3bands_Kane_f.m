%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Schrodinger solver on uniform grid with m(z,E)!!! %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% With the non-parabolic band 3x3k.p Kane model %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Ec,psi_c]=Schrod_3bands_Kane_f(z,Vc,Eg,EP,Dso,n,ac,av,bv,exx,ezz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant J.s
hbar=h/(2*pi);
e=1.602176487E-19;              %% charge de l electron Coulomb
m0=9.10938188E-31;              %% electron mass kg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nz=length(z);
dz = z(2)-z(1);

eyy = exx;
DCBO   = -abs(ac).*(exx+eyy+ezz) ; % shift of the CB due to strain
DVBOLH = +abs(av).*(exx+eyy+ezz) + abs(bv).*(exx-ezz) ; % shift of the VB due to strain
DVBOSO = +abs(av).*(exx+eyy+ezz) ; % shift of the VB due to strain

Vc=Vc+DCBO;
Vc(1)=5;
Vc(end)=5;

shift=min(Vc);
Vc=Vc-shift;

Vv=Vc-Eg+DVBOLH;
Vso=Vc-Dso-Eg+DVBOSO;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Building of the operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First derivative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DZ1c = (0.5)*diag(ones(1,Nz-1),+1) + (-0.5)*diag(ones(1,Nz-1),-1)  ;
DZ1b =   (1)*diag(ones(1,Nz)  ,0 ) +   (-1)*diag(ones(1,Nz-1),-1)  ;
DZ1f =   (1)*diag(ones(1,Nz-1),+1) +   (-1)*diag(ones(1,Nz)  ,0 )  ;

%DZ1c=DZ1c/dz;
DZ1b=DZ1b/dz;
DZ1f=DZ1f/dz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second derivative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DZ2 =(-2)*diag(ones(1,Nz)) + (1)*diag(ones(1,Nz-1),-1) + (1)*diag(ones(1,Nz-1),1);
DZ2=DZ2/dz^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vc  = [ (Vc(1:end-1)    + Vc(2:end)) / 2     Vc(end)  ];
Vv  = [ (Vv(1:end-1)    + Vv(2:end)) / 2     Vv(end)  ];
Vso = [ (Vso(1:end-1)   + Vso(2:end))/ 2     Vso(end) ];
EP  = [ (EP(1:end-1)    + EP(2:end)) / 2     EP(end)  ];

H0=(-(hbar^2)/(2*m0)) *  DZ2   ;

H11 =  H0 + diag( Vc *e ) ;
H22 = -H0 + diag( Vv *e ) ;
H33 = -H0 + diag( Vso*e ) ;

% Xunpeng Ma et al. JAP, 114, 063101 (2013)
% "Two-band finite difference method for the bandstructure calculation with nonparabolicity effects in quantum cascade lasers"
% ==> It s seems to be by far the most accurate method

H12 = +sqrt(2/3)*hbar/sqrt(2*m0) * (  diag(sqrt(EP*e),0) + diag(sqrt(EP(1:end-1)*e),-1)  ) .* DZ1b ;
H13 = -sqrt(1/3)*hbar/sqrt(2*m0) * (  diag(sqrt(EP*e),0) + diag(sqrt(EP(1:end-1)*e),-1)  ) .* DZ1b ;

H21 = -sqrt(2/3)*hbar/sqrt(2*m0) * (  diag(sqrt(EP*e),0) + diag(sqrt(EP(1:end-1)*e),+1)  ) .* DZ1f ;
H31 = +sqrt(1/3)*hbar/sqrt(2*m0) * (  diag(sqrt(EP*e),0) + diag(sqrt(EP(1:end-1)*e),+1)  ) .* DZ1f ;

H23 = H11*0;
H32 = H11*0;

H3x3=[
H11  H12  H13
H21  H22  H23
H31  H32  H33
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Diagonalisation of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H3x3=sparse(H3x3);
[psi_3x3,Energy] = eigs(H3x3,n,'SM');
E_3x3 = diag(Energy)/e ;

psi_c=[];
Ec=[];

for i=1:n
  if E_3x3(i) > min(Vc)
    Ec=[Ec abs(E_3x3(i))];
    psi_c = [psi_c  psi_3x3(1:length(z),i) ];
  end
end

Ec=Ec'+shift;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Normalization of the Wavefunction %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(Ec)
    psi_c(:,i)=psi_c(:,i)/sqrt(trapz(z',abs(psi_c(:,i)).^2)); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here is a small patch due to differences between Octave and Matlab
% Matlab order the eigen values while Octave reverse it

if length(Ec)>1
if Ec(1)>Ec(2)
  psi_c=psi_c(:,end:-1:1);
  Ec=Ec(end:-1:1);
end
end

end