%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Schrodinger solver on uniform grid with m(z,E)!!! %%%%%%%%%%%%%%%
%%%%% With the non-parabolic band 3x3k.p Luttinger model for Ec, LH and SO %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Ec,psi_c,Ev,psi_v]=Schrod_3bands_Luttinger_Kohn_f(z,Vc,Eg,EP,Dso,F,g1,g2,g3,n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant J.s
hbar=h/(2*pi);
e=1.602176487E-19;              %% charge de l electron Coulomb
m0=9.10938188E-31;              %% electron mass kg
h0=hbar^2/(2*m0) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nz=length(z);
dz = z(2)-z(1);

shift=min(Vc);
Vc =Vc-shift;
Vv =Vc-Eg;
Vso=Vc-Dso-Eg;

Vv(1)=-10;
Vv(end)=-10;
Vso(1)=-10;
Vso(end)=-10;

Vc  = [ (Vc(1:end-1)    + Vc(2:end)) / 2     Vc(end)  ];
Vv  = [ (Vv(1:end-1)    + Vv(2:end)) / 2     Vv(end)  ];
Vso = [ (Vso(1:end-1)   + Vso(2:end))/ 2     Vso(end) ];

gc = 1+2*F;
g1 = g1-EP/(3*Eg);
g2 = g2-EP/(6*Eg);
g3 = g3-EP/(6*Eg);

%gc=gc(round(Nz/2))*ones(1,Nz);
%g1=g1(round(Nz/2))*ones(1,Nz);
%g2=g2(round(Nz/2))*ones(1,Nz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Building of the operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First derivative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DZ1c = (0.5)*diag(ones(1,Nz-1),+1) + (-0.5)*diag(ones(1,Nz-1),-1)  ;
DZ1b =   (1)*diag(ones(1,Nz)  ,0 ) +   (-1)*diag(ones(1,Nz-1),-1)  ;
DZ1f =   (1)*diag(ones(1,Nz-1),+1) +   (-1)*diag(ones(1,Nz)  ,0 )  ;

DZ1c=DZ1c/dz;
DZ1b=DZ1b/dz;
DZ1f=DZ1f/dz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Building of the operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear DZ2 a b c

g1p = [         (g1(1:end-1) + g1(2:end)) / 2   g1(end) ];
g1m = [   g1(1) (g1(1:end-1) + g1(2:end)) / 2           ];

b = (g1p + g1m) .* ones(1,Nz) ;
a =  g1m(2:end) .* ones(1,Nz-1) ;
c =  g1m(2:end) .* ones(1,Nz-1) ;

DZ2 = (-1)*diag(b)  +  (1)*diag(a,-1)  +  (1)*diag(c,+1) ;
DZ2 = DZ2 / dz^2;

P = + h0 * DZ2 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear DZ2 a b c

g2p = [         (g2(1:end-1) + g2(2:end)) / 2   g2(end) ];
g2m = [   g2(1) (g2(1:end-1) + g2(2:end)) / 2           ];

b = (g2p + g2m) .* ones(1,Nz) ;
a =  g2m(2:end) .* ones(1,Nz-1) ;
c =  g2m(2:end) .* ones(1,Nz-1) ;

DZ2 = (-1)*diag(b)  +  (1)*diag(a,-1)  +  (1)*diag(c,+1) ;
DZ2 = DZ2 / dz^2;

Q = -2*h0 * DZ2 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear DZ2 a b c

gcp = [         (gc(1:end-1) + gc(2:end)) / 2   gc(end) ];
gcm = [   gc(1) (gc(1:end-1) + gc(2:end)) / 2           ];

b = (gcp + gcm) .* ones(1,Nz) ;
a =  gcm(2:end) .* ones(1,Nz-1) ;
c =  gcm(2:end) .* ones(1,Nz-1) ;

DZ2 = (-1)*diag(b)  +  (1)*diag(a,-1)  +  (1)*diag(c,+1) ;
DZ2 = DZ2 / dz^2;

A = + h0 * DZ2 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chineese method s
%
% Xunpeng Ma, Kangwen Li, Zuyin Zhang, Haifeng Hu, Qing Wang, Xin Wei,and Guofeng Song
% "Two-band finite difference method for the bandstructure calculation with nonparabolicity effects in quantum cascade lasers"
% JAP, 114, 063101 (2013)

%EP=EP(round(Nz/2));
%EP=EP*0.7;
H00 =  zeros(Nz,Nz);

U12 = +sqrt(1/3) * hbar/sqrt(2*m0) * (  diag(sqrt(EP*e),0) + diag(sqrt(EP(1:end-1)*e),-1)  ) .* DZ1b ;
U13 = U12;
U21 = -sqrt(1/3) * hbar/sqrt(2*m0) * (  diag(sqrt(EP*e),0) + diag(sqrt(EP(1:end-1)*e),+1)  ) .* DZ1f ;
U31 = U21;

%U12 = H00 ;
%U13 = H00 ;
%U21 = H00 ;
%U31 = H00 ;


%U =  sqrt(1/3) * P0 * kz;
%U =  sqrt(1/3) * hbar * sqrt(EP*e/2/m0) .* DZ1c ;
%V =  H00;
%W =  H00;
%T =  H00;
%R =  H00;
%S =  H00;

H11 = A + diag(Vc*e);
%H22 = A + diag(Vc*e);
%H33 = P+Q + diag(Vv*e);
H44 = P-Q + diag(Vv*e);
%H55 = P-Q + diag(Vv*e);
%H66 = P+Q + diag(Vv*e);
H77 = P   + diag(Vso*e);
%H88 = P   + diag(Vso*e);

%Hdiag = blkdiag( H11 , H22 , H33 , H44 , H55 , H66 , H77 , H88 );
%
%Hdiag = [A A P+Q P-Q P-Q P+Q -Dso+P -Dso+P];
% Ec- Ec+     HH+               LH+             LH-            HH-            SO-            SO+
%H=[
%  0   0        0               T'+V'       sqrt(2)*(W-U)  -sqrt(3)*(T-V)     (W-U)     sqrt(2)*(T'+V')  % Ec-
%  0   0  -sqrt(3)*(T'+V')  sqrt(2)*(W-U)       (T-V)            0       -sqrt(2)*(T-V)      W'+U        % Ec+
%  0   0        0                -S               R              0       -sqrt(2)  *R     sqrt(1/2)*S    % HH+
%  0   0        0                 0               0              R        sqrt(3/2)*S     sqrt(2)  *Q    % LH+
%  0   0        0                 0               0              S       -sqrt(2)  *Q     sqrt(3/2)*S'   % LH-
%  0   0        0                 0               0              0        sqrt(1/2)*S'    sqrt(2)  *R'   % HH-
%  0   0        0                 0               0              0              0              0         % SO-
%  0   0        0                 0               0              0              0              0         % SO+
%];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hdiag = blkdiag( H11 , H44 ,  H77 );

H=[
      H00            -sqrt(2)*U12           -U13       % Ec+
 -sqrt(2)*U21            H00             sqrt(2)*Q     % LH+
     -U31             sqrt(2)*Q'             H00       % SO-
];

H =  + H + Hdiag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Diagonalisation of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H3x3=sparse(H);
[psi_3x3,Energy] = eigs(H3x3,n,'SM');
%[psi_3x3,Energy] = eig(H);
E_3x3 = diag(Energy)/e ;

psi_c=[];
Ec=[];
psi_v=[];
Ev=[];

for i=1:n%length(H(1,:))
  if E_3x3(i) > min(Vc)
    Ec=[Ec E_3x3(i)]; 
    psi_c = [psi_c  psi_3x3(1:length(z),i) ];
  end
  
  if E_3x3(i) < max(Vv)
    E_3x3(i);
    Ev=[Ev E_3x3(i)];
%    psi_v = [psi_v  psi_3x3(:,i) ];
%    psi_v = [psi_v psi_3x3(1:Nz,i) + psi_3x3(Nz+1:2*Nz,i) + psi_3x3(2*Nz+1:end,i)]; % HERE, I am not sure it is correct!!!
    psi_v = [psi_v  psi_3x3(Nz+1:2*Nz,i) + psi_3x3(2*Nz+1:end,i)]; % HERE, I am not sure it is correct!!!
  end
end

Ec=Ec'+shift;
Ev=Ev'+shift;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Normalization of the Wavefunction %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(Ec)
    psi_c(:,i)=psi_c(:,i)/sqrt(trapz(z',abs(psi_c(:,i)).^2)); 
end

for i=1:length(Ev)
%    psi_v(:,i)=psi_v(:,i)/sqrt(trapz(z',abs(psi_v(:,i)).^2));
%    psi_v(:,i)=psi_v(:,i)/sqrt(trapz([z z z z]',abs(psi_v(:,i)).^2));

%size(psi_3x3(1:Nz,i))
%size(psi_3x3(Nz+1:2*Nz,i))
%size(psi_3x3(2*Nz+1:end,i))

%  psi_v(:,i)= psi_3x3(1:Nz,i) + psi_3x3(Nz+1:2*Nz,i) + psi_3x3(2*Nz+1:end,i); % HERE, I am not sure it is correct!!!
  psi_v(:,i)= psi_v(:,i)/sqrt(trapz(z',abs(psi_v(:,i)).^2)); 
  %psi_v(:,i)= psi_v(:,i)/sqrt(trapz([z z z z z z]',abs(psi_v(:,i)).^2)); 
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

if length(Ev)>1
if Ev(1)<Ev(2)
  psi_v=psi_v(:,end:-1:1);
  Ev=Ev(end:-1:1);
end
end

end