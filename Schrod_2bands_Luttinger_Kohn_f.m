%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Schrodinger solver on uniform grid with m(z,E)!!! %%%%%%%%%%%%%%%
%%%%%%%% With the non-parabolic band 2x2k.p Luttinger model for LH and SO %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is valid only if the growth direction is along z-axis
% In that case (and even with strain), HH is not coupled to the other band LH 
% and SO and therefore can be computed separatly
%
function[Ev,psi_v]=Schrod_2bands_Luttinger_Kohn_f(z,Vc,Eg,Dso,g1,g2,g3,n,av,bv,exx,ezz)

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

Vv =Vc-Eg;
shift=max(Vv);
Vv=Vv-shift;
Vso=Vv-Dso;

%Vv(1)=-10;
%Vv(end)=-10;
%Vso(1)=-10;
%Vso(end)=-10;

Vv  = [ (Vv(1:end-1)    + Vv(2:end)) / 2     Vv(end)  ];
Vso = [ (Vso(1:end-1)   + Vso(2:end))/ 2     Vso(end) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Building of the operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eyy = exx;
%exy = 0;
%ezx = 0;
%eyz = 0;
%ee  = exx+eyy+ezz;
av  = abs(av)*e;
bv  = abs(bv)*e;
%dv  = abs(dv)*e;

Pe =  diag( +av .* (exx+eyy+ezz) ) ;
Qe =  diag( -bv .* (exx-ezz)     ) ;
%Re =  sqrt(3)/2 * bv*(exx-eyy) - 1i*dv*exy;
%Se = -dv * (ezx - 1i*eyz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g1p = [         (g1(1:end-1) + g1(2:end)) / 2   g1(end) ];
g1m = [   g1(1) (g1(1:end-1) + g1(2:end)) / 2           ];

b = (g1p + g1m) .* ones(1,Nz) ;
a =  g1m(2:end) .* ones(1,Nz-1) ;
c =  g1m(2:end) .* ones(1,Nz-1) ;

DZ2 = (-1)*diag(b)  +  (1)*diag(a,-1)  +  (1)*diag(c,+1) ;
DZ2 = DZ2 / dz^2;

P = + h0 * DZ2 + Pe ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g2p = [         (g2(1:end-1) + g2(2:end)) / 2   g2(end) ];
g2m = [   g2(1) (g2(1:end-1) + g2(2:end)) / 2           ];

b = (g2p + g2m) .* ones(1,Nz) ;
a =  g2m(2:end) .* ones(1,Nz-1) ;
c =  g2m(2:end) .* ones(1,Nz-1) ;

DZ2 = (-1)*diag(b)  +  (1)*diag(a,-1)  +  (1)*diag(c,+1) ;
DZ2 = DZ2 / dz^2;

Q = -2*h0 * DZ2 + Qe ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H00 =  zeros(Nz,Nz);

%R =  H00;
%S =  H00;

%H11 = P+Q + diag(Vv*e);
H22 = P-Q + diag(Vv*e)  ;
%H33 = P-Q + diag(Vv*e);
%H44 = P+Q + diag(Vv*e);
H55 = P   + diag(Vso*e) ;
%H66 = P   + diag(Vso*e);


%P = -H0 * g1 * k^2;
%Q = -H0 * g2 *(kx^2 + ky^2 - 2*kz^2);
%R =  H0 * sqrt(3)  * (g2*(kx^2-ky^2) - 2i*g3*kx*ky );
%S =  H0 *2*sqrt(3) * g3*(kx-1i*ky)*kz;

%Hdiag = blkdiag( H11 , H22 , H33 , H44 , H55 , H66 );
%
%%   HH     LH       LH        HH           SO              SO
%
%H=[
%   H00     -S        R        H00    -sqrt(1/2)*S       sqrt(2)  *R  % HH
%   H00     H00      H00        R     -sqrt(2)  *Q       sqrt(3/2)*S  % LH
%   H00     H00      H00        S      sqrt(3/2)*S'      sqrt(2)  *Q  % LH
%   H00     H00      H00       H00    -sqrt(2)  *R'     -sqrt(1/2)*S' % HH
%   H00     H00      H00       H00         H00              H00        % SO
%   H00     H00      H00       H00         H00              H00        % SO
%];


Hdiag = blkdiag( H22  , H55  );

%  LH         SO
H=[
   H00       sqrt(2)  *Q  % LH
   H00        H00         % SO
];


H=H'+H+Hdiag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Diagonalisation of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H2x2=sparse(H);
[psi_2x2,Energy] = eigs(H2x2,n,'SM');
Ev = diag(Energy)/e ;
Ev=Ev+shift;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Normalization of the Wavefunction %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(Ev)
%    psi_v(:,i)=psi_v(:,i)/sqrt(trapz(z',abs(psi_v(:,i)).^2));
%    psi_v(:,i)=psi_v(:,i)/sqrt(trapz([z z z z]',abs(psi_v(:,i)).^2));

  psi_v(:,i)= psi_2x2(1:Nz,i) + psi_2x2(Nz+1:end,i); % HERE, I am not sure it is correct!!!
  psi_v(:,i)= psi_v(:,i)/sqrt(trapz(z',abs(psi_v(:,i)).^2)); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here is a small patch due to differences between Octave and Matlab
% Matlab order the eigen values while Octave reverse it

if length(Ev)>1
if Ev(1)<Ev(2)
  psi_v=psi_v(:,end:-1:1);
  Ev=Ev(end:-1:1);
end
end

end