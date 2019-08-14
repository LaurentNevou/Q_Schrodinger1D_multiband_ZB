%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Layers Structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first column is the material used from the "library"
% second column is the length of the layer in nm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%substrate=GaAs;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
%M=[
%GaAs      10
%InGaAs20  10
%GaAs      10
%];

%substrate=GaAs;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
%M=[
%AlGaAs40      5
%GaAs          5
%AlGaAs40      5
%];

% Sirtori, PRB, 50, 8663 (1994)
substrate=InP;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
M=[
AlInAs   6
InGaAs   5.2
AlInAs   6
];

% Sirtori, PRB, 50, 8663 (1994)
%substrate=InP;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
%M=[
%AlInAs   6
%InGaAs   5.9
%AlInAs   1.3
%InGaAs   2.4
%AlInAs   6
%];

% Sirtori, PRB, 50, 8663 (1994)
%substrate=InP;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
%M=[
%AlInAs   6
%InGaAs   4.6
%AlInAs   1.0
%InGaAs   2.0
%AlInAs   1.0
%InGaAs   1.9
%AlInAs   6
%];

%substrate=InP;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
%M=[
%GaAsSb   10
%InGaAs   10
%GaAsSb   10
%];

%substrate=GaAs;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
%M=[
%GaAs      10
%InGaAs30  10
%GaAs      10
%InGaAs20  10
%GaAs      10
%InGaAs10  10
%GaAs      10
%InGaAs05  10
%GaAs      10
%];


%M=[
%AlAs      3
%GaAs          4.8
%AlAs      3
%];


%https://www.nextnano.de/nextnano3/tutorial/1Dtutorial12.htm
%substrate=GaAs;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
%M=[
%GaAs      10
%InGaAs30  7
%GaAs      10
%];

% M. Schneirt PhD thesis: "Optical Pumping: A Possible Approach towards a SiGe Quantum Cascade Laser", page 24
%https://core.ac.uk/download/pdf/20642119.pdf
%substrate=SiGe50;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
%M=[
%SiGe50  10
%Si      2
%SiGe50  3.5
%Si      2
%SiGe50  10
%];

% M. Schneirt PhD thesis: "Optical Pumping: A Possible Approach towards a SiGe Quantum Cascade Laser", page 86
%https://core.ac.uk/download/pdf/20642119.pdf
%substrate=SiGe25;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
%M=[
%Si      5
%SiGe50  3.6
%Ge      1.2
%Si      5
%];

% M. El-Kurdi PhD thesis: "Dispositifs a ilots de GeSi pour la microphotonique proche infrarouge sur silicium", page 54
%substrate=Si;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
%M=[
%Si      5
%SiGe50  5
%Si      5
%];

% Paul Harrison book: 
% Quantum Wells, Wires and Dots.
% 4th edition (2016),
% Chap14: "Multiband envelope function (k.p) method", page 503
%substrate=SiGe30;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
%M=[
%Si      5
%SiGe40  11.1
%Si      5
%];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
