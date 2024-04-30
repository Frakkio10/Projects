%%%%%%%%%%%%%%%%%%
%%%   BEAM3D   %%%
%%%%%%%%%%%%%%%%%%

function [varargout]=beam2010(varargin)

%fprintf(' Trace de rayon @YM,LC,STF \n')
%fprintf(' Trace de rayon 2D pour micro-ondes Yannick MICHELOT (1993) \n')
%fprintf(' Adaptation Difpol : Laurent COLAS (1994/95), G.U.I. :  SEAK Teng-Fong (1998/99)\n')
%fprintf(' Approximation de l''Optique G?om?trique, faisceau de rayons
%d?pandents\n')
%fprintf(' Equations de Hamilton : dx/dt = -(dD/dk)/(dD/dw) dk/dt = +(dD/dx)/(dD/dw)\n')
% Auteur du ray tracing et des fonctions ci-dessous : Y. MICHELOT 1993
% -> plasma froid, shift de shafranov, profils th?oriques
% Modifications : Y. MICHELOT & L.COLAS
% -> 06/05/94 : profils exp?rimentaux de densit?
% -> 21/06/94 : Propagation  dans le vide, nombre de rayons quelconques
% -> 08/07/94 : Franchissement de l'interface plasma-vide : descartes
% -> 09/95 : Simulation de chocs fictifs, p?n?tration maximale des rayons 
% -> 11/95 :  d(psi)=dsha0(1-psi^(dshap/2)) ; dshap quelconque, experimental
% -> corrections sur B : paramagn?tisme / diamagn?tisme / ripple / B_polo?dal
% Modifications : SEAK Teng-Fong  1999-02-02
% -> Changements trop importants pour tout citer, e.g. sim_decor.m
% Modifications : Cyrille Honor? 2002
% -> Adaptation pour la reflectometrie Doppler
% -> 3D avec geometrie Ripple
% -> Trace de faisceau gaussien
% PH input exp & output graph
% Modifications : PH 2008 separation calcul N2 etc...
% varargin: 1: diagnostic configuration diag
% varargin: 2: beam frequency frei
%           3: beam angle pol
%           4: modex : O mode / X mode
%           5: beam angle tor
%           6: structure for ne profile
%           7: structure for plasma equilibrium and magnetic field
clear global

global rr0 a0r zpos el tau dsha0 dshap
global B0 RR j0 jp Te0 Tepiq
global phigola ripple
global rhobord nx2test
global ne_spline
global modex w
global imageik nbrayr nbraya
global rko hmin hmax tol testklim
global debog

% Constantes physiques
c = 3e8;  % [m/s]
epsi0 = 8.8541e-12;
me = 9.1e-31;
qe = 1.6e-19;

% Deboggage
debog = 0;  % Constante de deboggage (0=non)

if ~exist('username','var'),username=getenv('LOGNAME');end;
dateplot=date;

% initialisation of parameters (a priori not redefined after)
% unless NaN ->  read by the function "entree" 
%%%%%%%%%%%%%%%%%%%%
% Display Parameters
figureprofil = 1; % Display density radial profile 
figurewaist = 0;  % Display graph of waist & wavefront curvature as function s s (1='yes')
figurepoloidale = 1;  % Display poloidal figure  (1='oui')
figuretoroidale = 1;  % Display toroidal figure (1='oui')
figureprint = 0;  % figure automatic print (1='oui')
figure3d = 0;  % Display of 3D figure (1='oui')
%afficheanalytique = 1;  % Display text for analytical results  (1) 
afficherayons = 0;  % analytique text for ray tracing results  (1)
figurearretrayon = 0;  % stops ray tracing at the reflexion (turning point) 
set(0,'DefaultTextFontSize',14)  % initialisations pour les figures
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultLineLineWidth',1)
set(0,'DefaultAxesLineWidth',1)
set(0,'DefaultLineMarkerSize',8)

% Ray tracing parameters
raydim = 3;  % 2d (2) ou 3d [3]

nbrayr = 2;  % nombre de rayons radialement >=1 [2 9]
nbraya = 4;  % nombre de rayons par tour 4m [4 16]

% Limite du plasma
rhobord = 1.25;  % plasma extreme extent i.e. ne=0 & nx2=1 beyond [1.1]
rhobord = entree(rhobord,'Rho limite du plasma (ne=0,nx2=1 au dela) ?',1.1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% varargin parameters
% if no argument: input each parameter for simulation 

%%%%%%%%%%%%%%%%%%%%
% Diagnostic Parameters : which diagnostic configuration? 
if nargin>0
    diag = varargin{1}; 
else
    diag = 'difdop';  % diagnostic : ['difdop'], 'dreve' 'drefluc' ou 'difpol'
    % Possibilit� d'utiliser beam pour d'autres diagnostics
    diag = lower(entree(diag,'"Difdop" ou "Difpol" ou "" ?',1,'s'));
end
% Affichage : Titre
titrechoc = upper(diag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diag parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beam frequency
if nargin>1
    frei = varargin{2};  
else
    frei =  entree(NaN,'frequence sonde (GHz) ?',60);% frequence de la source (GHz) [60]
end
lambda = c/frei/1e9;  % longueur d'onde dans le vide [m]
ko = 2*pi/lambda;  % nombre d'onde dans le vide [m^-1]
w = 2*pi*frei*1e9;  % pulsation de la source

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% difdop poloidal antenna angle
if nargin>2
    difdoptheta = varargin{3};  
else
    difdoptheta = entree(NaN,'angle antenne ?',5);%  
    %difdoptheta = NaN;%5;%   % visee poloidale (0?=vers centre plasma) [?] pour difdop -7.5
    %difdoptheta = entree(difdoptheta,'theta poloidal (0?=visee centre plasma) (?) ?',0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beam polarisation
% mode O (modex=0) ou mode X avec B_vide et ripple (1) ou mode X avec B_para/diamagn. (2)
if nargin>3
    modex = varargin{4};  % polarisation
else
    modex = entree(NaN,'Polarisation (O) ?',0);          %  polarisation mode 0 par d?fault
end
%    afficheanalytique = 0;  % pas d'affichage texte des resultats analytiques (1) 
if modex==0
    afficheanalytique = 1;  % Affichage texte des resultats analytiques (1)
else
    afficheanalytique = 0 ; %
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% difdop antenna toroidal orientation / torus axis
if nargin>4
    phidirgoladeg = varargin{5};  
else
    if raydim==3
        phidirgoladeg = entree(NaN,'orientation toroidale en phi de l''antenne ?',1.85);%
    else
        phidirgoladeg = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
% PLASMA PARAMETERS    %
%%%%%%%%%%%%%%%%%%%%%%%%
% if input arguments: experimental profile, otherwise simulation
if nargin>5
    simulation=0;
    ne_struct=varargin{6};
    nefit=ne_struct.ne;
    rhofit=ne_struct.rhofit;
    %nebord ?

else
    % Affichage : Titre
    simulation=1;
    titrechoc = sprintf('%s   |  Simulation ',titrechoc);
    fprintf('\n%s\n',titrechoc);
    % ask values for analytic
    ne0 = NaN;  %5.5;  % central density (10^19 m^-3)
    ne0 = entree(ne0,'densit? centrale (10^19 m^-3) ?',5);
    ne0 = ne0*1e19;
    ne10 = 0.15; %0.15;  % density at rho=1 / density at rho=0
    ne10 = entree(ne10,'rapport ne(1)/ne(0) ?',.15);
    nepiq = NaN;% .6;  % density profile peaking  [0.6]
    nepiq = entree(nepiq,'piquage profil densit� (nepiq=1 => profil parabolique) ?',.6);
    nelambdan = 2; % density decreasing length (exponential) (in the SOL) (cm)
    nelambdan = entree(nelambdan,'longueur de decroissance de ne dans la SOL (cm) ?',2);
    nelambda = nelambdan/70;  % normalisation au petit rayon

    % Display which appears on figure
    textpla = sprintf('ne(0)=%.1f 10^{19}m^{-3}\n',ne0/1e19);
    textpla = sprintf('%s piq_{ne}=%.1f \n ne_a/ne_0=%.2f\n',textpla,nepiq,ne10);
    textpla = sprintf('%s \\lambda_{ne}=%.0f cm ',textpla,nelambda*a0r);
    fprintf('%s\n',textpla)
    
    % Analytical density profile
    rhofit = [0:.05:.9 .94:.02:1.24 1.3 1.4 1.5 1.7 2 2.5 3];
    necoeur = ne0 * ( ne10 + (1-ne10) * (1-rhofit.^2).^nepiq ) .*(rhofit>=0).*(rhofit<=1);
    nesol = ne0*ne10 * exp(-(sqrt(rhofit)-1)/nelambda) .* (rhofit>1).*(rhofit<=rhobord);
    nefit = necoeur + nesol;clear necoeur nesol;
end
% Figure du profil de densite  
% aussi pour calculer ne_spline
ne_spline=nerhofit(rhofit,nefit);
if figureprofil
    rho = [0:.05:.9 .94:.04:1.24 1.3 1.4 1.5 1.7 2 2.5 3];%[0:.002:.95 .96:.001:1.25 1.26:.005:3];
    %        ne = ppval(ne_spline,rho');
    figure(100)
    plot(rhofit,nefit,'o',rho,nerho(rho))
    xlabel('rho')
    ylabel('ne')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% structure for equilibrium
if nargin>6
    plasma=varargin{7};
    a0r = plasma.a0r;  % minor radius  (cm) [70]
    rr0 = plasma.rr0;  % major radius (cm) [240]
    zpos = plasma.zpos;  % vertical position of the plasma center (cm) [0] (not used yet)
    el = plasma.e; % ellipticity [0]
    tau = plasma.tau; % triangularity [0]
    % Shafranov shift
    dsha0 = plasma.dsha0;
    dshap = plasma.dshap;
    % Magnetic field
    B0 = plasma.B0;
    if modex==2
        j0=plasma.j0;
        jp=plasma.jp;
        Ip=  plasma.Ip;
        qpsi = 4; % edge safey factor
        Te0 = 3;  % central electonic temperature
        Tepiq = 2;% .6;  % temperature profile peaking [2]
    end
    paroiz=plasma.paroiz;
    paroix=plasma.paroix;

else
    % minor and major radius
    a0r = NaN;  % minor radius  (cm) [70]
    rr0 = NaN;  % major radius (cm) [240]
    zpos = NaN;  % vertical position of the plasma center (cm) [0] 
    el = NaN;  % ellipticity
    tau = NaN;  %triangularity
    a0r = entree(a0r,'petit rayon (cm) ?',70);
    rr0 = entree(rr0,'grand rayon (cm) ?',240);
    zpos = entree(zpos,'position verticale du centre du plasma (cm) ?',0);
    el = entree(el,'ellipticity (vertical elongation =b/a) ?',1);
    tau = entree(tau,'triangularity (cm) ?',0);
    % Shafranov shift
    dsha0 = NaN;  %  Shafranov shift  (central, rho=0) (cm) [7]
    dsha0 = entree(dsha0,'Shafranov shift (cm) ?',7);
    dshaexp = 0;  % type de piquage de decalage de shafranov: experimental(1) ou dshap=2 (0)
    dshaexp = entree(dshaexp,'Type de piquage du Shafranov : dshap=2 (0), dshap issu des donnees (1)',0);
    dshap=2;
    % Magnetic field
    B0 = NaN;  %  magnetic field at toroidal coil center (i.e at RR=2.36) (T)
    B0 = entree(B0,'Champ magn?tique au centre (T)?',3.8);
    if modex==2
        disp('mode X avec corrections para/dia')
        qpsi = 4; % edge safey factor
        Ip = NaN; % Plasma current
        Ip = entree(Ip,'courant plasma ?',1);
        Te0 = 3;  % central electonic temperature
        Te0 = entree(Te0,'Te0 ?',3);
        Tepiq = 2;% .6;  % temperature profile peaking [2]
    end
end
RR = 236;   % toroidal coil center for magnetic field calculation (cm)

textEqui = sprintf('R=%.2f m, Z= %.0f, a=%.0f cm, e=b/a=%.2f,  tau=%.1f cm, B0=%.1f T\n',rr0/100,zpos,a0r,el,tau,B0);
textEqui = sprintf('%s dsha0=%.1f cm  piq_{dsha}=%.1f\n ',textEqui,dsha0,dshap);
%textEqui = sprintf('%s e=%.1f cm  tau=%.1f cm',textEqui,el,tau);
fprintf('%s\n',textEqui)

% 7th argument: different diag?

%%%%%%%%%%%%%%%%%%%%%%%%
ripple = 1;  % ripple taken into account (for Tore Supra geometry, and for  B field) ([1]=oui)
if raydim==3
  ripple = entree(ripple,'prise en compte du ripple (1) ou pas (0) ?',1);
else
  ripple = 0;
end
nx2test = [];  % valeurs d'indices imposees dans le plasma
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametres du diagnostic 

difdopdivcst = 1;  % divergence angulaire en 1/freq selon Millitech (0) ou independante de la frequence (1)
Rgola = [];  % position en R de l'antenne (cm)
zgola = [];  % position en z de l'antenne (cm) 
phigoladeg = [];  % position en phi de l'antenne (0=plan median entre bobines)
dwaist = 50;  % distance de l'antenne au waist
waist = 40;  % taille du faisceau (mm)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position initiale de l'antenne %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%thetadirgoladeg = [];  % orientation poloidale de l'antenne (?)

% Differents diagnostics
if strcmp(diag,'difdop')
    % reflectometrie doppler
    Rgola = 430;  % 430=240+70+120
    difdopdivcst = entree(difdopdivcst,'divergence constante (1) ou selon Millitech (0) ?',1);
    if difdopdivcst
      focale = 120;      % (mm)
      waistgola = 4.7;  % (mm)
      waist = focale*lambda*1000 /pi/waistgola;
      dwaist = focale/10;  % (cm)
    else
      waist = 39.7; % Millitech : dtheta = 826/(DF) ["][GHz]->[?]
      dwaist = 12;  % (cm)
    end
    phigoladeg = -3;% translation de l'antenne par rapport au plan median du queusot
    zgola = -27;
    %orientation poloidale de l'antenne (0?=visee centre plasma) !! centre
    %plasma - > centre de la chambre!!
    %thetadirgoladeg = angle(-(Rgola-rr0)-i*(zgola))*180/pi-difdoptheta;
    thetadirgoladeg = 180-7.5-difdoptheta;
    % Antenne
    dessingolax = waist/10*[   0,   0:-.6:-4.8,     -4.8:.6:0, .34,.54,.6,.54,.34,  0];
    dessingolar = waist/10*[-1.8,1.8*ones(1,9),-1.8*ones(1,9),-1.2,-.6, 0, .6,1.2,1.8];
elseif strcmp(diag,'dreve')
    % 3eme voie reflectometrie doppler visee haute
    Rgola = 240+15;  % queusot haut
    % meme antenne difdop
    difdopdivcst = entree(difdopdivcst,'divergence constante (1) ou selon Millitech (0) ?',1);
    if difdopdivcst
      focale = 120;      % (mm)
      waistgola = 4.7;  % (mm)
      waist = focale*lambda*1000 /pi/waistgola;
      dwaist = focale/10;  % (cm)
    else
      waist = 39.7; % Millitech : dtheta = 826/(DF) ["][GHz]->[?]
      dwaist = 12;  % (cm)
    end
    phigoladeg = 0;% dans le plan median du queusot
    zgola = 220;
    %orientation poloidale (absolue) de l'antenne 
    %fixe, 2� vers le plasma en mars 2010 
    thetadirgoladeg = 270-2;
    % Antenne
    dessingolax = waist/10*[   0,   0:-.6:-4.8,     -4.8:.6:0, .34,.54,.6,.54,.34,  0];
    dessingolar = waist/10*[-1.8,1.8*ones(1,9),-1.8*ones(1,9),-1.2,-.6, 0, .6,1.2,1.8];
    
elseif strcmp(diag,'drefluc')
    % drefluc
    modex = 1;%2;
    Rgola = 428;
    waist = 39.7/2;
    dwaist = 12/2;
    phigoladeg = -2.7;
    zgola = 9.3/2;
    thetadirgoladeg = 1.5+180; 
elseif strcmp(diag,'difpol')
    % difpol
    frei = 60;
    Rgola = entree(Rgola,'Position en R des antennes (cm) ?',246.2);
    x_waist = Rgola-rr0;
    difpolgolahaut = 1;  % antenne haute (1) ou basse (0) de difpol
    difpolgolahaut = entree(difpolgolahaut,'Antenne : Haute (1) ou Basse (0) ?',1);
    %% caract?ristiques de l'optique gaussienne DIFPOL en 1994 ; LC 94
    if difpolgolahaut==0,
      zgola = -135.3;
      z_waist = -94.9;
      thetadirgoladeg = 90;
      waist = 35;  % dtheta = 2.58 a 60GHz / 5mm
    elseif difpolgolahaut==1,
      zgola = 135.3;
      z_waist = 73.3;
      thetadirgoladeg = -90;
      waist = 78;  % dtheta = 1.17;
    end;
    dwaist = abs(zgola-z_waist);
    phigoladeg = 0;
else
    % cas general
    Rgola = entree(Rgola,'position en R de l''antenne (cm) ?',0);
    zgola = entree(zgola,'position en z de l''antenne (cm) ?',0);
    dwaist = entree(dwaist,'distance avant le waist (cm) ?',12);
    waist = entree(waist,'waist du faisceau (mm) ?',30);
    phigoladeg = entree(phigoladeg,'position en phi de l''antenne (0=plan median entre bobines) ?',0);  
end

% Parametres de depart
xgola = Rgola-rr0;
ygola = 0;
phigola = phigoladeg * pi/180;
phidirgola = phidirgoladeg * pi/180;
dtheta = lambda/pi/waist*1000; % dtheta: demi-ouverture du faisceau
thetadirgola = thetadirgoladeg * pi/180;
thetarel = angle(-xgola-i*(zgola-zpos))-thetadirgola;
bimpact = abs(xgola+i*(zgola-zpos)) * tan(thetarel);
alphaa = asin(bimpact/a0r);

% Affichage des parametres de depart de l'onde
if modex, textmodex='X';textonde='Mode X';else textonde='Mode O';textmodex='O';end;
textonde = sprintf('%s\n F = %.1f GHz\n',textonde,frei);
if strcmp(diag,'difdop')
    textonde = sprintf('%s \\theta = %.1f�\n',textonde,difdoptheta);
else
    textonde = sprintf('%s \\theta = %.1f�\n',textonde,thetadirgoladeg);
end
textonde = sprintf('%s \\alphaa = %.1f�\n',textonde,alphaa*180/pi);
fprintf('%s\n',textonde)
  
if raydim==3,
    textDiagConfig = sprintf('  \\phi = %.3f�',phidirgoladeg);
end
textDiagConfig = sprintf('%s d\\theta = %.1f� ',textDiagConfig,dtheta*180/pi);
  
if figureprint==1
    textF(:)=sprintf('%g_%g',fix(frei),round(10*(frei-fix(frei))));
end
  
    
% Parametres non initialises
% Parametres du Runge-Kutta
rko = 2;  % ordre du Runge-Kutta
hmax = 2;  % pas maximum (distance en cm)
hmin = 1e-2; % pas minimum (cm pour le vide)
tol = 1e-4;  % tolerance pour 1 cm (dans le vide) optimum : 1e-4
testklim = .01;  % tolerance sur l'ecart entre k^2 et nx2 k_0^2

% Top depart...
chrono = now;
   
%%%%%%%%%%%%%%%%%%%%%%%
% VALEURS ANALYTIQUES %
%%%%%%%%%%%%%%%%%%%%%%%

if afficheanalytique
    % D'apres une communication de L. Colas (2002)
    fprintf(' Resultats analytiques comparables uniquement pour :\n ')
    fprintf('mode O, ')
    fprintf('pas de shift Safranov\n ')
    fprintf('profil de densit� parabolique, ')
    ne0=max(nefit);
    fprintf('ne(0)= %.2f 10^{19} m^{-3}\n ',ne0/1e19)
    fprintf('g�om�trie 2D (phi=0) ')
    fprintf('\n')
    nec = epsi0*me*w^2/qe^2;
    rhodifana = sqrt( (ne0/nec-1+sqrt(1+(ne0/nec)^2-2*ne0/nec*cos(2*alphaa))) / (2*ne0/nec) );
    kthetadifana = ko/100 * bimpact/(rhodifana*a0r);
    nx2difana = (sin(alphaa)/rhodifana)^2;
    fprintf(' Analytique : r/a=%.2f  N2=%.2f 2k_theta=%.2f cm^{-1}  \n\n',rhodifana,nx2difana,2*kthetadifana);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRACE DE FAISCEAU %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nombres de rayon (radialement et angulairement)
nbrayr = entree(nbrayr,'Nombre de rayons, radialement ? ',2);
%nbrayr = nbrayr - mod(nbrayr-1,-2);
if nbrayr==1, nbraya = 1; end
nbraya = entree(nbraya,'Nombre de rayons par tour ? (divisible par 4)',4);
if nbrayr>1, nbraya = ceil(nbraya/4)*4; end
rwmax = 1.6;  % rayon max du faisceau, relativement au waist
inrayrw = 1 + ceil((nbrayr-1)/rwmax);
if nbrayr>1
    [rayrad,rayang] = ndgrid((0:nbrayr)/(inrayrw-1),(0:nbraya-1)/nbraya*2*pi);
else
    [rayrad,rayang] = ndgrid([0,0],(0:nbraya-1)/nbraya*2*pi);
end
imageik = -rayrad.^2/(ko/100); % partie imaginaire de l'eikonale
if nbrayr>1
    inrayapol = [nbraya/4+1,nbraya*3/4+1];
    inrayator = [1,nbraya/2+1];
end

% Initialisation au waist
% Rotation du front d'onde d'angle phi et theta
rotphidirgola = [cos(phidirgola) sin(phidirgola) 0 ; -sin(phidirgola) cos(phidirgola) 0 ; 0 0 1];
rotthetadirgola = [cos(thetadirgola) 0 -sin(thetadirgola) ; 0 1 0 ; sin(thetadirgola) 0 cos(thetadirgola) ];
rotdirgola = rotphidirgola*rotthetadirgola;
xyzgola = [xgola;ygola;zgola];
uxyzwaist = rotdirgola*[1;0;0];  % direction de k
xyzwaist = xyzgola + dwaist*uxyzwaist;  % position du waist
% Rayon central rai(intemps) = [x,y,z,kx,ky,kz,E];
abscurv = 0;
rayrad = rayrad(1:nbrayr,:);
rayang = rayang(1:nbrayr,:);
waistbeam = zeros(nbrayr,nbraya,3);
waistbeam(:,:,2) = rayrad.*cos(rayang)*(waist/10);
waistbeam(:,:,3) = rayrad.*sin(rayang)*(waist/10);
% Faisceau : ray(intemps,inrad,inang) = [x,y,z,kx,ky,kz,E]
for inrayr = 1:nbrayr,
    for inraya = 1:nbraya;
        ray(1,inrayr,inraya,1:3) = xyzwaist + rotdirgola*squeeze(waistbeam(inrayr,inraya,:));
        ray(1,inrayr,inraya,4:6) = -uxyzwaist;
    end
end

% Chemin arriere vers l'antenne
[ray,abscurv] = equatray(ray,abscurv,0,floor(dwaist),inf);
abscurv = -abscurv(end:-1:1);
ray = ray(end:-1:1,:,:,:); % t -> -t
ray(:,:,:,4:6) = -ray(:,:,:,4:6); % k -> -k
inabscurvw = find(abscurv==0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trace vers le plasma
[ray,abscurv] = equatray(ray,abscurv,0,500,2.75*rhobord);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Valeur des parametres au rayon minimum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% recherche de l'indice ou diff(rho) change de signe
rhoy = rhox(ray);
nx2y = nx2x(ray);
for ir = 1:nbrayr
    for ia = 1:nbraya
        dnx2y = diff(nx2y(:,ir,ia));
        indif = max(find(nx2y(:,ir,ia)==min(nx2y(:,ir,ia))));
        if indif<length(nx2y(:,ir,ia))
            indifmid = ((indif-1/2)*dnx2y(indif)-(indif+1/2)*dnx2y(indif-1)) / (dnx2y(indif)-dnx2y(indif-1));
            dray = (ray(indif+1,ir,ia,:)-ray(indif-1,ir,ia,:))/2;
            ddray = (ray(indif+1,ir,ia,:)+ray(indif-1,ir,ia,:))/2-ray(indif,ir,ia,:);
            raydif(ir,ia,:) = ray(indif,ir,ia,:) + dray*(indifmid-indif) + ddray/2*(indifmid-indif)^2;
        else
            raydif(ir,ia,:) = ray(end,ir,ia,:);
        end
        inrdif(ir,ia) = indif;
    end
end
%raydif = squeeze(raydif);
xdif(:,:) = raydif(:,:,1);
ydif(:,:) = raydif(:,:,2);
zdif(:,:) = raydif(:,:,3);
kxdif(:,:) = ko/100 * raydif(:,:,4);
kydif(:,:) = ko/100 * raydif(:,:,5);
kzdif(:,:) = ko/100 * raydif(:,:,6);
[rhodif(:,:),kspher] = rhox(raydif,'kspher');
krdif(:,:) = ko/100 * kspher(:,:,1);
kphidif(:,:) = ko/100 * kspher(:,:,2);
kthetadif(:,:) = ko/100 * kspher(:,:,3);
kdif(:,:) = ko/100 * sqrt(sum(kspher.^2,3));
nx2dif(:,:) = nx2x(raydif);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Varargout
dif.ripple=ripple;
dif.nbrayr=nbrayr;
dif.nbraya=nbraya;
dif.waist=waist;
dif.rhodif=rhodif; 
dif.xdif=xdif;dif.ydif=ydif;dif.zdif=zdif;
dif.kdif=kdif; dif.kxdif=kxdif; dif.kydif=kydif; dif.kzdif=kzdif;
dif.krdif=krdif; dif.kthetadif=kthetadif; dif.kphidif=kphidif;
dif.nx2dif=nx2dif;
opt.rko=rko; opt.hmin=hmin; opt.hmax=hmax; opt.tol=tol; opt.testklim=testklim;
varargout{1}=dif;
varargout{2}=opt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Divergence et courbure du faisceau
if figurewaist & nbrayr>1 & nbboucle==1
    eclat = pi * (waist/10)^2 / (lambda*100);
    waistzth = (waist/10)*sqrt(1+abscurv{end}.^2/eclat^2);
    rcourbzth = (abs(abscurv{1})+eps) .* (1 + eclat^2./(abscurv{end}+eps).^2 );
    waistz{1} = sqrt(sum(squeeze(ray(:,inrayrw,inrayapol(1),1:3)-ray(:,1,1,1:3)).^2,2))./waistzth;
    waistz{2} = sqrt(sum(squeeze(ray(:,inrayrw,inrayapol(2),1:3)-ray(:,1,1,1:3)).^2,2))./waistzth;
    waistz{3} = sqrt(sum(squeeze(ray(:,inrayrw,inrayator(1),1:3)-ray(:,1,1,1:3)).^2,2))./waistzth;
    waistz{4} = sqrt(sum(squeeze(ray(:,inrayrw,inrayator(2),1:3)-ray(:,1,1,1:3)).^2,2))./waistzth;
    for inz = 1:size(ray,1)
        rpola = squeeze(ray(inz,inrayrw,inrayapol(1),1:3)-ray(inz,1,1,1:3));
        rpolb = squeeze(ray(inz,inrayrw,inrayapol(2),1:3)-ray(inz,1,1,1:3));
        matpol = [cross(rpola,rpolb)';rpola';rpolb'];
        if det(matpol)<1e-10, matpol = matpol+1e-5*eye(3); end
        vecpol = [0;norm(rpola)^2/2;norm(rpolb)^2/2];
        veccourb = matpol^(-1)*vecpol;
        rcourbz{1}(inz) = norm(veccourb);  %/rcourbzth(inz);
        if abs(norm(veccourb-rpola)-norm(veccourb))>1e-5*norm(veccourb), warning('erreur de calcul sur la courbure'), end
        if abs(norm(veccourb-rpolb)-norm(veccourb))>1e-5*norm(veccourb), warning('erreur de calcul sur la courbure'), end
        rtora = squeeze(ray(inz,inrayrw,inrayator(1),1:3)-ray(inz,1,1,1:3));
        rtorb = squeeze(ray(inz,inrayrw,inrayator(2),1:3)-ray(inz,1,1,1:3));
        mattor = [cross(rtora,rtorb)';rtora';rtorb'];
        if det(mattor)<1e-10, mattor = mattor+1e-5*eye(3); end
        vector = [0;norm(rtora)^2/2;norm(rtorb)^2/2];
        rcourbz{2}(inz) = norm(mattor^(-1)*vector);  %/rcourbzth(inz);
    end
    figure
    subplot(2,1,1)
    plot(abscurv{1},waistz{1,1},'r',abscurv{1},waistz{1,2},'r',abscurv{1},waistz{1,3},'b',abscurv{1},waistz{1,4},'b')
    ylim = get(gca,'ylim');set(gca,'ylim',[min(ylim(1),0),max(ylim(2),0)]);
    set(gca,'ylim',[.9,1.5]);
    set(gca,'xlim',[-20,300]);
    xlabel('abscisse curviligne (cm)');
    ylabel('divergence');
    title('Divergence relative [waist(z) / (w0*sqrt(1+(z/e)^2))]')
    subplot(2,1,2)
    plot(abscurv{1},rcourbzth,'k',abscurv{1},rcourbz{1,1},'r',abscurv{1},rcourbz{1,2},'b')
    ylim = get(gca,'ylim');set(gca,'ylim',[min(ylim(1),0),max(ylim(2),0)]);
    set(gca,'ylim',[0,2000]);
    set(gca,'xlim',[-20,300]);
    xlabel('abscisse curviligne (cm)');
    ylabel('rayon de courbure');
    title('Rayon de courbure relatif [Rcourbure(z) / z(1+e^2/z^2)]')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Affichage des resultats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
textfin =  sprintf('r/a = %.2f ',rhodif(1,1));
if nbrayr>1
    textfin = sprintf('%s(%+.2f,%+.2f) ',textfin,rhodif(inrayrw,inrayapol)-rhodif(1,1));
end
textfin = sprintf('%s\n 2k_{s} = %.2f cm-1',textfin,2*kdif(1,1));
if nbrayr>1
    dkdif = 2*kdif(inrayrw,inrayapol)-2*kdif(1,1);
    textfin = sprintf('%s (%+.2f,%+.2f)',textfin,dkdif);
end
textfin = sprintf('%s\n 2k_{\\theta} = %.2f cm-1',textfin,2*kthetadif(1,1));
if nbrayr>1
    dkthetadif = 2*kthetadif(inrayrw,inrayapol)-2*kthetadif(1,1);
    textfin = sprintf('%s (%+.2f,%+.2f)',textfin,dkthetadif);
end
if modex==1,
    textfin = sprintf('%s\n 2k_{r} = %.2f cm-1',textfin,2*krdif(1,1));
    if nbrayr>1
        dkrdif = 2*krdif(inrayrw,inrayapol)-2*krdif(1,1);
        textfin = sprintf('%s (%+.2f,%+.2f)',textfin,dkrdif);
    end
end

textfin2 = sprintf('2k_{\\phi} = %.2f cm-1',2*kphidif(1,1));
if nbrayr>1
    dkphidif = 2*kphidif(inrayrw,inrayator)-2*kphidif(1,1);
    textfin2 = sprintf('%s (%+.2f,%+.2f)',textfin2,dkphidif);
end
kxidif = angle(kthetadif(1,1)+i*kphidif(1,1)) * 180/pi;
textfin2 = sprintf('%s\n \\xi_k = %.1f?',textfin2,kxidif);
if nbrayr>1
    dkxidif =  angle(kthetadif(inrayrw,inrayator)+i*kphidif(inrayrw,inrayator)) * 180/pi - kxidif;
    textfin2 = sprintf('%s (%+.1f,%+.1f)',textfin2,dkxidif);
end
fprintf(' %s\n',textfin);
fprintf(' %s\n',textfin2);

textindice = sprintf('N^2 = %.2f ',nx2dif(1,1));
fprintf(' %s\n',textindice);

alphas = asin(bimpact/(rhodif(1,1)*a0r)); % angle au niveau de la surface de r?flexion
textonde = sprintf('%s \\alpha_s = %.1f�',textonde,alphas*180/pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   AFFICHAGE  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parois
if simulation
    paroix=[  3.1257    3.1252    3.1238    3.1214    3.1179    3.1136    3.1083 ...
    3.1020    3.0947    3.0866    3.0775    3.0675    3.0566    3.0448 ... 
    3.0321    3.0186    3.0042    2.9890    2.9729    3.1192    3.0997 ... 
    3.0794    3.0583    3.0365    3.0165    2.9931    2.9689    2.9441 ... 
    2.9187    2.8905    2.8640    2.8370    2.8094    2.7814    2.7529 ... 
    2.7240    2.6948    2.6652    2.6353    2.6052    2.5365    2.5084 ... 
    2.4802    2.4519    2.4234    2.3950    2.3666    2.3381    2.3098 ... 
    2.2816    2.2535    2.2256    2.1978    2.1704    2.1432    2.1163 ... 
    2.0897    2.0635    2.0377    2.0124    1.9875    1.9631    1.9393 ... 
    1.9160    1.8932    1.8711    1.8497    1.8289    1.8087    1.7893 ... 
    1.7707    1.7528    1.7357    1.7193    1.7038    1.6892    1.6754 ... 
    1.6625    1.6505    1.6393    1.6292    1.6199    1.6116    1.6042 ... 
    1.5978    1.5924    1.5879    1.5845    1.5820    1.5805    1.5800 ... 
    1.5805    1.5820    1.5845    1.5879    1.5924    1.5978    1.6042 ... 
    1.6116    1.6199    1.6292    1.6393    1.6505    1.6625    1.6754 ... 
    1.6892    1.7038    1.7193    1.7357    1.7528    1.7707    1.7893 ...
    1.8087    1.8289    1.8497    1.8711    1.8932    1.9160    1.9393 ... 
    1.9631    1.9875    2.0124    2.0377    2.0635    2.0897    2.1163 ...
    2.1472    2.1768    2.2067    2.2368    2.2671    2.2977    2.3100 ...
    2.3100    2.2747    2.2260    2.2200    2.2200    2.2205    2.2220 ...
    2.2244    2.2275    2.2311    2.2350    2.2530    2.2860    2.3190 ...
    2.3520    2.3850    2.4180    2.4510    2.4840    2.5170    2.5500 ...
    2.5830    2.6160    2.6490    2.6820    2.7150    2.7150    2.7150 ...
    2.7240    2.7529    2.7814    2.8094    2.8370    2.8640    2.8905 ...
    2.9164 2.9418 2.9665 2.9905 3.0138 3.0365 3.0583 3.0794 3.0997 3.1192 2.9729 2.9890 3.0042 3.0186 ...
    3.0321 3.0448 3.0566 3.0675 3.0775 3.0866 3.0947 3.1020 3.1083 3.1136 3.1179 3.1213 3.1238 3.1252 3.125];

    paroiz= [0    0.0279    0.0558    0.0836    0.1113    0.1389    0.1663 ...
    0.1935    0.2205    0.2472    0.2736    0.2997    0.3254    0.3507 ...
    0.3756    0.4000    0.4239    0.4474    0.4702    0.5455    0.5695 ...
    0.5928    0.6155    0.6373    0.6614    0.6818    0.7013    0.7200 ...
    0.7378    0.7514    0.7673    0.7823    0.7963    0.8094    0.8215 ...
    0.8326    0.8426    0.8517    0.8597    0.8666    0.8169    0.8214 ...
    0.8250    0.8275    0.8290    0.8295    0.8290    0.8275    0.8250 ...
    0.8214    0.8169    0.8114    0.8049    0.7974    0.7889    0.7795 ...
    0.7691    0.7578    0.7455    0.7324    0.7184    0.7035    0.6877 ...
    0.6711    0.6537    0.6354    0.6164    0.5967    0.5762    0.5550 ...
    0.5332    0.5107    0.4876    0.4639    0.4396    0.4147    0.3894 ...
    0.3636    0.3374    0.3107    0.2837    0.2563    0.2286    0.2007 ...
    0.1725    0.1440    0.1154    0.0867    0.0579    0.0289   -0.0000 ...
   -0.0289   -0.0579   -0.0867   -0.1154   -0.1440   -0.1725   -0.2007 ...
   -0.2286   -0.2563   -0.2837   -0.3107   -0.3374   -0.3636   -0.3894 ...
   -0.4148   -0.4396   -0.4639   -0.4876   -0.5107   -0.5332   -0.5550 ...
   -0.5762   -0.5967   -0.6164   -0.6354   -0.6537   -0.6711   -0.6877 ...
   -0.7035   -0.7184   -0.7324   -0.7455   -0.7578   -0.7691   -0.7795 ...
   -0.8426   -0.8517   -0.8597   -0.8666   -0.8725   -0.8774   -0.7946 ...
   -0.7446   -0.7446   -0.7372   -0.7372   -0.7316   -0.7277   -0.7241 ...
   -0.7210   -0.7186   -0.7171   -0.7166   -0.7166   -0.7166   -0.7166 ...
   -0.7166   -0.7166   -0.7166   -0.7166   -0.7166   -0.7166   -0.7166 ...
   -0.7166   -0.7166   -0.7166   -0.7166   -0.7166   -0.7446   -0.7946 ...
   -0.8326   -0.8215   -0.8094   -0.7963   -0.7823   -0.7673   -0.7514 ...
   -0.7345   -0.7168   -0.6982   -0.6787   -0.6584   -0.6373   -0.6155 ...
   -0.5928   -0.5695   -0.5455   -0.4702   -0.4474   -0.4239   -0.4000 ... 
   -0.3756   -0.3507   -0.3254   -0.2997   -0.2736   -0.2472   -0.2205 ...
   -0.1935   -0.1663   -0.1389   -0.1113   -0.0836   -0.0558   -0.0279    0.0000];

% position de paroi avant CIEL
%  paroiz = [0 34 68 102 136 169 202 235 267 299 358 391 423 454 485 514 542 569 595 620 643 665 686 704 722 733 748 ...
%    761 772 781 789 794 798 799 798 796 833 822 809 794 778 760 741 721 698 675 650 624 596 567 537 505 472 438 402 ...
%    365 327 288 249 208 168 126 84 42 0 -42 -84 -126 -168 -208 -249 -288 -327 -365 -402 -438 -472 -505 -537 -567 -596 ...
%    -624 -650 -675 -698 -721 -741 -760 -778 -794 -809 -822 -833 -788 -791 -791 -790 -786 -781 -774 -765 -754 -741 ...
%    -726 -722 -704 -686 -665 -643 -620 -595 -569 -542 -514 -485 -454 -423 -391 -358 -299 -267 -235 -202 -169 -136 ...
%    -102 -68 -34 0]/1000;
%  paroix = 2.35 + [760 760 757 753 748 741 733 723 712 699 736 718 699 679 656 633 608 582 554 526 496 465 432 399 365 ...
%   328 293 257 220 182 144 105 66 27 -13 -52 -100 -140 -179 -218 -256 -294 -330 -366 -401 -435 -467 -499 -530 -560 -588 ...
%    -615 -641 -664 -686 -706 -724 -740 -754 -766 -776 -784 -789 -793 -794 -793 -789 -784 -776 -766 -754 -740 -724 ...
%    -706 -686 -664 -641 -615 -588 -560 -530 -499 -467 -435 -401 -366 -330 -294 -256 -218 -179 -140 -100 -51 -12 27 66 ...
%    105 143 181 218 255 291 326 365 399 432 465 496 526 554 582 608 633 656 679 699 718 736 699 712 723 733 741 748 ...
%    753 757 760 760]/1000;
%paroix = 100*paroix-rr0;
%paroiz = paroiz*100;
end

paroixa = 100*paroix;
paroiza = 100*paroiz;% en cm

% Limites pour les figures
% PH08: affichage sur grille fixe centr�e sur la chambre
xlim = [150 450];
ylim = [-100 100];
zlimA = [-100 100];
couleurs = [ones(256,1) (0:255)'/300*ones(1,2)];  % du rouge clair a fonce
%couleurs = (1:256)'/384*[1 1 1];  % du noir au gris

% Limites pour evaluer Nx2
pasp = 1;
pasn = 2.5;
xlimn = [min(paroix),ceil(rhobord*a0r)+rr0+pasn];
ylimn = ylim;
zlimn = [-ceil(rhobord*a0r)-pasn,ceil(rhobord*a0r)+pasn];

% Figure poloidale x0z
if figurepoloidale 
    figure;set(gcf,'position',[100 10 550 700])
    if figuretoroidale
        subplot(211)
    end
    colormap(couleurs);
    hold on
    % Trace des decors
    plot(paroixa,paroiza,'b')
    % Trace des antennes
    dessingolapol = rotdirgola*[dessingolax;[0;1]*dessingolar];
    dessingolapolx = xgola + dessingolapol(1,:)+rr0;
    dessingolapolz = zgola + dessingolapol(3,:);%+zpos;
    plot(dessingolapolx,dessingolapolz,'k')
    
    [xp,zp] = meshgrid(xlimn(1)-rr0:pasp:xlimn(2)-rr0,zlimn(1):pasp:zlimn(2));
    xyz = zeros([size(xp,1),1,size(xp,2),3]);
    xyz(:,1,:,1) = xp;
    xyz(:,1,:,2) = 0;
    xyz(:,1,:,3) = zp;
    Nx2p = reshape(nx2x(xyz),size(xp));
    if isempty(nx2test)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Last Closed Magnetic Surface
        cercletheta = linspace(0,2*pi,400);
        x = a0r*cos(cercletheta+tau/a0r*(sin(cercletheta)).^3)+rr0;
        xe = a0r*cos(cercletheta)+rr0;% ellipticity alone
        z = (a0r*el)*sin(cercletheta)+zpos;
        %plot(x,z,'--','color',[.7 0 0],'linewidth',1);
        plot(xe,z,'r-.',x,z,'k--','linewidth',2);
        % Centre du plasma
        plot(0,zpos,'r.');
        plot(dsha0+rr0,zpos,'r+')
        % Contours des iso-indices
        if modex
            maxNx2=max(max(Nx2p));
            %if maxNx2>1
            %    [cs,h] = contour(xp,zp,Nx2p,[1 maxNx2]);
            %    clabel(cs,h,'color','r');
            %else
            [cs,h] = contour(xp+rr0,zp,Nx2p,[0:.2:1]);
            clabel(cs,h,'color','r','labelspacing',1000);
            %end
        else
            [cs,h] = contour(xp+rr0,zp,Nx2p,[0:.2:1]);
            clabel(cs,h,'color','r','labelspacing',1000);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    elseif nx2test~=1
        contour(xp+rr0,zp,Nx2p,1);
    end
    
    %trace de la couche fce & 2fce
    Rfce=27.99*B0*RR/frei;R2fce=2*(Rfce);
    plot([Rfce Rfce],100*[-1 1]);
    plot([R2fce R2fce],100*[-1 1],'b--');
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Faisceaux de rayons
    if figurearretrayon, tfin = inrdif(1,1); else tfin = size(ray,1); end
    plot(ray(1:tfin,1,1,1)+rr0,ray(1:tfin,1,1,3),'color',.2*ones(1,3));
    %  plot(ray(1:inrdif(1,1),1,1,1),ray(1:inrdif(1,1),1,1,3),'color',.2*ones(1,3));
    plot(xdif(1,1)+rr0,zdif(1,1),'xk')
    quiver(xdif(1,1)+rr0,zdif(1,1),3*kxdif(1,1),3*kzdif(1,1),0,'k')
    %if nbboucle==1 & nbrayr>1
    if  nbrayr>1
        inrayr = 2:nbrayr;
        for inra = inrayapol
            for inrr = inrayr
                if figurearretrayon, tfin = inrdif(inrr,inra); else tfin = size(ray,1); end
                plot(ray(1:tfin,inrayr,inra,1)+rr0,ray(1:tfin,inrayr,inra,3),'color',.4*ones(1,3));
                %	  plot(ray(1:inrdif(inrr,inra),inrayr,inra,1),ray(1:inrdif(inrr,inra),inrayr,inra,3),'color',.4*ones(1,3));
            end
            plot(xdif(inrayrw,inra)+rr0,zdif(inrayrw,inra),'xk')
            quiver(xdif(inrayrw,inra)+rr0,zdif(inrayrw,inra),3*kxdif(inrayrw,inra),3*kzdif(inrayrw,inra),0,'k')
        end
            plot(squeeze(ray(inabscurvw,end,inrayapol,1)+rr0),squeeze(ray(inabscurvw,end,inrayapol,3)),'k');
    end
    % end %%%%
    % Titre
    title(sprintf('%s',titrechoc),'HorizontalAlignment','right')
    %xlabel(['R(cm) [radial]   |   ',textDiagConfig],'HorizontalAlignment','left');
    ylabel('Z(cm) [vertical] | Plan poloidal');
    axis equal
    set(gca,'xlim',xlim,'ylim',zlimA);
    grid
    % Affichage des param?tres plasma et des r?sultats sur le faisceau
    text(150,100,textEqui,'HorizontalAlignment','left','VerticalAlignment','top');
    if simulation
        text(155,-110,textpla,'HorizontalAlignment','right','VerticalAlignment','top');
    end
    text(380,80,textonde,'HorizontalAlignment','left','VerticalAlignment','top')
    text(155,-80,textindice,'HorizontalAlignment','left','VerticalAlignment','bottom')
    text(480,-99,textfin,'HorizontalAlignment','right','VerticalAlignment','bottom')
    text(155,-100,sprintf('R(cm) [radial]'),'HorizontalAlignment','left','VerticalAlignment','bottom')
end
pause(.5);
if figuretoroidale==0
    if figureprint
        orient tall;
        print('-dpsc2',[num2str(chocnum),'Beam',textmodex,'F',textF(:)]);
    end;
end;
pause(.5);

% Figure toroidale x0y
if figuretoroidale  & raydim==3
    if figurepoloidale==1
        subplot(212)
    else
        figure
    end
    colormap(couleurs);
    hold on
    % Trace des decors
    arcphi = linspace(-1,1,100);
    arcphipp = arcphi(find(abs(arcphi)<asin(100/(rr0+max(paroixa)))));
    plot((max(paroixa))*cos(arcphipp),(max(paroixa))*sin(arcphipp),'b')
    arcphipm = arcphi(find(abs(arcphi)<asin(100/(rr0+min(paroixa)))));
    plot((min(paroixa))*cos(arcphipm),(min(paroixa))*sin(arcphipm),'b')
    % Trace des antennes
    dessingolator = rotdirgola*[dessingolax;[1;0]*dessingolar];
    dessingolatorx = xgola + dessingolator(1,:)+rr0;
    dessingolatory = ygola + dessingolator(2,:);
    plot(dessingolatorx,dessingolatory,'k')
    [xt,yt] = meshgrid(xlimn(1):pasp:xlimn(2),ylimn(1):pasp:ylimn(2));
    xyz = zeros([size(xt),1,3]);
    xyz(:,:,1,1) = xt;
    xyz(:,:,1,2) = yt;
    xyz(:,:,1,3) = 0;
    Nx2t = nx2x(xyz);
    if isempty(nx2test)
        % Trace de la derni?re surface magn?tique
        arcphimp = arcphi(find(abs(arcphi)<asin(100/(rr0+a0r))));
        plot((rr0+a0r)*cos(arcphimp),(rr0+a0r)*sin(arcphimp),'--','color',[.7 0 0],'linewidth',2);
        arcphimm = arcphi(find(abs(arcphi)<asin(100/(rr0-a0r))));
        plot((rr0-a0r)*cos(arcphimm),(rr0-a0r)*sin(arcphimm),'--','color',[.7 0 0],'linewidth',2);
        % Centre du plasma
        arcphicm = arcphi(find(abs(arcphi)<asin(100/rr0)));
        plot(rr0*cos(arcphicm),rr0*sin(arcphicm),'r.');
        arcphics = arcphi(find(abs(arcphi)<asin(100/(rr0+dsha0))));
        plot((rr0+dsha0)*cos(arcphics),(rr0+dsha0)*sin(arcphics),'r+');
        % Axe de symetrie du ripple
        if ripple
            phigolap = phigola+2*pi/18; phigolam = phigola-2*pi/18;
            plot([rr0-a0r,rr0+a0r]*cos(-phigola),[rr0-a0r,rr0+a0r]*sin(-phigola),'r-.','linewidth',.5)
            rlimp = [rr0-a0r,min(rr0+a0r,100/abs(sin(-phigolap)))];
            plot(rlimp*cos(-phigolap),rlimp*sin(-phigolap),'r-.','linewidth',.5)
            rlimm = [rr0-a0r,min(rr0+a0r,100/abs(sin(-phigolam)))];
            plot(rlimm*cos(-phigolam),rlimm*sin(-phigolam),'r-.','linewidth',.5)
        end
        % Contours des iso-indices
        if modex
            maxNx2t=max(max(Nx2t));
            %if maxNx2t>1
            %    [cs,h] = contour(xt,yt,Nx2t,[1 maxNx2t]);
            %    clabel(cs,h,'color','r');
            %else
            [cs,h] = contour(xt+rr0,yt,Nx2t,[0:.2:1]);
            clabel(cs,h,'color','r','labelspacing',1000);
            %end
        else
            [cs,h] = contour(xt+rr0,yt,Nx2t,[0:.2:1]);
            clabel(cs,h,'color','r','labelspacing',1000);
        end
    elseif nx2test~=1
        [cs,h] = contour(xt+rr0,yt,Nx2t,1);
    end
    % Faisceaux de rayons
    if figurearretrayon, tfin = inrdif(1,1); else tfin = size(ray,1); end
    plot(ray(1:tfin,1,1,1)+rr0,ray(1:tfin,1,1,2),'color',.2*ones(1,3));
    %    plot(ray(1:inrdif(1,1),1,1,1),ray(1:inrdif(1,1),1,1,2),'color',.2*ones(1,3));
    plot(xdif(1,1)+rr0,ydif(1,1),'xk')
    quiver(xdif(1,1)+rr0,ydif(1,1),3*kxdif(1,1),3*kydif(1,1),0,'k')
    if  nbrayr>1
        inrayr = 2:nbrayr;
        for inra = inrayator
            for inrr = inrayr
                if figurearretrayon, tfin = inrdif(inrr,inra); else tfin = size(ray,1); end
                plot(ray(1:tfin,inrayr,inra,1)+rr0,ray(1:tfin,inrayr,inra,2),'color',.4*ones(1,3));
                plot(ray(1:inrdif(inrr,inra),inrayr,inra,1)+rr0,ray(1:inrdif(inrr,inra),inrayr,inra,2),'color',.4*ones(1,3));
            end
            plot(xdif(inrayrw,inra)+rr0,ydif(inrayrw,inra),'xk')
            quiver(xdif(inrayrw,inra)+rr0,ydif(inrayrw,inra),3*kxdif(inrayrw,inra),3*kydif(inrayrw,inra),0,'k')
        end
        plot(squeeze(ray(inabscurvw,end,inrayator,1)),squeeze(ray(inabscurvw,end,inrayator,2)),'k');
    end
    % end
    % Titre
    %title(sprintf('Toroidal Plane'),'HorizontalAlignment','left','VerticalAlignment','bottom')
    xlabel(['R(cm) [radial]  ',username,'   ',dateplot])
    ylabel('y(cm) [toroidal] | Plan toroidal');
    axis equal
    set(gca,'xlim',xlim,'ylim',ylim)
    grid
    % Affichage des param?tres plasma et des r?sultats sur le faisceau
    text(460,-99,textfin2,'HorizontalAlignment','right','VerticalAlignment','bottom')
    text(460,-50,textDiagConfig,'HorizontalAlignment','right');
end
pause(.5);
    


if figureprint
    orient tall;
    print('-dpsc2',[num2str(chocnum),'Beam',textmodex,'-','F',textF(:)]);
end;
pause(.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Figure 3D
if figure3d & nbboucle==1 & raydim==3
  figure
  colormap(couleurs);
  hold on
  [x3,y3,z3] = meshgrid(xlimn(1):pasn:xlimn(2),ylimn(1):pasn:ylimn(2),zlimn(1):pasn:zlimn(2));
  xyz = zeros([size(x3),3]);
  xyz(:,:,:,1) = x3; 
  xyz(:,:,:,2) = y3; 
  xyz(:,:,:,3) = z3; 
  Rho3 = rhox(xyz);
  Nx23 = nx2x(xyz);
  % Dessin de l'antenne
  tour = 0:pi/8:pi;
  dessingolax = ones(length(tour),1)*dessingolax;
  dessingolay = cos(tour')*dessingolar;
  dessingolaz = sin(tour')*dessingolar;
  dessingola3 = rotdirgola * [dessingolax(:)';dessingolay(:)';dessingolaz(:)'];
  dessingola3x = xgola + reshape(dessingola3(1,:),size(dessingolax));
  dessingola3y = ygola + reshape(dessingola3(2,:),size(dessingolax));
  dessingola3z = zgola + reshape(dessingola3(3,:),size(dessingolax));
  hgola = surf(dessingola3x,dessingola3y,dessingola3z);
  set(hgola,'facecolor',[.3,.3,.3],'edgecolor','black');
  if isempty(nx2test)
    % Centre du plasma
    plot(-rr0+rr0*cos(arcphicm),rr0*sin(arcphicm),'r.');
    plot(-rr0+(rr0+dsha0)*cos(arcphics),(rr0+dsha0)*sin(arcphics),'r+');
    % Trace de la couche de coupure
    hcoup = patch(isosurface(x3,y3,z3,Nx23,0));
    set(hcoup,'facecolor',[1,0,0],'edgecolor',[1,0,0]);
    if str2num(strtok(version,'.'))>=6,
      set(hcoup,'facealphaa',.2,'edgealpha',.05);
    end
  end
  if isempty(nx2test) | nx2test~=1
    % Trace de la limite du plasma
    if str2num(strtok(version,'.'))>=6,
      ha0r = patch(isosurface(x3,y3,z3,Rho3,rhobord-eps));
      set(ha0r,'facecolor',[1,0,0],'edgecolor',[1,0,0]);
      set(ha0r,'facealpha',.05,'edgealpha',.02);
    end
  end
  % Faisceaux de rayons
  plot3(ray(:,1,1,1),ray(:,1,1,2),ray(:,1,1,3),'color',.2*ones(1,3));
  plot3(xdif(1,1),ydif(1,1),zdif(1,1),'xk')
  quiver3(xdif(1,1),ydif(1,1),zdif(1,1),3*kxdif(1,1),3*kydif(1,1),3*kzdif(1,1),0,'k')
  %if nbboucle==1 & nbrayr>1
  if  nbrayr>1
      for inrr = 2:nbrayr
          inra = [1:nbraya,1];
          color = [(inrr-1)/(nbrayr-1) 0 (nbrayr-inrr)/(nbrayr-1)];
          hfais = surf(squeeze(ray(:,inrr,inra,1)),squeeze(ray(:,inrr,inra,2)),squeeze(ray(:,inrr,inra,3)));
          set(hfais,'facecolor',.8*color,'edgecolor',.1*color);
          if str2num(strtok(version,'.'))>=6,
              set(hfais,'facealpha',1,'edgealpha',1);
          end
    
          
          for inra = 1:nbraya
              quiver3(ray(end,inrr,inra,1),ray(end,inrr,inra,2),ray(end,inrr,inra,3), ...
                  3*ray(end,inrr,inra,4),3*ray(end,inrr,inra,5),3*ray(end,inrr,inra,6),0,'k')
          end
      end
      global iplr ipla
      if ~isempty(iplr)
          quiver3(ray(end,iplr,ipla,1),ray(end,iplr,ipla,2),ray(end,iplr,ipla,3), ...
              3*ray(end,iplr,ipla,4),3*ray(end,iplr,ipla,5),3*ray(end,iplr,ipla,6),0,'g')
      end
  end
axis('square')
axis equal
axis([xlim,ylim,zlimA])
view(-15,5)
% Titre
title(sprintf('%s figure 3D  %s',titrechoc))
xlabel('x(cm) [radial]');
ylabel('y(cm) [toroidal]');
zlabel('z(cm) [vertical]');
end



if afficherayons & nbboucle==1 
  if nbray==3, fprintf('                 -dtheta dtheta -dphi   dphi\n'); end
  fprintf('rho     : %s\n',sprintf('%6.2f ',[rhodif(1,1),rhodif(inrayrw,[inrayadop,inrayator])]));
  fprintf('nx2     : %s\n',sprintf('%6.2f ',nx2dif));
  fprintf('kr      : %s cm-1\n',sprintf('%6.2f ',krdif));
  fprintf('2kphi   : %s cm-1\n',sprintf('%6.2f ',2*kphidif));
  fprintf('2ktheta : %s cm-1\n',sprintf('%6.2f ',2*kthetadif));
end

% Top fin
temps2str(now-chrono);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rho,out,Yriprel] = rhox(xyz,modecal)

% Calcul de rho a partir de la position cartesienne
% Modifi? par Seak, Teng-Fong  1998-11-24, CH 2002
% rho = rhox(xyz)
%   calcul de rho, avec correction ripple et shafranov
% [rho,xphi,Yriprel] = rhox([x,y,z],'ripple')
%   xphi, Yriprel : corrections ripple sur x et B
% [rho,kspher] = rhox([x,y,z,kx,ky,kz],'kspher') 
%   kpher=[kr,kphi,ktheta] par rotation de kx,ky,kz
% [rho,kcart] = rhox([x,y,z,kr,kphi,ktheta],'kcart')
%   kcart=[kx,ky,kz] par rotation inverse

global rr0 a0r zpos el tau dsha0 dshap
global phigola rhobord ripple

dim = size(xyz);
rho = ones([dim(1:end-1),1]);
xyz = reshape(xyz,[prod(dim(1:end-1)),dim(end)]);

x = xyz(:,1); y = xyz(:,2); z = xyz(:,3); 

% Rotation toroidale : x,y -> xphi,0
%phi = angle((x+rr0)+i*y);
phi = atan2(y,(x+rr0));
xphi = sqrt((x+rr0).^2+y.^2)-rr0;  

% Correction ripple pour xphi et z
rr0rip = 220;
rrip = sqrt((xphi+rr0-rr0rip).^2+z.^2);
%thetarip = angle(xphi+rr0-rr0rip+i*z);
thetarip = atan2(z,xphi+rr0-rr0rip);
% Formule inverse approch?e issue de [TS 20.91.02 P. Hertout & A. Seigneur]
drrip = 5.8e-3 * exp(5.5e-2*rrip-1.4e-1*thetarip.^2); % 1er ordre
%drrip = 5.8e-3 * exp(5.5e-2*(rrip-drrip.*(cos(18*(phi+phigola))-1))-1.4e-1*thetarip.^2); % 2eme ordre
drrip = ( sqrt(xphi.^2+z.^2) <= a0r*rhobord ) .* drrip;  % correction seulement si rho<=rhobord
drrip = ripple .* drrip;  % pas de correction si ripple=0
rrip0 = rrip - drrip .* (cos(18*(phi+phigola))-1);
xphi = rrip0.*cos(thetarip)+rr0rip-rr0;
z = rrip0.*sin(thetarip);

% Rho
rho(:) = sqrt((xphi./a0r).^2+((z-zpos)./(el*a0r)).^2); % hors champ magnetique
% Correction Shafranov pour rho<1
inrho = find(rho<1);
if (dsha0~=0) & (~isempty(inrho)),
  xin = xphi(inrho);
  zin = (z(inrho)-zpos)./el;
  % Formule nanalytique si dshap=2 :       ______________________________
  %             [  2                      / 4      2                 2  2   ]1/2 / [  _    ]
  % rho(x,z) =  [ a  - 2 d  (x - d ) - \ / a  - 4 a d  (x - d ) - 4 d  z    ]   /  [ V2 d  ]
  %             [  0      0       0     \   0      0 0       0       0      ]  /   [     0 ]
  rho(inrho) = sqrt( a0r^2 - 2*dsha0*(xin-dsha0) - sqrt(a0r^4-4*a0r^2*dsha0*(xin-dsha0)-4*dsha0^2*zin.^2) )  ...
               /  sqrt(2)/dsha0;
  % Correction si l'exposant du Shafranov est different de 2
  if dshap~=2,
    % suite convergeant toujours vers rho recherche @CH 2002
    rho1 = rho(inrho); rhoi = inf;
    while max(abs(rho1-rhoi))>2e-3 
      rhoi = rho1;  % on conserve le rho precedent
      dsha = dsha0*(1-rhoi.^dshap);  % on evalue le shafranov avec le rho precedent
      rho1 = sqrt((xin-dsha).^2 + zin.^2) / a0r;  % on recalcule rho avec ce shafranov
    end
    rho(inrho) = rho1;
  end;
end;

if nargout>=2 & modecal=='ripple'
  % Ripple, Correction sur B (formule empirique pour TS)
  if length(dim)>2
    xphi = reshape(xphi,dim(1:end-1)); 
    phi = reshape(phi,dim(1:end-1)); 
    z = reshape(z,dim(1:end-1));
  end
  uu = sqrt( 1 + 5.2e-3*(xphi+rr0-204) - sqrt(1+1.04e-2*(xphi+rr0-204)-2.7e-5*z.^2) );
  Yriprel = 2.2e-4 * exp(13.6*uu+11.83*uu.^2).*cos(18.*(phi+phigola));% ajout PH cos(18*(phi+phigola))-1) ???
  %Yriprel = 2.2e-4 * exp(13.6*uu+11.83*uu.^2);
  out = xphi;
elseif nargout>=2
  % Calcul de la matrice de transformation : x y z -> rho yphi ztheta
  if modecal=='kspher'
    exprot = 1;
  elseif modecal=='kcart'
    exprot = -1;
  end  
  % Angle de correction toroidale due au ripple : tan(dphirip) = 1/R dR/dphi
  dphirip = atan( cos(thetarip) / (rrip0.*cos(thetarip)+rr0rip) * 18*drrip.*sin(18*(phi+phigola)) );
  % Angle theta poloidal
  dsha = dsha0 * (1-rho(:).^dshap) .* (rho(:)<=1);
  theta = angle(xphi-dsha+i*z);
  for in = 1:length(x)
    % Matrice de rotation toroidale : x+rr0 y z -> xphi+rr0 yphi z
    Rphi = [cos(phi(in)) sin(phi(in)) 0 ; -sin(phi(in)) cos(phi(in)) 0 ; 0 0 1];
    % Matrice de rotation toroidale pour la correction Ripple : xphi+rr0 yphi z -> xphirip+rr0 yphi zrip
    Rphirip = [cos(dphirip(in)) sin(dphirip(in)) 0 ; -sin(dphirip(in)) cos(dphirip(in)) 0 ; 0 0 1];
    % Matrice de transformation poloidale pour la correction Ripple
    %Rthetarip = ......... ...... .. 
    % Matrice de rotation poloidale : xphi-dsha yphi z -> rho yphi ztheta
    Rtheta = [cos(theta(in)) 0 sin(theta(in)) ; 0 1 0 ; -sin(theta(in)) 0 cos(theta(in)) ];
    % Matrice de transformation totale : x y z -> rho yphi ztheta
    Rot = Rtheta * Rphirip * Rphi;
    % exprot=1 : kx ky kz -> kr kphi ktheta
    % exprot=-1 : kr kphi ktheta -> kx ky kz
    % kr kphi ktheta : perpendiculaire et tangentiel au champ magnetique
    %                  dirige vers les r, phi et theta croissants
    kin = xyz(:,4:6);
    kout(in,:) = (Rot^exprot*kin(in,:)')'; 
  end
  out = reshape(kout,[dim(1:end-1),3]);
end


function ne_spline = nerhofit(rhofit,nefit)

% Fitting the density profile with splines
% Densite (m-3) nefit as a function of rhofit
% 

% methodes d'interpolation par des polynomes d'ordre 3 entre donnees
% spline : garantit la continuite de la derivee seconde entre segments
% pchip : garantit la monotonie entre les points (cf >>doc pchip) 
interpol = 'spline'; 


% On calcule l'interpolation entre points 

% Determination des coefficients du spline cubique (ou par pchip)
ne_spline = feval(interpol,rhofit,nefit);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ne = nerho(rho)

% Densite (m-3) en fonction de rho : nefit d�fini sur rhofit
% Vecteurs de taille quelconque
% 

global ne_spline
% methodes d'interpolation par des polynomes d'ordre 3 entre donnees
% spline : garantit la continuite de la derivee seconde entre segments
% pchip : garantit la monotonie entre les points (cf >>doc pchip) 
interpol = 'spline'; 


% Version simlplifiee de ppval de matlab ...
% if necessary, sort rho 
lx = length(rho);rho=reshape(rho,1,lx);
l = length(ne_spline.breaks)-1;
[rho,inrho] = sort(rho);
% for each data point, compute its breakpoint interval
[ignored,index] = sort([ne_spline.breaks(1:end-1) rho]);
index = max([find(index>l)-(1:lx);ones(1,lx)]);
% now go to local coordinates ...
rho = rho-ne_spline.breaks(index);
% ... and apply nested multiplication:
ne = ne_spline.coefs(index,1);
for i=2:4
  ne = rho'.*ne + ne_spline.coefs(index,i);
end
ne(inrho) = ne;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nx2 = nx2x(xyz)

% Calcul de l'indice a partir de la position

global modex ripple rhobord rr0 a0r nx2test
global B0 RR j0 jp Te0 Tepiq
global w
% Constantes physiques
c = 3e8;  % [m/s]
epsi0 = 8.8541e-12;
me = 9.1e-31;
qe = 1.6e-19;
mu0 = 4*pi*1e-7;

%if modex
  [rho,xphi,Yriprel] = rhox(xyz,'ripple');
%else
%  rho = rhox(xyz);
%end

% exterieur au plasma
% nx2(find(rho>rhobord)) = 1;
dim = size(xyz); 
dim(end) = 1;
nx2 = ones(dim);

% interieur au plasma
inp = find(rho<=rhobord);
if isempty(inp), return, end

neoc = epsi0*me*w^2/qe^2;
v = nerho(rho(inp)) / neoc;
nx2(inp) = 1-v;
ttt=(qe/me /w)^2;
save v_ts v ttt

if modex,
  Bvide=B0.*RR./(rr0+xphi(inp));% Attention, B0 en RR=2.36 pas en RR0
                                   %Bvide=B0./(1+xphi(inp)/rr0);% old
                                   
  %% Ripple
  Brip = Bvide .* Yriprel(inp); % Yriprel from function rho
  if ~ripple, Brip = 0*Brip; end
                                  
  B2 = (Bvide-Brip).^2;
if modex == 2 % B avec corrections dia et paramagn�tique, Bpol
    
    %% Champ poloidal : mu0*j0*a^2/(2*(jp+1)) * (1-(1-rho^2)^(jp+1))/rho
    %Bpol= mu0*j0*a0r.^2/(2*(jp+1)) * (1-(1-rho(inp).^2)^(jp+1))./ (rho(inp)+eps) .* (rho(inp)<1); 
    Bpol= mu0*j0*a0r./(2*(jp+1)) * (1-(1-rho(inp).^2).^(jp+1))./ (rho(inp)+eps) .* (rho(inp)<1); 
    Bpol_bord= mu0*j0*a0r./(2*(jp+1)) ./ (rho(inp)) .* (rho(inp)>1) ;
    Bpol=Bpol+Bpol_bord;
       % ajouter formule rho>1
    
    %% Diamagnetisme: correction � Bz^2: -2*mu0*P
    Te_K =zeros(size(rho(inp)));
    Te_K = Te0 * 1.6e-19 * (1-rho(inp).^2).^(Tepiq) .* (rho(inp)<1);
    CorDia= -2*mu0* nerho(rho(inp)) .* Te_K *2;
    
    %% Paramagnetisme: correction a Bz^2 :
    %%                 int_r^a (2*mu0*j*Bpol) = mu0^2*j0^2 *a^2/(jp+1)int_r^a  1/rho *(1-rho^2)^jp * (1-(1-rho^2)^(jp+1))
    CorPara=zeros(size(rho(inp)'));
    n_pas_para=100;drho_para=1/(n_pas_para-1);x_para=(0:drho_para:1)';
    rho_para=x_para*ones(size(rho(inp)'));  % matrice n_pas_para x l_para (size rho <1)
    k_rho=find((rho_para-ones(n_pas_para,1)* abs(rho(inp)') )<0);
    CorPara_temp= ( (1-rho_para.^2).^jp  - (1-rho_para.^2).^(2.*jp+1) )./(rho_para +eps); % signe formule AS/RS
    CorPara_temp(k_rho)=zeros(size(k_rho));
    CorPara_int = trapz(x_para,CorPara_temp);
    CorPara= mu0^2.*j0.^2./(jp+1).*a0r.^2.*CorPara_int';
        
    %% correction de B: devrait etre: B^2= (Bvide^2+corPara+corDia)+(Bpol^2)
    B2 = (Bvide-Brip).^2 + CorPara+CorDia + Bpol.^2;
  end
  u = (qe/me /w)^2 .* B2; % wce^2 / w^2
  nx2(inp) = ((1-v).^2-u) ./ (1-u-v);
end

% Indice nx2 test
if ~isempty(nx2test),
  nx2(inp) = nx2test; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function raydt = dhraysnell(rays,h)

% DHRAY  @ YM, LC 22-06-94, CH 2002
% Calcul de l'int?gration du trace de rayon sur un pas h 
% Cas limite quand l'indice varie beaucoup sur une faible distance
% On traite le cas comme un saut d'indice ? travers une surface
% On applique alors les lois de Descartes de reflexion / transmission

global debog

% ray : position initiale
ray(1:3) = rays(1:3);
% indice initial
nx2 = nx2x(ray);
% rayk : k initial
rayk(1:3) = rays(4:6);

% dray : deplacement attendu
dray = h*rayk;
% indice attendu apres le pas
nx2out = nx2x(ray+dray);
% Oups... l'indice varie en fait peu -> easy
if abs(nx2out-nx2)<1e-5
  % On avance d'un pas, k ne change pas
  raydt = [dray,0,0,0];
  return
end

% ray+part*dray : position juste avant l'interface (a h/1000 pres)
part = 0; 
for inp = 1:10
  partin = part + 1/2^(inp);
  nx2in = nx2x(ray+partin*dray);
  if abs(nx2in-nx2out) > abs(nx2in-nx2)
    part = partin;
  end
end
% On cherche la normale ? l'interface...
% Pour cela on va chercher 2 tangeantes
% drn,dro : d'abord, il nous faut 2 vecteurs orthogonaux ? dray
drphi = angle(dray(1)+i*dray(2));
drtheta = angle(norm(dray(1:2))+i*dray(3));
rotphi = [cos(drphi) -sin(drphi) 0 ; sin(drphi) cos(drphi) 0 ; 0 0 1];
rottheta = [cos(drtheta) 0 -sin(drtheta) ; 0 1 0 ; sin(drtheta) 0 cos(drtheta) ];
drn = norm(dray) * rotphi*rottheta * [0;1;0];
dro = norm(dray) * rotphi*rottheta * [0;0;1];
% drayn : dans le plan (dray,drn), on cherche la tangente ? l'interface
dirn = 0;
for ina = 1:12
  dirin = dirn + pi/2^ina;
  drayn = -cos(dirin)*dray + sin(dirin)*drn';
  nx2in = nx2x(ray+part*dray+drayn);
  if abs(nx2in-nx2out) > abs(nx2in-nx2)
    dirn = dirin;
  end
end
% drayo : dans le plan (dray,dro), on cherche la tangente ? l'interface
diro = 0;
for ina = 1:12
  dirin = diro + pi/2^ina;
  drayo = -cos(dirin)*dray + sin(dirin)*dro';
  nx2in = nx2x(ray+part*dray+drayo);
  if abs(nx2in-nx2out) > abs(nx2in-nx2)
    diro = dirin;
  end
end
% La normale est orthogonale aux 2 tangentes
normale = cross(drayn,drayo);
normale = sign(dot(rayk,normale)) * normale / norm(normale);

% Composante de k suivant la normale
raykn = dot(rayk,normale);
% Loi de Descartes - Snell
rayknout = sqrt(nx2out/nx2*norm(rayk)^2 - (norm(rayk)^2-raykn^2));  % transmission
if ~isreal(rayknout)
  rayknout = -raykn;  % reflexion 
end
% Variation de k
raydt(4:6) = (rayknout-raykn) * normale';
% Deplacement total
raydt(1:3) = h*rayk + (1-part)*h*raydt(4:6);
if debog>=2, fprintf('(%.1f,%.1f,%.1f)',normale), end
if debog>=1, fprintf('?'), end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function raydt = dhray(rays,h)

% DHRAY  @ YM, LC 22-06-94, CH 2002
% faisceau gaussien [Nowak & Orefice, Phys. Plasmas _1_(5), pp1242-1250 (1994)]

global nbrayr nbraya imageik 
global debog

raydt = zeros(size(rays));

% Variation de dx = k Cdt
% Pas = h * Cdt = h
raydt(:,:,1:3) = h*rays(:,:,4:6);

% Variation de dk = 1/2 (grad(nx2) + grad(I^2)) Cdt
% Effet d'indice
drn = .3*h;
raysra = rays(:,:,1:3);
dray1 = zeros(nbrayr,nbraya,3);
dray1(:,:,1) = drn;
dray2 = zeros(nbrayr,nbraya,3);
dray2(:,:,2) = drn;
dray3 = zeros(nbrayr,nbraya,3);
dray3(:,:,3) = drn;
raysd(:,:,1,:) = raysra+dray1/2;
raysd(:,:,2,:) = raysra+dray2/2;
raysd(:,:,3,:) = raysra+dray3/2;
raysd(:,:,4,:) = raysra-dray1/2;
raysd(:,:,5,:) = raysra-dray2/2;
raysd(:,:,6,:) = raysra-dray3/2;
nx2d = nx2x(raysd);
gradnx2 = (nx2d(:,:,1:3)-nx2d(:,:,4:6))/drn;
raydt(:,:,4:6) = h*gradnx2/2;

if nbrayr>1
  % Matrice de la partie imaginaire de l'eikonale
  dri = norm(squeeze(raydt(1,1,1:3)))/10/norm(squeeze(diff(raysra(1:2,1,1:3))));
  imageikra = imageik(2:nbrayr,:);
  imageikrpa = imageik(3:nbrayr+1,:);
  imageikrma = imageik(1:nbrayr-1,:);
  % Matrices des positions adjacentes pour calculer les gradients
  raysra = raysra(2:nbrayr,:,:);
  rayksra = rays(2:nbrayr,:,4:6);
  raysrma(2:nbrayr-1,:,1:3) = raysra(1:nbrayr-2,:,:);
  raysrma(1,1:nbraya,1:3) = rays(1,1:nbraya,1:3);
  raysrpa(1:nbrayr-2,:,1:3) = raysra(2:nbrayr-1,:,:);
  raysrpa(nbrayr-1,1:nbraya,1:3) = 2*raysra(nbrayr-1,1:nbraya,:)-raysrma(nbrayr-1,1:nbraya,:);
  raysrap = raysra(:,[2:nbraya,1],:);
  raysram = raysra(:,[nbraya,1:nbraya-1],:);
  % Gradient de la partie imaginaire de l'eikonale : gradI
  dprojsrpa = sum((raysrap-raysram).*(raysrpa-raysra),3) ./ sqrt(sum((raysrap-raysram).^2,3));
  dperpsrpa = sqrt(sum((raysrpa-raysra).^2,3)-dprojsrpa.^2);
  dprojsrma = sum((raysrap-raysram).*(raysra-raysrma),3) ./ sqrt(sum((raysrap-raysram).^2,3));
  dperpsrma = sqrt(sum((raysra-raysrma).^2,3)-dprojsrma.^2);
  gradipsra = (imageikrpa-imageikra) ./ dperpsrpa;
  gradimsra = (imageikra-imageikrma) ./ dperpsrma;
  gradisra = (dperpsrpa.*gradimsra+dperpsrma.*gradipsra) ./ (dperpsrpa+dperpsrma);
  gradisrma(1:nbrayr-1,:) = [ zeros(1,nbraya) ; gradisra(1:nbrayr-2,:) ];
  gradisrpa(1:nbrayr-2,:) = gradisra(2:nbrayr-1,:);
  gradisrpa(nbrayr-1,:) = 2*gradisra(nbrayr-1,:)-gradisrma(nbrayr-1,:);
  gradisrap = gradisra(:,[2:nbraya,1]);
  gradisram = gradisra(:,[nbraya,1:nbraya-1]);
  % Idem pour s+ds
  rayspra = raysra+raydt(2:nbrayr,:,1:3)/dri;
  raysprma(2:nbrayr-1,:,1:3) = rayspra(1:nbrayr-2,:,:);
  raysprma(1,1:nbraya,1:3) = rays(1,1:nbraya,1:3)+raydt(1,1:nbraya,1:3)/dri;
  raysprpa(1:nbrayr-2,:,1:3) = rayspra(2:nbrayr-1,:,:);
  raysprpa(nbrayr-1,1:nbraya,1:3) = 2*rayspra(nbrayr-1,1:nbraya,:)-raysprma(nbrayr-1,1:nbraya,:);
  raysprap = rayspra(:,[2:nbraya,1],:);
  rayspram = rayspra(:,[nbraya,1:nbraya-1],:);
  dprojsprpa = sum((raysprap-rayspram).*(raysprpa-rayspra),3) ./ sqrt(sum((raysprap-rayspram).^2,3));
  dperpsprpa = sqrt(sum((raysprpa-rayspra).^2,3)-dprojsprpa.^2);
  dprojsprma = sum((raysprap-rayspram).*(rayspra-raysprma),3) ./ sqrt(sum((raysprap-rayspram).^2,3));
  dperpsprma = sqrt(sum((rayspra-raysprma).^2,3)-dprojsprma.^2);
  gradipspra = (imageikrpa-imageikra) ./ dperpsprpa;
  gradimspra = (imageikra-imageikrma) ./ dperpsprma;
  gradispra = (dperpsprpa.*gradimspra+dperpsprma.*gradipspra) ./ (dperpsprpa+dperpsprma);
  gradisprma(1:nbrayr-1,:) = [ zeros(1,nbraya) ; gradispra(1:nbrayr-2,:) ];
  gradisprpa(1:nbrayr-2,:) = gradispra(2:nbrayr-1,:);
  gradisprpa(nbrayr-1,:) = 2*gradispra(nbrayr-1,:)-gradisprma(nbrayr-1,:);
  gradisprap = gradispra(:,[2:nbraya,1]);
  gradispram = gradispra(:,[nbraya,1:nbraya-1]);
  % Idem pour s-ds
  raysmra = raysra-raydt(2:nbrayr,:,1:3)/dri;
  raysmrma(2:nbrayr-1,1:nbraya,1:3) = raysmra(1:nbrayr-2,1:nbraya,:);
  raysmrma(1,1:nbraya,1:3) = rays(1,1:nbraya,1:3)-raydt(1,1:nbraya,1:3)/dri;
  raysmrpa(1:nbrayr-2,1:nbraya,1:3) = raysmra(2:nbrayr-1,1:nbraya,:);
  raysmrpa(nbrayr-1,1:nbraya,1:3) = 2*raysmra(nbrayr-1,1:nbraya,:)-raysmrma(nbrayr-1,1:nbraya,:);
  raysmrap = raysmra(:,[2:nbraya,1],:);
  raysmram = raysmra(:,[nbraya,1:nbraya-1],:);
  dprojsmrpa = sum((raysmrap-raysmram).*(raysmrpa-raysmra),3) ./ sqrt(sum((raysmrap-raysmram).^2,3));
  dperpsmrpa = sqrt(sum((raysmrpa-raysmra).^2,3)-dprojsmrpa.^2);
  dprojsmrma = sum((raysmrap-raysmram).*(raysmra-raysmrma),3) ./ sqrt(sum((raysmrap-raysmram).^2,3));
  dperpsmrma = sqrt(sum((raysmra-raysmrma).^2,3)-dprojsmrma.^2);
  gradipsmra = (imageikrpa-imageikra) ./ dperpsmrpa;
  gradimsmra = (imageikra-imageikrma) ./ dperpsmrma;
  gradismra = (dperpsmrpa.*gradimsmra+dperpsmrma.*gradipsmra) ./ (dperpsmrpa+dperpsmrma);
  gradismrma(1:nbrayr-1,:) = [ zeros(1,nbraya) ; gradismra(1:nbrayr-2,:) ];
  gradismrpa(1:nbrayr-2,:) = gradismra(2:nbrayr-1,:);
  gradismrpa(nbrayr-1,:) = 2*gradismra(nbrayr-1,:)-gradismrma(nbrayr-1,:);
  gradismrap = gradismra(:,[2:nbraya,1]);
  gradismram = gradismra(:,[nbraya,1:nbraya-1]);
  % Gradient(gradI^2) (connu pour 3 directions)
  dray(:,:,1,:) = raysrap-raysram;
  dray(:,:,2,:) = raysrpa-raysrma;
  dray(:,:,3,:) = rayspra-raysmra;
  dgradi2(:,:,1) = gradisrap.^2-gradisram.^2;
  dgradi2(:,:,2) = gradisrpa.^2-gradisrma.^2;
  dgradi2(:,:,3) = gradispra.^2-gradismra.^2;
  detdray = dot(dray(:,:,:,1),cross(dray(:,:,:,2),dray(:,:,:,3),3),3);
  % Recherche des cas ou les 3 vecteurs sont coplanaires
  plan = abs(detdray)./sqrt(sum(dray(:,:,1,:).^2,4))./sqrt(sum(dray(:,:,2,:).^2,4))./sqrt(sum(dray(:,:,3,:).^2,4));
  global iplr ipla
  [iplr,ipla] = find(plan<1e-2);
  for in = 1:length(iplr)
    % vecteurs coplanaires : on fixe un gradient nul suivant la normale
    inplr = iplr(in); inpla = ipla(in);
    vec1 = raysrap(inplr,inpla,:)-raysram(inplr,inpla,:);
    vec2 = raysrpa(inplr,inpla,:)-raysrma(inplr,inpla,:);
    dray(inplr,inpla,3,:) = cross(vec1,vec2,3);
    dgradi2(inplr,inpla,3) = 0;
    detdray(inplr,inpla) = dot(dray(inplr,inpla,:,1),cross(dray(inplr,inpla,:,2),dray(inplr,inpla,:,3),3),3);
    planin = abs(detdray(inplr,inpla))./sqrt(sum(dray(inplr,inpla,1,:).^2,4)) ...
	     ./sqrt(sum(dray(inplr,inpla,2,:).^2,4))./sqrt(sum(dray(inplr,inpla,3,:).^2,4));
    if min(min(planin))<1e-2,
      if debog>=2, fprintf('[%d:%d:%.0f]',inplr,inpla,10*log10(planin)); end
    end
  end
  % la matrice dgradi2 est inversee 
  gradgradi2 = zeros(nbrayr,nbraya,3);
  gradgradi2(2:nbrayr,:,1) = dot(dgradi2(:,:,:),cross(dray(:,:,:,2),dray(:,:,:,3),3),3)./detdray;
  gradgradi2(2:nbrayr,:,2) = dot(dgradi2(:,:,:),cross(dray(:,:,:,3),dray(:,:,:,1),3),3)./detdray;
  gradgradi2(2:nbrayr,:,3) = dot(dgradi2(:,:,:),cross(dray(:,:,:,1),dray(:,:,:,2),3),3)./detdray;
  % grad  = dray^-1 * dgrad^2
  raydt(:,:,4:6) = raydt(:,:,4:6) + h*gradgradi2/2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dray,matprecision] = dhrayrk(rays,h)

% Integration par Runge-Kutta (cf: numerical recepies)

global rko tol
% Parametres Cash-Karp pour le Runge-Kutta 5eme ordre
rka2 = 1/5; rka3 = 3/10; rka4 = 3/5; rka5 = 1; rka6 = 7/8;
rkb21 = 1/5; 
rkb31 = 3/40; rkb32 = 9/40; 
rkb41 = 3/10; rkb42 = -9/10; rkb43 = 6/5;
rkb51 = -11/54; rkb52 = 5/2; rkb53 = -70/27; rkb54 = 35/27;
rkb61 = 1631/55296; rkb62 = 175/512; rkb63 = 575/13824; rkb64 = 44275/110592; rkb65 = 253/4096;
rkc1 = 37/378; rkc2 = 0; rkc3 = 250/621; rkc4 = 125/594; rkc5 = 0; rkc6 = 512/1771;
rkd1 = 2825/27648; rkd2 = 0; rkd3 = 18575/48384; rkd4 = 13525/55296; rkd5 = 277/14336; rkd6 = 1/4;

if rko==5
  dray1 = dhray(rays,h);
  dray2 = dhray(rays+rkb21*dray1,h);
  dray3 = dhray(rays+rkb31*dray1+rkb32*dray2,h);
  dray4 = dhray(rays+rkb41*dray1+rkb42*dray2+rkb43*dray3,h);
  dray5 = dhray(rays+rkb51*dray1+rkb52*dray2+rkb53*dray3+rkb54*dray4,h);
  dray6 = dhray(rays+rkb61*dray1+rkb62*dray2+rkb63*dray3+rkb64*dray4+rkb65*dray5,h);
  dray = rkc1*dray1+rkc3*dray3+rkc4*dray4+rkc5*dray5+rkc6*dray6;
  drayom = rkd1*dray1+rkd3*dray3+rkd4*dray4+rkd5*dray5+rkd6*dray6;
elseif rko==2
  dray1 = dhray(rays,h);
  dray2 = dhray(rays+dray1,h);
  dray3 = dhray(rays+dray1/4+dray2/4,h);
  dray = (dray1+4*dray3+dray2)/6;
  drayom = (dray1+dray2)/2;
end
raysp = rays + dray;
testks = sum(rays(:,:,4:6).^2,3)./nx2x(rays);
testksp = sum(raysp(:,:,4:6).^2,3)./(nx2x(raysp)+eps);
%matprecision = max(max(abs(dray-drayom)/tol/h,[],3),abs(testks-testksp)/10/tol/h);
matprecision = max(abs(dray-drayom)/tol/h,[],3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ray,abscurv] = equatray(ray,abscurv,dt,tfin,rhofin)

% EQUATRAY   @ YM, LC 10/11/94, STF 19990204
% CH 2002 erreur calculee composante par composante
% ray = [x(n) y(n) z(n) kx(n) ky(n) kz(n)] matrice du rayon
% tfin : nombre max de periodes d'onde sur lequel est integre le trace 
% dt : nombre max de periodes d'onde entre 2 points du trace

global nbrayr nbraya rko hmax hmin testklim
global debog

% Initialisation de l'integration
rays = reshape(ray(end,:,:,:),[nbrayr,nbraya,6]); % faisceau au dernier pas calcule
rhos = rhox(rays(1,1,1:3)); % rho du rayon central au dernier pas calcule
if size(ray,1)>1
  rhosm = rhox(ray(end-1,1,1,1:3));  % rho du rayon central a l'avant dernier pas
else
  rhosm = nan;
end
abscurvs = abscurv(end);  % temps initial
tvec = abscurvs+(dt:dt:tfin);  % temps intermediaires pour lesquels on veut connaitre ray
if dt==0, tvec = abscurvs+tfin; end  % si dt==0  pas de temps intermediaires
inlimits = 1; % Variable test sur les limites du trace 
h = hmax;  % pas de calcul initial

% Boucle sur les temps intermediaires a atteindre
for tnext = tvec
  % Boucle sur les pas de calcul
  while abscurvs<tnext & h>=hmin & inlimits 
    % Le pas ne doit pas depasser le point intermediaire tnext
    if abscurvs+h>tnext,
      h = tnext-abscurvs; 
    end
    [dray,matprecision] = dhrayrk(rays,h);
    precision = max(max(matprecision));  % Precision relative
    if (precision<=1) | (h<=hmin)
      dnx2s = nx2x(rays+dray)-nx2x(rays);
      [inr,ina] = find(matprecision>1);
      for in = 1:length(inr)
        dray(inr(in),ina(in),:) = dhraysnell(rays(inr(in),ina(in),:),h);
	precision = max(max(matprecision(find(matprecision<1))));
	if isempty(precision), precision = 1; end
      end
      % Si la precision est <=1, on avance d'un pas
      abscurvs = abscurvs + h;
      abscurv = [abscurv;abscurvs];
      rays = rays + dray;
      ray(end+1,:,:,:) = rays;
      rhosm = rhos;
      rhos = min(min(rhox(rays(:,:,1:3))));
      inlimits = 1;~((rhos>=rhofin) & (rhosm<rhofin)) & (rhos<2.75);
      if debog>=1, fprintf('#'); end
      if debog>=2
	fprintf('%.0f',h/hmin);
	testk = sum(rays(1,1,4:6).^2,3)/nx2x(rays(1,1,:));
	if abs(1-testk)>testklim, fprintf(' k^2/N^2k0^2:%.2f',testk); end
	fprintf(' %.0f:%.0f:%.0f %.2f:%.2f:%.2f %.2f\n',rays(1,1,:),nx2x(rays(1,1,:))); 
      end
      % Adaptation du pas pour garder une precision proche de 1 ou un pas proche de hmax
      h = min(hmax,.9*h*precision^(-1/(rko-1)));  % rko-1 si tol proportionnel a h
      h = max(h,hmin);
    else
      if debog>=1, fprintf('@'), end
      h = min(hmax,.8*h*max(precision^(-1/(rko-1)),.1));
      h = max(h,hmin);
    end
  end
end
% Si le trace n'atteint ni tfin ni les limites, c'est qu'il y a une singularite : h->hmin 
if abscurvs<tvec(end) & inlimits
  fprintf('Singularit? sur le trace ! ')
  if debog>=1
    fprintf('(0,0) : %.2f \n',matprecision(1,1));
    for inrr = 1:nbrayr-1
      fprintf('(%d,:) : ',inrr);
      fprintf('%.1f ',matprecision(inrr+1,:));
      fprintf('\n');
    end
  end
end
% Test sur la condition : k^2 = N^2 k0^2 pour le rai central
testk = norm(squeeze(rays(1,1,4:6)))^2/nx2x(rays(1,1,:));
if abs(1-testk)>testklim
  fprintf('rai central k^2/N^2k0^2 = %.2f ~= 1 \n',testk)
  if debog>=1
    testky = squeeze(sum(rays(:,:,4:6).^2,3))./nx2x(rays);
    for inrr=1:nbrayr-1
      fprintf('(%d,:) : ',inrr);
      fprintf('%.2f ',testky(inrr+1,:));
      fprintf('\n');
    end
  end
end
if debog>=1, fprintf('\n'), end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rayout = propagvide(rayin,rholim);

xyz = rayin(1:3); 
kxyz = rayin(4:6); 
% Recherche a 5 cm pres sur 3 m
in = (0:5:300)';
xyz = ones(size(in)) * xyz + in * kxyz;
rho = rhox(xyz);
inbord = min(find(diff(rho<rholim)~=0));
if isempty(inbord)
  rayout = [];
else
  % Recherche a 1 mm pres
  in = (0:.1:5)';
  xyz = ones(size(in)) * xyz(inbord,:) + in * kxyz;
  rho = rhox(xyz);
  inbord = min(find(diff(rho<rholim)~=0));
  % Recherche a 1/10 mm pres
  in = (0:.01:.12)';
  xyz = ones(size(in)) * xyz(inbord,:) + in * kxyz;
  rho = rhox(xyz);
  inbord = find(diff(rho<rholim)~=0)+1;
  rayout = [xyz(inbord,:),kxyz];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Out = entree(In,Message,Defaut,TypeIn)

%ENTREE demande ? l'utilisateur ou donne une valeur par d?faut ?
%une variable
%   Out = entree(In,Message,Defaut,TypeIn)
% In: Donn?e ? traiter
% Message: Message pour la demande (si In==nan ou '?')
% Defaut: Valeur par d?faut (si In == [])
% TypeIn: type de donn?e 's'=string

if nargin == 3
  TypeIn = 'n';
end
if strcmp(TypeIn,'n')&strcmp(TypeIn,'s')
  error('Type de donn?e ''n'' ou ''s''')
end
if (length(Defaut)>1)&strcmp(TypeIn,'n')
%  Defaut = Defaut(1);
end  
if isempty(Defaut)
  if strcmp(TypeIn,'s')
    Defaut ='';
  else
    Defaut = nan;
  end
end

quest = 0;
if ~isempty(In)
  if strcmp(TypeIn,'s')&strcmp(In,'?')
    In = input(sprintf('%s [%s] ',Message,Defaut),'s');
    quest = 1;
  end
  if strcmp(TypeIn,'n')&isnan(In)
    In = input(sprintf('%s [%g] ',Message,Defaut));
    quest = 1;
  end
end

if isempty(In)
  Out = Defaut;
else
  Out = In;
end

if 0& ~quest & ~isempty(inputname(1))
  if strcmp(TypeIn,'s')
    fprintf('?  %s = %s\n',inputname(1),Out);
  else
    fprintf('?  %s = %g\n',inputname(1),Out);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TempsStr = temps2str(Temps,Prec);
%TEMPS affiche la duree d'un calcul
% string = temps2str(duree,prec);

if nargin < 2
  Prec = 1;
end  
%  duree: temps (en sec)
if length(Temps)==1
  TempsVec = datevec(Temps);
else
  TempsVec = Temps;
end
Annee = TempsVec(1);
Mois = TempsVec(2);
Jour = TempsVec(3);
Heure = TempsVec(4);
Mn = TempsVec(5);
Sec = floor(TempsVec(6)/Prec)*Prec;
TempsStr = '';
if Annee > 1000
  TempsStr = sprintf('%s%02d',TempsStr,Jour);
  TempsStr = sprintf('%s/%02d',TempsStr,Mois);
  TempsStr = sprintf('%s/%04d ',TempsStr,Annee);
elseif (Jour ~= 0) & (Jour ~= 31)
  TempsStr = sprintf('%s%dj ',TempsStr,Jour);
end  
if Mn ~= 0
  if Heure ~= 0
    TempsStr = sprintf('%s%02dh',TempsStr,Heure);
  end  
  TempsStr = sprintf('%s%02d''',TempsStr,Mn);
end
if Prec >= 1
  TempsStr = sprintf('%s%02d''''',TempsStr,Sec);
else  
  TempsStr = sprintf('%s%02f''''',TempsStr,Sec);
end  
if nargout == 0
  fprintf(1,'%s\n',TempsStr);
  clear TempsStr
end
