function [neprofile]=difdop_ne4beam(varargin)
% set the density profile used for beam tracing
% in the calculation of the optical index
% uses experimental measurements :
%       (1) analytic form using central ne and peaking factor from
%       interferometry
%       (2) interpolation using gne from interferometry
%       (3) TPROF 
%       (4) reflectometry with EFIT for converting Rc to rhofit
%       (5) reflectometry with TMAG for converting Rc to rhofit
%
%       (0) for viewing all profiles type before choosing
%
% varargin: 1: numchoc
%           2: times (t0FA & t0Fp)
%           3: type_profile
%           4: filename (and path)

couleur=['b' 'c' 'r' 'm' 'k' 'g' 'y' 'b' 'c' 'r' 'm'];

% shot number
if nargin>0
    numchoc = varargin{1};  % numero de choc
else
    numchoc = input(' choc ');
end
disp(['Choc ' int2str(numchoc)]);

% time vector
if nargin>1
    time = varargin{2};
    if isstruct(time);
        t0FA=time.t0FA;
        t0Fp=time.t0Fp;
    else
        t0FA=time;t0Fp=time;
    end
    nbtA=length(t0FA);nbtimes=length(t0Fp);
else
    filemat=['Sdifdop',int2str(numchoc),'.mat'];
    load(filemat);clear Skw* dtDifdop nfft taille_seg* Polarisation difdop_anglepol F_pal*
    nbtA=length(t0FA);nbtimes=length(t0Fp);
end

% how is evaluated the density profile
% 1/analytic representation from experimental central density (interferometry)
%            and peaking factor (from <ne> and nec): edge profile not well known
% 2/interferometry
% 3/from TPROF but only for rho<1, still pb in the SOL
% 4-5/ Laure: from reflectometry
if nargin>2
    type_profile = varargin{3};
else
    disp('density profile ? : (0) all profiles available for comparison; ')
    disp('                    (1) analytic');
    disp('                    (2) interp with interferometry');
    disp('                    (3) TPROF ');
    disp('                    (4) experimental with EFIT, (5) with TMAG  ');
    type_profile=input('density profile ? : [default : 0] ?');
    if isempty(type_profile),type_profile=0;end
end

if type_profile==0
    type_profile0=input('avec reflectometry (0/1) ? : [default : 0] ?');
    if isempty(type_profile0),type_profile0=0;end
else
    type_profile0=0;
end
    
if type_profile>3 || type_profile0==1
    if nargin>3
        fileneref = varargin{4};
    else
        [fileneref,PathName,FilterIndex] = uigetfile('*.mat');
    end
end

if type_profile==1
    prof_type='A';
elseif type_profile==2
    prof_type='I';
elseif type_profile==3
	prof_type='T';
elseif type_profile==4
	prof_type='EXE';
%    load([fileneref])
    load([PathName,fileneref])
    disp(['Fichier utilise ', fileneref]);
    % le fichier doit contenir une info normalis�e
    % nefit rhofit 
    % nefit matrice: 2eme dimension est le temps ->moyenne � int�grer dans la boucle temps
    % est-ce que la sortie est d�j� "smooth"?
elseif type_profile==5
	prof_type='EXT';
    load([PathName,fileneref])
    %load([fileneref])
    disp(['Fichier neprofile utilise ', fileneref]);
    if ~exist('nefit'),error('fichier profil de densit� reflec ne contient pas nefit'),end;
end

% rho for the interpolation
rhop0 = [0:.05:.85 .88:.02:1.24 1.3 1.4 1.5 1.7 2 2.5 3];rhobord0=1.25;
ksep=find(rhop0==1);% faire un test sur le bord

%%%%% reading the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gne used by all profile determination?
if type_profile==0 | prof_type=='A' | prof_type=='I'| prof_type=='T'
    [vgne,tn,rhone] = tsbase(numchoc,'gne');
    vne0 =  vgne(:,1);vnea =  vgne(:,end);
    [vnmoy,tn] = tsbase(numchoc,'snmoy');
    nelambda=2;
    nelambda0=2/70;
    %nelambda = nelambdan/a0r;  % normalisation au petit rayon
    %nepiq = ne0/ne-1; %1/(ne0/ne-1); calculated after, for each time
end

if type_profile==0 | prof_type=='T'
    % recherche de toutes les occurences TPROF disponibles
    ocurlist=[];
    filePROF=['Prof',int2str(numchoc),'-'];
    for ocur=0:3 % test de l'existence d'un traitement TPROF
        if exist([filePROF,int2str(ocur),'.mat'])==2,
            profilne=1;ocurlist=[ocurlist ocur];
            numchoc1=numchoc+.0;numchoc1=numchoc1+ocur/10;
            [GPROFNEFIT,GPROFt,GPROFrho]=tsbase(numchoc1,'GPROFNEFIT');
            [GPROFNE,GPROFtd]=tsbase(numchoc1,'GPROFNE');
            [GPROFRHONE,GPROFtd]=tsbase(numchoc1,'GPROFRHONE');            
            disp(sprintf('\n l ocurrence %d existe pour fichier Prof%d',ocur,numchoc));
            disp(sprintf('temps :%2.2f � %2.2f ',min(GPROFt),max(GPROFt)));
            disp(sprintf('nombre de points :%2.2f  ',length(GPROFt)));
            disp(sprintf('diff(t)  :%2.2f  ',mean(diff(GPROFt))));
            
            PROFNEFIT{ocur+1}=GPROFNEFIT;
            PROFNE{ocur+1}=GPROFNE;
            PROFRHONE{ocur+1}=GPROFRHONE;
            tPROF{ocur+1}=GPROFt;
            tPROFd{ocur+1}=GPROFtd;
        end;
    end
    nplot=length(ocurlist);

    % choose the "occurrence" of TPROF or analytic profile (display all)
    % after displaying each trigger time (each angle)
    for k_iter=1:nbtA
        figure(k_iter+1000);
        % for analytic profile
        ne0= mean(vne0(tn>=t0FA(k_iter) &tn<=t0FA(k_iter)+.1));
        nea= mean(vnea(tn>=t0FA(k_iter) &tn<=t0FA(k_iter)+.1));
        nmoy= mean(vnmoy(tn>=t0FA(k_iter) &tn<=t0FA(k_iter)+.1));
        nepiq = ne0/nmoy-1; %1/(ne0/ne-1);
        ne1r=0.15;nea=ne1r*ne0;%nea/ne0; % tmp
        necoeur = ne0 * ( ne1r + (1-ne1r) * (1-rhop0.^2).^nepiq ) .*(rhop0>=0).*(rhop0<=1);
        nesol = ne0*ne1r * exp(-(sqrt(rhop0)-1)/nelambda0) .* (rhop0>1).*(rhop0<=rhobord0);
        nep0 = necoeur + nesol;
        plot(rhop0,nep0,'r--');hold on;

        % interferometry
        ne= mean(vgne(tn>=t0FA(k_iter) &tn<=t0FA(k_iter)+.1,:));
        neinterp = interp1(rhone,ne,rhop0).* (rhop0>=0).*(rhop0<=1.0);
        nebord=neinterp(ksep);neinterp(isnan(neinterp))=0;
        nesolI = nebord * exp(-(sqrt(rhop0)-1)/nelambda0) .* (rhop0>1).*(rhop0<=rhobord0);
        nepI = neinterp + nesolI;
        
        plot(rhone,ne,'m+',rhop0,nepI,'m')

        % reflectometry
        if type_profile>3 || type_profile0==1
            if(sizene(1)~=nbtA & sizene(1)==nbtA),
                plot(rhofit',nefit')
            else
                plot(rhofit,nefit)
            end
        end
            
        % 2/ toutes les occurences TPROF
        for ocur=ocurlist
            % TPROF times at the acquisition time?
            
            if t0FA(k_iter)>min(tPROF{ocur+1}) & (t0FA(k_iter)+.1)<max(tPROF{ocur+1})
                neprofdata=mean(PROFNE{ocur+1}(tPROFd{ocur+1}>=t0FA(k_iter)&tPROFd{ocur+1}<=t0FA(k_iter)+.1,:));
                nerhodata=mean(PROFRHONE{ocur+1}(tPROFd{ocur+1}>=t0FA(k_iter)&tPROFd{ocur+1}<=t0FA(k_iter)+.1,:));
                plot(nerhodata,neprofdata,'b*');
            
                nefitdata = mean(PROFNEFIT{ocur+1}(tPROF{ocur+1}>=t0FA(k_iter)&tPROF{ocur+1}<=t0FA(k_iter)+.1,:));
                plot(GPROFrho,nefitdata,'k--');                
                necoeure=interp1(GPROFrho,nefitdata,rhop0).* (rhop0>=0).*(rhop0<=1.0);
                nebord=necoeure(ksep);
                necoeure(isnan(necoeure))=0;
                nesole = nebord * exp(-(sqrt(rhop0)-1)/nelambda0) .* (rhop0>1).*(rhop0<=rhobord0);
                nepe{ocur+1}(k_iter,:) = necoeure + nesole;
                plot(rhop0,nepe{ocur+1}(k_iter,:),[couleur(ocur+1),'--']);xlabel(['rhofit, TRPOF occurrence ',num2str(ocur)]);
            else
                disp('Caution !!! acquisition time different from TPROF time')
                disp(sprintf('acquisition times : ts = %d te %d',t0FA(k_iter),t0FA(k_iter)))
                disp(sprintf('TPROF times : t1 = %d t2 = %d',min(tPROF{ocur+1}),max(tPROF{ocur+1})))
                ag=input('Use it anyway ? [0 or 1]');
                if ag==1
                    if min(tPROF{ocur+1})>=t0FA(k_iter) & (t0FA(k_iter)+.1)<=max(tPROF{ocur+1})
                        neprofdata=mean(GPROFNE(tPROF{ocur+1}>=min(tPROF{ocur+1})+1,:));
                        nerhodata=mean(GPROFRHONE(tPROF{ocur+1}>=min(tPROF{ocur+1})+1,:));
                        plot(nerhodata,neprofdata,'b*');
                        
                        nefitdata = mean(GPROFNEFIT(tPROF{ocur+1}>=min(tPROF{ocur+1})+1,:));
                        necoeure=interp1(nerhodata,nefitdata,rhop0).* (rhop0>=0).*(rhop0<=1.0);
                        nebord=necoeure(rhop0==1);necoeure(isnan(necoeure))=0;
                        nesole = nebord * exp(-(sqrt(rhop0)-1)/nelambda0) .* (rhop0>1).*(rhop0<=rhobord0);
                        nepe{ocur+1}(k_iter,:) = necoeure + nesole;
                        plot(rhop0,nepe{ocur+1}(k_iter,:),couleur(ocur+1));xlabel(['rhofit, TPROF occurrence ',num2str(ocur)]);                    
                    elseif t0FA(k_iter)>=max(tPROF{ocur+1}) &  min(tPROF{ocur+1})<=t0FA(k_iter)
                        neprofdata=mean(GPROFNE(tPROF{ocur+1}>=max(tPROF{ocur+1})-1,:));
                        nerhodata=mean(GPROFRHONE(tPROF{ocur+1}>=max(tPROF{ocur+1})-1,:));
                        plot(nerhodata,neprofdata,'b*');
                        nefitdata = mean(GPROFNEFIT(tPROF{ocur+1}>=max(tPROF{ocur+1})-1,:));
                        necoeure=interp1(nerhodata,nefitdata,rhop0).* (rhop0>=0).*(rhop0<=1.0);
                        nebord=necoeure(rhop0==1);necoeure(isnan(necoeure))=0;
                        nesole = nebord * exp(-(sqrt(rhop0)-1)/nelambda0) .* (rhop0>1).*(rhop0<=rhobord0);
                        nepe{ocur+1}(k_iter,:) = necoeure + nesole;
                        plot(rhop0,nepe{ocur+1}(k_iter,:),couleur(ocur+1));xlabel(['rhofit, TRPOF occurrence ',num2str(ocur)]);
                    end
                end
            end
        end;
    title([sprintf('TS# %d, d %d, t=%2.3g s',numchoc,k_iter,t0FA(k_iter))],'Fontsize',20)%,'Fontweight','bold')
    end

    % choice here
    disp('density profile ? : (1) analytic')
    disp('                    (2) interp with interferometry');
    disp('                    (3) TPROF ');
    disp('                    (4) experimental with EFIT, (5) with TMAG  ');
    type_profile_choix=input('density profile ? : [default : 2] ?');
    if isempty(type_profile_choix);type_profile_choix=2;end
else
    type_profile_choix=type_profile;
end   
    
if type_profile_choix==1
    neprofile.prof_type='A';
    neprofile.times=t0Fp;nbtimes=length(t0Fp);
    for k_trigger=1:nbtimes
        ne0= mean(vne0(tn>=neprofile.times(k_trigger) &tn<=neprofile.times(k_trigger)+0.01));
        nea= mean(vnea(tn>=neprofile.times(k_trigger) &tn<=neprofile.times(k_trigger)+0.01));
        nmoy= mean(vnmoy(tn>=neprofile.times(k_trigger) &tn<=neprofile.times(k_trigger)+0.01));
        nepiq = ne0/nmoy-1; %1/(ne0/ne-1);
        ne1r=0.15;nea=ne1r*ne0;%nea/ne0; % tmp otherwise ne1r=nea/ne0;
        necoeur = ne0 * ( ne1r + (1-ne1r) * (1-rhop0.^2).^nepiq ) .*(rhop0>=0).*(rhop0<=1);
        nesol = ne0*ne1r * exp(-(sqrt(rhop0)-1)/nelambda0) .* (rhop0>1).*(rhop0<=rhobord0);
        nep0 = necoeur + nesol;
        neprofile.ne(k_trigger,:)=nep0;
    end
    neprofile.rho=rhop0;
    neprofile.nelambda=nelambda;
elseif type_profile_choix==2
    neprofile.prof_type='I';
    neprofile.times=t0Fp;nbtimes=length(t0Fp);  
    for k_trigger=1:nbtimes
        ne= mean(vgne(tn>=neprofile.times(k_trigger) &tn<=neprofile.times(k_trigger)+0.01,:));
        necoeur = interp1(rhone,ne,rhop0).* (rhop0>=0).*(rhop0<=1.0);
        nebord=necoeur(ksep);necoeur(isnan(necoeur))=0;
        nesol = nebord * exp(-(sqrt(rhop0)-1)/nelambda0) .* (rhop0>1).*(rhop0<=rhobord0);
        nep0 = necoeur + nesol;
        neprofile.ne(k_trigger,:)=nep0;
    end
    neprofile.rho=rhop0;
    neprofile.nelambda=nelambda;
elseif type_profile_choix==3
    disp('ok')
    choix_ocur=input('quelle occurence TPROF utiliser pour le profil de densit? ?');
    if isempty(choix_ocur);choix_ocur=0;end
	neprofile.prof_type='T';
	neprofile.tprof=choix_ocur;
    neprofile.times=t0FA;   
    neprofile.ne=nepe{ocur+1};
    neprofile.rho=rhop0;
    neprofile.nelambda=nelambda;
end


if prof_type=='EXE'
    neprofile.prof_type='EXE';
    neprofile.times=t0FA;
    neprofile.ne=nefit;neprofile.rho=rhofit;
    neprofile.nelambda=[];
    figure(100)
    title('Profile from reflecto using EFIT equilibrium')
    plot(neprofile.rho',neprofile.ne')
    xlabel('rho')
    ylabel('ne')
end

if prof_type=='EXT'            
    neprofile.prof_type='EXT';neprofile.ne=nefit;neprofile.rho=rhofit;
    neprofile.nelambda=[];
    neprofile.times=t0FA;nbtimes=length(t0FA); 
    %  alignement sur  rhop0 = [0:.05:.85 .88:.02:1.24 1.3 1.4 1.5 1.7 2 2.5 3];rhobord0=1.25;
    rhopref = [0:.02:.6 .61:.01:1.24 1.3 1.4 1.5 1.7 2 2.5 3];
    for k_trigger=1:nbtimes-3
%        knan=find(isnan(nefit(1,:)));
%        rhofit1=rhofit(k_trigger,(max(knan)+1):end);
%        nefit1=nefit(k_trigger,(max(knan)+1):end);
        knan=find(~isnan(nefit(1,:)));
        rhofit1=rhofit(k_trigger,knan);
        nefit1=nefit(k_trigger,knan);
        if rhofit1(1)~=0
            k0=find(rhopref<min(rhofit1));
            rhofit1=[rhopref(k0) rhofit1];
            nefit1=[nefit1(1)*ones(size(rhopref(k0))) nefit1];
            
        end
        if rhofit1(end)<3
            kmax=find(rhopref>max(rhofit1));
            rhofit1=[rhofit1 rhopref(kmax)];
            nefit1=[nefit1 zeros(size(rhopref(kmax)))];
        end
        neprofile.ne(k_trigger,:)=interp1(rhofit1,nefit1,rhopref);neprofile.rho=rhopref;
        clear rhofit1 nefit1
    end
    figure(1000)
    title('Profile from reflecto using TMAG equilibrium')
    plot(neprofile.rho',neprofile.ne')
    xlabel('rho')
    ylabel('ne')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% density profile %%%%%%%%%%%%%%%%%%%
