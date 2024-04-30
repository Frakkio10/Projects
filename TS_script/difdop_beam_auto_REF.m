function difdop_beam_auto(varargin)
% Trace de rayon/faisceau pour les experiences TS 
% varargin: 1: numchoc
%           2: polarisation pour difdop ou config dreve
%           3: type_profile
%           4: filename (and path) [if ne from DREFLEC]


%--------------------------------------
% lecture  parametres
% a adapter pour une lecture des parametres Top 
% et signaux lents (inclinometre)
%
% ce qu'il faut passer comme argument au coeur de programme
% DIFDOP Frequence
% DIFDOP Angle pol
% DIFDOP angle tor
% DIFDOP Polarisation
% PLASMA density profile
% PLASMA equilibrium
% PLASMA magnetic field
%--------------------------------------
diag='difdop';
close all

if nargin>0
    numchoc = varargin{1};
else
    numchoc = input(' choc ');
end
disp(['Choc ' int2str(numchoc)]);


if ~exist('userdir','var'),userdir=getenv('HOME');end;
if ~exist('username'),username=getenv('LOGNAME');end;
if isempty(username),username=userdir(max(strfind(userdir,'/'))+1:end);
end

RayDir=[userdir,'/DIFDOP_LV/processed/TraceRay'];
FigureDir=[RayDir,'/',int2str(numchoc)];
if exist(FigureDir,'dir')~=7,eval(['mkdir ',FigureDir]);end

% required for experimental analysis:
% DIFDOP Frequency
% DIFDOP Angle pol
% DIFDOP angle tor
%   how much trigger times?

% here, from Sdifdop file
% which has already retrieved Top parameters, slow and fast data
% but can be done independantly (no need of fast data except times (date)
filemat=[userdir,'/DIFDOP_LV/data/Sdifdop',int2str(numchoc),'.mat'];
load(filemat);clear Skw* dtD* nfft;
% first: choose which diag and which polarisation:

if ~exist('Polarisation','var')
    if numchoc>41200 %campagne 2008 2 voies simultan�es mode O/X
        Polarisation=2;
    elseif numchoc>39090
        Polarisation=input('Polarisation (0/1)?[0]');
        if isempty(Polarisation),Polarisation=0;end;
    elseif numchoc<39090
        Polarisation=0;
    end
end
%if length(difdop_anglepol)~=nb_trig;error('Pb donnees: pas le meme nombre iteration et angle antenne');end

% choose which diag and which channel (O/X mode) to analyse


if nargin>1 
    if varargin{2}==3 % dreve
        diag='dreve'
        modex = 0
        nb_freq=length(F_pal_dreve)
        nb_trig=length(t0FA_dv)
        tdecl=reshape(t0Fp_dv',nb_freq,nb_trig)
        F_pal=F_pal_dreve;textmodex='VO';
    elseif Polarisation==2 % numchoc>41200 %campagne 2008 2 voies simultan�es mode O/X
        modex = varargin{2};
        diag='difdop';
        % calculates the number of Ray Tracing to be done
        % eventually proposes a choice on time (angle) then channel O/X?
        nb_trig=length(t0FA)
        if modex==0
            nb_freq=length(F_pal_O);tdecl=reshape(t0Fp',nb_freq,nb_trig);
            F_pal=F_pal_O;textmodex='_O';
        else
            nb_freq=length(F_pal_X);tdecl=reshape(t0Fp',nb_freq,nb_trig);
            F_pal=F_pal_X;textmodex='_X';
        end
    end
elseif nargin<2 && Polarisation==2
    modex=input('Analyse de quelle polarisation O/X [default 0->mode O]');
    if isempty(modex);modex=0;end
    if modex==0
        nb_trig=length(t0FA)
        nb_freq=length(F_pal_O);tdecl=reshape(t0Fp',nb_freq,nb_trig);
        F_pal=F_pal_O;textmodex='_O';
    else
        nb_trig=length(t0FA)
        nb_freq=length(F_pal_X);tdecl=reshape(t0Fp',nb_freq,nb_trig);
        F_pal=F_pal_X;textmodex='_X';
    end
elseif Polarisation==1 % campagne 2007 Polarisation choisie par switch
    diag='difdop';
    nb_freq=length(F_pal);tdecl=reshape(t0Fp',nb_freq,nb_trig);modex=1;textmodex='_X';
elseif Polarisation==0
    diag='difdop';
    nb_freq=length(F_pal);tdecl=reshape(t0Fp',nb_freq,nb_trig);modex=0;textmodex='_O';
end

disp(['Analyse du mode (mode O: 0 / mode X: 1) ->' int2str(modex)]);
disp(['pour ', diag]);


disp(sprintf('il y a %d declenches, aux temps  et angles',nb_trig))
disp(sprintf('      %2d            %1.3g\n',cat(2,[1:nb_trig]',t0FA)'))
choix_trigger=input(['quelles declenches ?[1-',num2str(nb_trig),']']);
%choix_trigger=1:nb_trig;

if numchoc>50000 & modex==0
    nb_freq=nb_freq_V;
elseif numchoc>50000 & modex==1
    nb_freq=nb_freq_W;
end

if question(['tous les paliers de frequence du faisceau sonde? [1-',num2str(nb_freq),']?'],'y')=='n'
            choix_freq=input(['choix des paliers de frequence du faisceau sonde? [1-',num2str(nb_freq),']?']);
        else
            choix_freq=1:nb_freq;
end
% corrections/calibration signal reconstruction
[difdop_phase,difdop_amplitude,phidirgoladeg]=difdop_get_cor(numchoc);
clear difdop_phase difdop_amplitude
%  % angle toroidal par rapport � la vis�e vers le centre du plasma
%  % compte tenu de l'indice du hublot: 2� � l'ext�rieur -> 0� dans la machine
%  % angle > 0 : boitier s'eloigne du diag thomson, faisceau vers le centre du queusot
%  phidirgoladeg d�fini suivant cahier exp�rimental
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLASMA equilibrium: analytic 

if numchoc>53000
fprintf(' Lecture des parametres plasma ... \n')   
[GE,t] = tsbase(numchoc,'GMAG_GEOM');

vaa=GE(:,3)/1000;
vb=GE(:,4)/1000; 
kappa=GE(:,5)/1000;GE

ve=kappa;
delta_u=GE(:,6);
delta_d=GE(:,7);
vtau=(delta_u+delta_d)./2000;

[BA,t]= tsbase(numchoc,'GMAG_BARY');
vRo=BA(:,1);
vZo=BA(:,2);

[SS,t]= tsbase(numchoc,'GMAG_SHAF');
dshaps=SS(:,3)/1000;
sd0=-SS(:,6)/1000;

[Itor,titor]= tsbase(numchoc,'GMAG_ITOR');
t=t(:,1);tq=t;
titor=titor(:,1);tdsha0=t;tp=t;dshaexp=1;
BR = (mean(Itor((titor>0) & (titor<max(t0Fp)))))*0.0073008;
RR=236;
plasma.B0=BR/RR*100;
out=imas_west_get(numchoc,'wall');
plasma.paroiz =out.description_2d{1}.limiter.unit{1}.outline.z;
plasma.paroix=out.description_2d{1}.limiter.unit{1}.outline.r;

else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLASMA equilibrium: analytic 
fprintf(' Lecture des parametres plasma ... \n')
[va,t] = tsbase(numchoc,'SAMIN');
[vaa,t] = tsbase(numchoc,'SAEQA');
[vRo,t] = tsbase(numchoc,'SRMAJ');
[vZo,t] = tsbase(numchoc,'SZPOS');
[ve,tq] = tsbase(numchoc,'SELLIP');
[vtau,tq] = tsbase(numchoc,'STRIANG');

[sd0,tdsha0] = tsbase(numchoc,'SD0MAG');
%  [dsha0,tdsha0] = tsbase(numchoc,'SD0MAG');
if isempty(sd0),
    [vsbeli,t] = tsbase(numchoc,'SBELI');
end;

dshaexp=1;
% dshaexp = entree(dshaexp,'Type de piquage du Shafranov : dshap=2 (0), dshap issu des donnees (1)',0);
plasma.dshaexp=dshaexp;
if dshaexp,
    [dshaps,tp] = tsbase(numchoc,'ssmag');
else
    dshap=2;
end

% Magnetic field
RR = 236;   % toroidal coil center for magnetic field calculation (cm)
[Itor,titor]=tsbase(numchoc,'STORI');
if isempty(Itor),[Itor,titor]=tsbase(numchoc, 'SITOR');end %avant ciel
BR = (mean(Itor((titor>0) & (titor<max(t0Fp)))))*0.0073008;
plasma.B0=BR/RR*100;
[plasma.paroiz,plasma.paroix] = tsbase(numchoc,'SPAROI');

end %if

%if length(tpo)>2
 %   plasma.paroix=mean(paroix);plasma.paroiz=mean(paroiz);
%else
 %   plasma.paroix=paroix;plasma.paroiz=paroiz;
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if X mode, use B with diamagnetic/paramagnetic corrections or not
if modex
    qmodeb=1;
    %qmodeb=input('use B with diamagnetic/paramagnetic corrections (2) or not (1) [default 1->without]');
    if qmodeb==2;
        modex=2;
        disp('mode X avec corrections para/dia')
        if plasma.B0~=0
            [vIp,tip] = tsmat(numchoc,'DPOLO-SIPMES');
            [vj0,tj] = tsbase(numchoc,'SJCENTRE');% en MA/m2: pb difference avec Ip/pi/a^2 *(v+1) ??
            [vpiqj,tj]= tsbase(numchoc,'SJPIQ');
              if isempty(vj0),
                  [vsdiam,tsd] = tsbase(numchoc,'SDIAM');
              end;
        else
        Ip = entree(Ip,'courant plasma ?',1);
        Te0 = entree(Te0,'Te0 ?',3);
        end
    end
end

% which density profile
% 1/analytic representation from experimental central density (interferometry)
%            and peaking factor (from <ne> and nec): edge profile not well known
% 2/from TPROF but only for rho<1, still pb in the SOL
% 3/from reflectometry...
% evaluate at all acquisition time
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
time.t0FA=t0FA(choix_trigger);
time.t0Fp=reshape(tdecl(choix_freq,choix_trigger),length(choix_trigger)*length(choix_freq),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin>3
    filene=varargin{4};
    [ne_struct]=difdop_ne4beam_b(numchoc,time,type_profile,filene);
    disp('coucou')
else
    [ne_struct]=difdop_ne4beam_b(numchoc,time,type_profile);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(ne_struct.ne)
    disp('beurk')
    return
end
save test ne_struct

nb_times_ne=length(ne_struct.times);
prof_type=ne_struct.prof_type
if prof_type=='T';ocur=ne_struct.tprof;end
%neprofile.rhofit=ne_struct.rho;
if prof_type=='T';nelambdad=ne_struct.nelambda;end
%krhoa=find(neprofile.rhofit==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% density profile %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop on all times (angles + frequencies) selected
for k_trigger=1:length(choix_trigger)
    angle_pol=difdop_anglepol(choix_trigger(k_trigger)) ;        
    Freq=F_pal(choix_freq);
    
    for k_freq=1:length(choix_freq)
        % experimental -> now averages plasma parameters at each texp (each beam frequency step)

        % % minor and major radius
        if isempty(ve),
            plasma.e=1;
            plasma.a0r = 100*mean(va(t>=tdecl(choix_freq(k_freq),choix_trigger(k_trigger)) & ...
                                 t<=tdecl(choix_freq(k_freq),choix_trigger(k_trigger))+0.0082));
        else
            plasma.e = mean(ve(tq>=tdecl(choix_freq(k_freq),choix_trigger(k_trigger)) &...
                           tq<=tdecl(choix_freq(k_freq),choix_trigger(k_trigger))+0.0082));
            plasma.a0r = 100*mean(vaa(t>=tdecl(choix_freq(k_freq),choix_trigger(k_trigger)) & ...
                                 t<=tdecl(choix_freq(k_freq),choix_trigger(k_trigger))+0.0082));
        end;
        plasma.rr0 = 100*mean(vRo(t>=tdecl(choix_freq(k_freq),choix_trigger(k_trigger)) & ...
                                  t<=tdecl(choix_freq(k_freq),choix_trigger(k_trigger))+0.0082));
        plasma.zpos = 100*mean(vZo(t>=tdecl(choix_freq(k_freq),choix_trigger(k_trigger)) & ...
                                   t<=tdecl(choix_freq(k_freq),choix_trigger(k_trigger))+0.0082));
        if isempty(ve),plasma.tau=0;else
        plasma.tau = mean(vtau(tq>=tdecl(choix_freq(k_freq),choix_trigger(k_trigger)) &...
                                   tq<=tdecl(choix_freq(k_freq),choix_trigger(k_trigger))+0.0082));end
        % Shafranov shift
        if isempty(sd0),
            sbeli = mean(vsbeli(t>=tdecl(choix_freq(k_freq),choix_trigger(k_trigger)) &t<=tdecl(choix_freq(k_freq),choix_trigger(k_trigger))+0.0082));
            plasma.dsha0 = plasma.a0r^2/2/rr0*sbeli*(1-plasma.a0r/rr0);
        else
            plasma.dsha0 = 100*mean(sd0(tdsha0>=tdecl(choix_freq(k_freq),choix_trigger(k_trigger)) &tdsha0<=tdecl(choix_freq(k_freq),choix_trigger(k_trigger))+0.0082));
        end;
        if isempty(plasma.dsha0)
            plasma.dsha0 = 1;
        end;
        if dshaexp,
            plasma.dshap = mean(dshaps(tp>=tdecl(choix_freq(k_freq),choix_trigger(k_trigger)) &tp<=tdecl(choix_freq(k_freq),choix_trigger(k_trigger))+0.0082));
        else
            plasma.dshap = 2;
        end
        
        if modex==2,
            if plasma.B0~=0
                if type_profile==4|5
                    plasma.j0 = mean(vj0(tj>=tdecl(choix_freq(k_freq),choix_trigger(k_trigger)) &tj<=tdecl(choix_freq(k_freq),choix_trigger(k_trigger))+0.0082))*1e6;
                    plasma.jp = mean(vpiqj(tj>=tdecl(choix_freq(k_freq),choix_trigger(k_trigger)) &tj<=tdecl(choix_freq(k_freq),choix_trigger(k_trigger))+0.0082));
                    plasma.Ip = mean(vIp(tip>=tdecl(choix_freq(k_freq),choix_trigger(k_trigger)) &tip<=tdecl(choix_freq(k_freq),choix_trigger(k_trigger))+0.0082));
                    if isempty(plasma.j0)
                        %   pij et j0 analytique (Eq pour aspect petit, Wesson)
                        %   valable pour des profils monotones parabol+piq
                        %     manque betapol :SDIAM ou SFDIAM
                        sdiam = mean(vsdiam(tsd>=tdecl(choix_freq(k_freq),choix_trigger(k_trigger)) &tsd<=tdecl(choix_freq(k_freq),choix_trigger(k_trigger))+0.0082));
                        if ~isempty(sdiam)
                            li=2*(sbeli-sdiam);
                            plasma.jp=(1/0.89)*(exp(li)-1.65);
                        else
                            plasma.jp=2;
                        end
                        plasma.j0=(plasma.Ip*1e6./pi./(plasma.a0r.^2)).*(1+ plasma.jp );
                    end;
                else
                    plasma.jp=2;
                    plasma.j0=(plasma.Ip*1e6./pi./(plasma.a0r.^2)).*(1+ plasma.jp );
                end
            end
        end
        % structure equilibre        plasma.a0r=a0r...

        % Density Profile        
         %if nb_times_ne==length(time.t0Fp)   neprofile.ne=ne_struct.ne; 
           %  neprofile.ne=ne_struct.ne((k_trigger-1)*length(choix_freq)+k_freq,:);   
            neprofile.ne=ne_struct.ne; 
            neprofile.rhofit=ne_struct.rho;
           % krhoa=find(neprofile.rhofit==1);
      % elseif nb_times_ne==length(time.t0FA) % cas avec reflec?
           % neprofile.ne=ne_struct.ne(k_trigger,:);
           % if size(ne_struct.rho,1)>size(neprofile.ne,1), 
            %    neprofile.rhofit=ne_struct.rho(k_trigger,:);
             %   krhoa=iround(neprofile.rhofit,1);
           % else
            %    neprofile.rhofit=ne_struct.rho;
             %   krhoa=find(neprofile.rhofit==1);
           % end
       % else
          %  neprofile.ne=ne_struct.ne;
          %  neprofile.rhofit=ne_struct.rho;
            krhoa=find(neprofile.rhofit==1);
     %  end
        % finally structure ne_struct with ne rhofit
       
        % Display which appears on figure 
         titrechoc = sprintf('Beam Tracing ');
         fprintf('\n%s\n',titrechoc);
    
         textexp=sprintf('%g_%g',fix(tdecl(choix_freq(k_freq),choix_trigger(k_trigger))),round(10*(tdecl(choix_freq(k_freq),choix_trigger(k_trigger))-fix(tdecl(choix_freq(k_freq),choix_trigger(k_trigger))))));

         textshot = sprintf('TS#%d @t=%.1f s\n',numchoc,tdecl(choix_freq(k_freq),choix_trigger(k_trigger)));
         fprintf('%s\n',textshot)

         if prof_type=='T'
             textpla = sprintf('ne TPROF oc-%d\n ne(0)=%.1f 10^{19}m^{-3}\n',ocur,neprofile.ne(1)/1e19);
             textpla = sprintf('%s ne_a=%.1f 10^{19}m^{-3}\n',textpla,neprofile.ne(krhoa)/1e19);
             textpla = sprintf('%s \\lambda_{ne}=%.0fcm  ',textpla,nelambdad);
         elseif prof_type=='A'
             textpla = sprintf('ne Analytique\n ne(0)=%.1f 10^{19}m^{-3}\n',neprofile.ne(1)/1e19);
             textpla = sprintf('%s ne_a =%.2f 10^{19}m^{-3}\n',textpla,neprofile.ne(krhoa)/1e19);
             textpla = sprintf('%s \\lambda_{ne}=%.0fcm  ',textpla,nelambdad);
         elseif prof_type=='I'
             textpla = sprintf('ne Interferometry\n ne(0)=%.1f 10^{19}m^{-3}\n',neprofile.ne(1)/1e19);
             textpla = sprintf('%s ne_a =%.2f 10^{19}m^{-3}\n',textpla,neprofile.ne(krhoa)/1e19);
             textpla = sprintf('%s \\lambda_{ne}=%.0fcm  ',textpla,nelambdad);
         elseif prof_type=='EXE'
             textpla = sprintf('ne Reflec EFIT \n ne(0)=%.1f 10^{19}m^{-3}\n',neprofile.ne(1)/1e19);
             textpla = sprintf('%s ne_a =%.2f 10^{19}m^{-3}\n',textpla,neprofile.ne(krhoa)/1e19);
             %textpla = sprintf('%s \\lambda_{ne}=%.0fcm  ',textpla,nelambdad);
         elseif prof_type=='EXT'
             textpla = sprintf('ne Reflec TMAG \n ne(0)=%.1f 10^{19}m^{-3}\n',neprofile.ne(1)/1e19);
             textpla = sprintf('%s ne_a =%.2f 10^{19}m^{-3}\n',textpla,neprofile.ne(krhoa)/1e19);
             %textpla = sprintf('%s \\lambda_{ne}=%.0fcm  ',textpla,nelambdad);
         end
         fprintf('%s\n',textpla)

         save test_jab plasma neprofile angle_pol phidirgoladeg
         %%%%%%%% trac� rayon
         [dif,opt]=beam2010(diag,Freq(k_freq),angle_pol,modex,phidirgoladeg,neprofile,plasma);
         %%%%%%%% trac� rayon
   
         ripple = dif.ripple;
         %global thetadirgola -> on se sert de angle_pol
         nbrayr = dif.nbrayr;
         nbraya = dif.nbraya;
         waist = dif.waist;
       
         rhodif(k_freq,:,:)= dif.rhodif; nx2dif(k_freq,:,:)= dif.nx2dif; 
         xdif(k_freq,:,:)= dif.xdif;ydif(k_freq,:,:)= dif.ydif;zdif(k_freq,:,:)= dif.zdif;
         kdif(k_freq,:,:)= dif.kdif; kxdif(k_freq,:,:)= dif.kxdif; kydif(k_freq,:,:)= dif.kydif; kzdif(k_freq,:,:)= dif.kzdif; 
         krdif(k_freq,:,:)= dif.krdif; kthetadif(k_freq,:,:)= dif.kthetadif; kphidif(k_freq,:,:)= dif.kphidif;
         rko =opt.rko; hmin=opt.hmin; hmax=opt.hmax; tol=opt.tol; testklim=opt.testklim;
         
         subplot(211)
         text(455,120,textshot,'HorizontalAlignment','right','VerticalAlignment','top');
%         text(455,100,textpla,'HorizontalAlignment','right','VerticalAlignment','top');
         subplot(212)
         text(150,100,textpla,'HorizontalAlignment','left','VerticalAlignment','bottom');
         %%%%% figure densite
         h=axes('position',[0.5  0.42 .36 .12]);
         %set(h,'visible','off')
         plot(neprofile.rhofit,neprofile.ne);%,'o',rho,nerho(rho))
         set(gca,'xlim',[0 1.8]);
         xlabel('rho')
         ylabel('ne')
         pause(0.5)
    
    end
    file=[RayDir,'/',num2str(numchoc),'RKdif',textmodex,...
        '_',prof_type,'_',num2str(choix_trigger(k_trigger)),'_',username(1:2)]

    save(file,'numchoc','t0FA*','angle_pol',...
          'Freq', 'phidirgoladeg','Polarisation','modex',...
          'prof_type','neprofile',...
          'nbrayr','nbraya','waist','ripple',...
          'rhodif','xdif','ydif','zdif',...
          'kxdif','kydif','kzdif','krdif','kthetadif','kphidif','kdif')

      figure(1);orient tall
    FigureName=[FigureDir,'/',num2str(numchoc),'Beam',textmodex,...
        '_',prof_type,'_',num2str(choix_trigger(k_trigger)),'_',username(1:2)];
    eval(['print -dpsc2 ',FigureName]);

    if nb_freq>1
        for k_pal=2:length(choix_freq)
            figure(k_pal);orient tall
            eval(['print -dpsc2 -append ',FigureName]);
            
        end
    end
    
    eval(['!ps2pdf ',FigureName,'.ps ',FigureName,'.pdf']);
    eval(['delete ',FigureName,'.ps']);

  
    if length(choix_trigger)~=1
        close all
    end
    
end