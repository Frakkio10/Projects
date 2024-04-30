% difdop.m
%
% Controle des declenches et signaux 
% Affichage du scenarios et des acquis lentes
% Calcul et affichage des spectres par declenche
% mode paliers egaux
%      paliers programmes
%      acqui longue F continue
% novembre 2005 PH
% Modif angles 2006
% modif nouvelle acqui 10/2007
 
clear all ;%close all;
if ~exist('username'),username=getenv('LOGNAME');end;
if ~exist('userdir'),userdir=getenv('HOME');end;
dateplot=date;

numchoc=input(' Numero du choc : ');
[date,heure]=tsdate(numchoc);

DirPlot=[userdir,'doppler/dataplot/'];
if exist(DirPlot,'dir') ~=7;DirPlot='./';end;

epsi0 = 8.8541e-12;
me = 9.1e-31;
qe = 1.6e-19;
%<> ~

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% controle des acquis lentes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[svideo,tlent,srms,tlent,sposverin,tlent,sfreq,tlent] = tsmat(numchoc,...
     'DIFDOP-SRMSVIDEO','DIFDOP-SRMSCOS','DIFDOP-SPOSVERIN','DIFDOP-SFREQ');

 % en volts
sposverin=sposverin*10/2048;
srms=srms*10/2048;
svideo=svideo*.014-8.3;  % x 10 dBm calibration AT
%sfreq=sfreq*10/2048;
% le signal sfreq est depuis 2006 devolu a l'inclinometre
sfreq=asin(sfreq/2048)*180/pi-7.5;% inclinometre asin(SV/10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  calibration signal verin en angle

if numchoc > 39100 %campagne 2007 (calibr CH inclinom?tre)
  sangle=sposverin.^2*0.00785+sposverin*1.09 + 4.64; % 2007
  %%%%% pb de perte de reference (zero) du verin depuis automne 2007
  if numchoc > 40050
      sangle1=sangle+7;
  end;
elseif numchoc > 34300 %campagne 2005 (mesure 30/05/2005 0.01966 1.18296 7.42788)
 %sangle=sposverin.^2*0.0196+sposverin*1.1829 +7.4279;
  %elseif numchoc > 32330 %campagne 2004
  %sangle=(sposverin+13.3157)/.9263 -7.5;
  %sangle1=sposverin*1.0773+6.87;
  sangle=sposverin.^2*0.0115+sposverin*1.1156 +6.7398;
elseif numchoc > 32040 %campagne 2003, novembre
  sangle=(sposverin+2.3427)/.3365 -7.5;
else
  sangle=[];%pas de calib ou de signal moteur
end

%%%%%%%%%%
% Affichage du Sc?nario du choc
%if (exist('svideo'))&(~isempty(svideo))
disp('lecture de Itor, nl et Ip')
[Itor,titor,ignitron] = tsmat(numchoc,'DMAG-STORI','IGNITRON|1');
if isempty(ignitron),ignitron=32;end;
[nl,tnl,Ip,tIp] = tsmat(numchoc,'DGAZ-SNLI','DPOLO-SIPMES');
tmin = min(tIp(Ip > 0.95*max(Ip)));
tmax = max(tIp(Ip > 0.95*max(Ip)));
[SAMIN,ta,ca,SQPSI,tq,cq,SRMAJ,tR,cR,SZPOS,tZ,cZ]=...
            tsbase(numchoc,'SAMIN','SQPSI','SRMAJ','SZPOS');

if ~isempty(Itor)&~isempty(tmax)~isempty(tmin);
Ipo = (mean(Ip((tIp>tmin) & (tIp<tmax))));
a = (mean(SAMIN((tq>tmin) & (tq<tmax))));
qpsi = (mean(SQPSI((tq>tmin) & (tq<tmax))));
Ro = (mean(SRMAJ((tq>tmin) & (tq<tmax))));
Bo = (mean(Itor((titor>tmin) & (titor<tmax))))*0.0073/2.36;
fce0=qe*Bo/me/2/pi/1e9;fcea=fce0*Ro/(Ro+a);
else
a=[];Ro=[];Bo=[];qpsi=[];fce0=[];fcea=[];
end

disp('lecture des divers (chauffages/divertor/glacons)')
[nmoy,tnmoy,cr] = tsbase(numchoc,'SNMOY');
[nepiq,tgne,cpiq]=tsbase(numchoc,'SPIQ');
[GNE,tgne,rhone,cpn]=tsbase(numchoc,'GNE');
[PFCE1,tFCE]=tsbase(numchoc,'SPIA1');
[PFCE2,tFCE]=tsbase(numchoc,'SPIA2');
PFCE=PFCE2+PFCE1;clear PFCE1 PFCE2
[hyb,thyb] = tsmat(numchoc,'DHYB-GPHYB%3');
[fci,tfci] = tsbase(numchoc,'gbilan%3');%'DFCI-SPUISS' ???
[idn,tidn] = tsbase(numchoc,'SIDN');

[smhd,tmhd]=tsbase(numchoc,'SMHD');
duree_ip=tsbase(numchoc,'RDUREE');
if isempty(duree_ip);duree_ip=tmax+3;end; % au lieu de 15
if isempty(thyb)
   [hyb,thyb] = tsmat(numchoc,'DHYB-GHYB%3'); 
end

% ECE traite
%if 0
[Tecenv,ttecenv,voies,certdon]=tsbase(numchoc,'GSHTENV');
[Tece,ttece,voies,certdon]=tsbase(numchoc,'GSHTE');
if isempty(ttece)
   disp('No slow ECE signals') 
   ttece=ttecenv;
end
%  [cert]=tsbase_cert(certdon);
[Rece,ttece]=tsbase(numchoc,'GSHR'); 
Receo = mean(Rece((ttece>tmin) & (ttece<tmax),:));
kece=find(abs(Receo-2.44)==min(abs(Receo-2.44)));
%  [Tecenv]=tsbase(numchoc,'GSHTENV');
if length(Tece)>0
   Nvoies=length(Tece(1,:));
else
  Nvoies=0.;
end
%end % Te

figure(1);clf;grid on;
subplot(211)
plot(tIp,Ip,'k',tnl,nl,'b',tfci,fci,'r',thyb,hyb,'m',tq,SQPSI,'y')
     if ~isempty(GNE);hold on;plot(tgne,GNE(:,1)/1e19,'c');end;
     if ~isempty(idn);hold on;plot(tidn,idn,'r:');end;
%     if exist('Tece','var')==1&~isempty(Tece);hold on;plot(ttece,Tece(:,kece),'m--');end;
set(gca,'xlim',[-1 duree_ip]);grid on;
set(gca,'ylim',[0 min(max([max(nl) max(fci) max(Ip)]),10)]);
title(sprintf('TS# %g, B_0=%5.2f T, Ip=%5.2f MA, %s-%s',...
                   numchoc,Bo,Ipo, date,heure));

subplot(212)
plot(tlent,svideo,tlent,srms,tlent,sfreq);hold on;

set(gca,'xlim',[-1 duree_ip]);
  

%%%%%%%%%%%%%%%%%%%%%%
TDIFDOP(1)=tsmat(numchoc,['TDIFDOP|',int2str(1)]);
[T,TRIG,Cr]=tsdeclenche(numchoc);
TDIFDOP=T(tsfindecl(TRIG,'TDIFDOP'))/1e3;
%%%%%%%%%%%%%%%%%%%%%%%
% temps des declenches TDIFDOP (chrono TS) pour comparaison avec vecteur temps
if numchoc>32305 & numchoc <32945;
  for id=1:iter_nb
   TDIFDOP(id)=tsmat(numchoc,['TDIFDOP|',int2str(id)]);
  end
elseif numchoc>32945% sept2004: declenches angles gerees par LW110
  [T,TRIG,Cr]=tsdeclenche(numchoc);
  TDIFDOP=T(tsfindecl(TRIG,'TDIFDOP'))/1e3;
end

nb_top=length(TDIFDOP); % nombre de declenche programmes dans TOP
%nb_top=1;
disp('trigger TDIFDOP-ignitron @ ')
disp(sprintf('%g\n',TDIFDOP-ignitron))
% sometimes wrong, indicates several TDIFDOP with less than 100ms
if nb_top>=2;
    for itop=2:nb_top;
        if TDIFDOP(itop)-TDIFDOP(itop-1)<.1,TDIFDOP(itop)=[];
        end;
    end
 else
   itop=1;
end

if itop==1;
TDIFDOP(itop)=TDIFDOP;
end
%nb_top=length(TDIFDOP); % nombre de declenche programmes dans TOP


%%%%%%%%%%%%%%%%
% Lecture des acquis rapides
%%%%%%%%%%%%%%%
if tsexist(numchoc,'FDATEDOP') == 0 | tsexist(numchoc,'FACQIRIS_DAT') == 0
    disp('OLD CASES')
	%%%%%%%%%%%%%%%%
	% Lecture des parametres DIFDOP
	%%%%%%%%%%%%%%%
    %campagne 2007 switch polarisation modeO/X
	if numchoc<39090;
		Polarisation=0;
	elseif numchoc>39090&numchoc<41200;
		Polarisation=tsmat(numchoc,'DIFDOP;RABOTEUR;Polarisation');
	elseif numchoc>41200
		Polarisation=2;
	end;

    % lecture des variables DEBOG : quelle source? verin ok?
    debog_lw=tsmat(numchoc,'DIFDOP;DEBOG;debog_lw');
    debog_c3=tsmat(numchoc,'DIFDOP;DEBOG;debog_c3');
    debog_mg=tsmat(numchoc,'DIFDOP;DEBOG;debog_mg');

    % PILOTAGE: VERIN
	Angle_Ini=tsmat(numchoc,'DIFDOP;VERIN;Angle_Ini');
	NbDecl_verin=tsmat(numchoc,'DIFDOP;VERIN;NbDecl_verin');
	Angle_Decl=tsmat(numchoc,'DIFDOP;VERIN;Angle_Decl');
	VitAng_Decl=tsmat(numchoc,'DIFDOP;VERIN;VitAng_Decl');
	%Angle_Top=[Angle_Ini Angle_Decl(1:NbDecl_verin)];
	Angle_Top=[Angle_Ini Angle_Decl(1:nb_top)];

%%%%%%%%%%
% correction du pb sur les angles ??? partir de spos verin: perte du zero
    if numchoc > 40050
         sangle_ref=zeros(size(Angle_Top));spos_ref=sangle_ref;
         diff_Angle_Top=diff(Angle_Top);
%         if nb_top==1&diff_Angle_Top==0, sangle=Angle_Top(1)*ones(size(sposverin));
         if diff_Angle_Top==0|NbDecl_verin==0, sangle=Angle_Top(1)*ones(size(sposverin));
         else 
         
            t_ref_Ini=[1 TDIFDOP'+diff_Angle_Top./VitAng_Decl(1:nb_top)-ignitron+0.2]
            t_ref_Fin=[TDIFDOP'-0.2-ignitron duree_ip-1]
                for itop=1:nb_top+1;
                sangle_ref(itop)=mean(sfreq(tlent>t_ref_Ini(itop)&tlent<t_ref_Fin(itop)));% valeur entre d???clenches
                disp(['valeur Angle_Top / Angle inclinom???tre ',num2str(Angle_Top(itop)),' / ',num2str(sangle_ref(itop))])
                spos_ref(itop)=mean(sangle1(tlent>t_ref_Ini(itop)&tlent<t_ref_Fin(itop)));% valeur spos entre d???clenches
                end
         
                if nb_top>1
                p=polyfit(spos_ref,sangle_ref,2);
                sangle=p(3)+p(2)*sangle1+p(1)*sangle1.^2;
                else
                    if diff(spos_ref)<=0.25 | diff(sangle_ref)>=0.3
                    sangle=mean(sangle_ref)/mean(spos_ref)*sangle1;
                    else
                    p=polyfit(spos_ref,sangle_ref,1);
                    sangle=p(2)+p(1)*sangle1;
                    end       
                end
         end
         plot(tlent,sangle,'k--')
    end %if numchoc > 40050
     
    % PILOTAGE: source: ANRITSU
	% strategie de pilotage de la frequence sonde
	if (numchoc >32330) % apr???s campagne 2004
        mod_egaux=tsmat(numchoc,'DIFDOP;ANRITSU;mod_egaux');
        nb_freq=tsmat(numchoc,'DIFDOP;ANRITSU:nb_freq');
        F_start=tsmat(numchoc,'DIFDOP;ANRITSU;F_start');
        F_end=tsmat(numchoc,'DIFDOP;ANRITSU;F_end');
        duree_pal=tsmat(numchoc,'DIFDOP;ANRITSU;duree_pal');
        mod_inegaux=tsmat(numchoc,'DIFDOP;ANRITSU;mod_inegaux');
        F_pal=tsmat(numchoc,'DIFDOP;ANRITSU;F_pal');
        mod_CW=tsmat(numchoc,'DIFDOP;ANRITSU;mod_cw');
        F_CW=tsmat(numchoc,'DIFDOP;ANRITSU;F_CW');
        mod_balay=tsmat(numchoc,'DIFDOP;ANRITSU;mod_balay');
        Prof_start=tsmat(numchoc,'DIFDOP;ANRITSU;Prof_start');
        Prof_end=tsmat(numchoc,'DIFDOP;ANRITSU;Prof_end');
        duree_rampe=tsmat(numchoc,'DIFDOP;ANRITSU;duree_rampe');
        puissance=tsmat(numchoc,'DIFDOP;ANRITSU;puissance');
    end %if (numchoc >32330) 


	%%%%%%%%%%%%%%%%%%%%%
	% Acquis carte ICV122 jusqu'a campagne automne 2007
	%%%%%%%%%%%%%%%%%%%%%
	if numchoc<40050
    status = tsrfile (numchoc,'FDATEDOP',['FAtime.dat']);
		%if status==0 %le fichier est cr??, test sur le temps
    fid1=fopen(['FAtime.dat'],'r');
    [Ftime]=fscanf(fid1,'%d');
            if ~isempty(Ftime)
        	status = tsrfile (numchoc,'FDOP1',['Fdop1.dat']);
            fid2=fopen(['Fdop1.dat'],'r','ieee-be');
            [v1]=fread(fid2,'short');
            status = tsrfile (numchoc,'FDOP2',['Fdop2.dat']);
            fid3=fopen(['Fdop2.dat'],'r','ieee-be');
            [v2]=fread(fid3,'short');
            end;
  	fclose all;

        % Strategie d'acquisition et de param acqui rapide ICV122
    Mode_acqui=tsmat(numchoc,'DIFDOP;RABOTEUR;Mode_acqui');
    Retard_Rab=tsmat(numchoc,'DIFDOP;RABOTEUR;Retard_Rab');
    NbCoups_Rab=tsmat(numchoc,'DIFDOP;RABOTEUR;NbCoups_Rab');
        
    taille_seg=NbCoups_Rab;

    repsd =tsmat(numchoc,'DIFDOP;MDIFDOP2;repsd');
	maxram=tsmat(numchoc,'DIFDOP;MDIFDOP2;maxram');
	nbvoies=tsmat(numchoc,'DIFDOP;MDIFDOP2;nbvoies');
	datation=tsmat(numchoc,'DIFDOP;MDIFDOP2;datation');
	nbmotifs=tsmat(numchoc,'DIFDOP;MDIFDOP2;nbmotifs');
        

		% vecteur Ftime: 1 date tous les repsd pts (default 512)
		% soit NbCoups_Rab/repsd (8192/512) pts de date par paliers 
		
		% sachant que repsd est pris par l'horloge interne de TS
		% sauf si datation=2, horloge interne 8MHz
		% soit Ftime/1e6= temps absolu en s
		% ou Ftime/8e6 temps a partir de UNE?

	dthorloge=min(diff(Ftime));
	FT=Ftime;
            if datation==2
			Ftime=Ftime/8e6;
            elseif(datation==1&dthorloge/repsd>1);
			textpb='probable probleme de datation, - 8MHz';
			disp('probable probleme de datation, - 8MHz');
			Ftime=Ftime/8e6;
            elseif datation==1
			Ftime=Ftime/1e6;textpb=' ';
            end

		% calcul du temps (indice) de debut de chaque enregistrement 
        % diff(Ftime/1e6) = repsd * .250e-6 (horloge raboteur)
		% sauf lorqu'on change de palier ou d'iteration
	kDecl=find(diff(Ftime)>2*repsd*250e-9); %2 pour les arrondis
        
        %old stuff
        %if Mode_acqui==1%&mod_egaux==1  % paliers
            %disp(' Mode acquisition 1, paliers de frequence sonde')
            
            %kiter=find(diff(Ftime(kDecl))>2*min(diff(Ftime(kDecl))));% les iterations
            %iter_ev=(length(kDecl)+1)/nb_freq;% ou bien length(kiter);
            % les temps de chaque iteration (angle)
           % t0FA=Ftime(kDecl(kiter)+1) ;t0FA=[Ftime(1);t0FA]-ignitron;

            %if debog_lw==1000 | numchoc <32930
             %   iter_nb=iter_ev;
            %else
    nb_iter=tsmat(numchoc,'DIFDOP;VCO_LW110:iter_nb');
              
                % nombre d'angles sur lesquels les paliers sont effectu?s
        if length(Ftime)*nb_iter
                    % coherence avec le nombre de declenche programmees TDIFDOP pilote
                    % TDIFDOP*iter_nb==iter_ev
        textDecl=sprintf('%d declenches TDIFDOP',length(TDIFDOP));
        disp(textDecl);
                    %iter_nb=iter_ev;
        else
        textDecl=sprintf('Pb nb declenches');
        disp(textDecl);
                    %iter_nb=iter_ev;
        end
            %end


            % les temps de chaque paliers
    t0Fp=Ftime(kDecl+1);t0Fp=[Ftime(1);t0Fp]-ignitron;
    if mod_egaux
       if F_end~=0
       F_end=min(F_end,18.75*4);
       F_start=max(F_start,12.25*4);
       F_pal=3/2*((0:nb_freq-1)*(F_end-F_start)/(nb_freq-1)+F_start);
       elseif F_end==0
       AGAZ_nl=tsmat(numchoc,'AGAZ;Nl;VALEURS');
       maxnlprog=max(AGAZ_nl(:,2))*1e19;
       F_max=floor(1.05*sqrt(0.90*.95*maxnlprog*qe^2/me/epsi0)/2/pi/1e9);
       F_max=max(F_max,60);F_max=min(F_max,75);
       F_pal=(0:nb_freq-1)*(F_max-F_start)/(nb_freq-1)+F_start;
       end
    end
    
        %end
        
    dtDifdop=.25e-6;
        
	%%%%%%%%%%%%%%%%
	% Acquis rapide carte Acqiris automne 2007

    elseif numchoc>40051
  	status = tsrfile(numchoc,'FACQIRIS_DAT','Fdate.gz');
  	!gzip -df *.gz
    fid1=fopen(['Fdate'],'r','ieee-le');
    [Ftime]=fread(fid1,'ulong');
            if ~isempty(Ftime)
            status = tsrfile(numchoc,'FACQIRIS0_V1','acq0_v1.gz');
            !gzip -fd acq0*.gz
            fid2=fopen('acq0_v1','r','ieee-le');
            [v1]=fread(fid2,'short');
            status = tsrfile(numchoc,'FACQIRIS0_V2','acq0_v2.gz');
            !gzip -fd acq0*.gz
            fid3=fopen('acq0_v2','r','ieee-le');
            [v2]=fread(fid3,'short');
            % 2 separate channels (mode O/X) from 2008d
                if numchoc>41200
                status = tsrfile(numchoc,'FACQIRIS1_V1','acq1_v1.gz');
                !gzip -fd acq1*.gz
                fid4=fopen('acq1_v1','r','ieee-le');
                [v3]=fread(fid4,'short');
                status = tsrfile(numchoc,'FACQIRIS1_V2','acq1_v2.gz');
                !gzip -fd acq1*.gz
                fid5=fopen('acq1_v2','r','ieee-le');
                [v4]=fread(fid5,'short');
                end;
            fclose all;
            
            ACQIRIS_0_freq = tsmat(numchoc,'DIFDOP;ACQIRIS_0;FREQ');  %freq d'acquisition en Hz
            ACQIRIS_0_tailleseg = tsmat(numchoc,'DIFDOP;ACQIRIS_0;TAILLE_SEG');
            ACQIRIS_0_nbseg = tsmat(numchoc,'DIFDOP;ACQIRIS_0;NB_SEG');
            ACQIRIS_0_umaxv1 = tsmat(numchoc,'DIFDOP;ACQIRIS_0;U_MAX_V1');  % [mV]
            ACQIRIS_0_umaxv2 = tsmat(numchoc,'DIFDOP;ACQIRIS_0;U_MAX_V2');
            TRIG_SOURCE=tsmat(numchoc,'DIFDOP;ACQIRIS_0;TRIG_SOURCE')
            taille_seg=ACQIRIS_0_tailleseg;
            dtDifdop=1/ACQIRIS_0_freq;
                if numchoc>46806
                taille_seg_W=ACQIRIS_0_tailleseg;
                dtDifdop_W=1/ACQIRIS_0_freq; 
                elseif numchoc>41200
                ACQIRIS_1_freq = tsmat(numchoc,'DIFDOP;ACQIRIS_1;FREQ');  %freq d'acquisition en Hz
                ACQIRIS_1_tailleseg = tsmat(numchoc,'DIFDOP;ACQIRIS_1;TAILLE_SEG');
                ACQIRIS_1_nbseg = tsmat(numchoc,'DIFDOP;ACQIRIS_1;NB_SEG');
                ACQIRIS_1_umaxv1 = tsmat(numchoc,'DIFDOP;ACQIRIS_1;U_MAX_V1');  % [mV]
                ACQIRIS_1_umaxv2 = tsmat(numchoc,'DIFDOP;ACQIRIS_1;U_MAX_V1');
                taille_seg_W=ACQIRIS_1_tailleseg;
                dtDifdop_W=1/ACQIRIS_1_freq;
                end
          

            
            Ftime=Ftime/1e6;textpb=' ';
            % le vecteur temps est ici constitue des temps de chaque declenche (iteration)
            kDecl=1:length(Ftime)-1;
            
            %nb_freq =ACQIRIS_0_nbseg; %verifier que nb_freq = nb_seg pour chacune des cartes
            
                if mod_egaux==1 | mod_inegaux ==1% paliers
                disp('mode paliers')
                kiter=find(diff(Ftime(kDecl))>2*min(diff(Ftime(kDecl))));% les iterations
                iter_ev=(length(kDecl)+1)/nb_freq;% ou bien length(kiter);
                % les temps de chaque angle
                t0FA=Ftime(kiter+1) ;t0FA=[Ftime(1);t0FA]-ignitron;
                % les temps de chaque palier
                t0Fp=Ftime(kDecl+1);t0Fp=[Ftime(1);t0Fp]-ignitron;
                
                nb_iter=tsmat(numchoc,'DIFDOP;VCO_LW110:iter_nb');
       
                % nombre d'angles sur lesquels la sequence des paliers est declenchee
                
                    if length(v1)==nb_iter*length(TDIFDOP)*nb_freq*taille_seg
                    % coherence avec le nombre de declenche programmees TDIFDOP pilote
                    % TDIFDOP*iter_nb==iter_ev
                	textDecl=sprintf('%d declenches TDIFDOP',length(TDIFDOP));
                	disp(textDecl);
                	%iter_nb=iter_ev;
                    else
                	textDecl=sprintf('Pb nombre de declenches');
                	disp(textDecl);
                	%iter_nb=iter_ev;
                    end

                end  %if mod_egaux==1 | mod_inegaux ==1% paliers

            end %  if ~isempty(Ftime) lecture des acquis rapides
        %%%%%% affichage
    figure(1)
    subplot(211)
    couleur=['b' 'g' 'c' 'r' 'm' 'k' 'y' 'b' 'g' 'c' 'r' 'm'];
    plot([t0FA t0FA]',[zeros(length(t0FA),1) 10*ones(length(t0FA),1)]','b--')
    subplot(212)
    title([sprintf('TS# %g, Polarisation %g, Angles =',numchoc,Polarisation),sprintf('%g-',Angle_Top)]);
    xlabel('time [s]')
    qfig=input('imprimer ?','s');
        if isempty(qfig),qfig='y';end;
        if qfig=='y'
        orient tall
        fileprintS=[DirPlot,num2str(numchoc),'S'];
        print('-dpsc2',fileprintS)
        end

        if mod_egaux
                if F_end~=0
                    if Polarisation==2 %numchoc>41200
                    F_end=min(F_end,18.75*4);
                    F_start=max(F_start,12.25*4);
                    F_pal_O=(0:nb_freq-1)*(F_end-F_start)/(nb_freq-1)+F_start;
                    F_pal_X=3/2*((0:nb_freq-1)*(F_end-F_start)/(nb_freq-1)+F_start);
                    elseif Polarisation==1
                    F_end=min(F_end,106/6*4);
                    F_pal=3/2*((0:nb_freq-1)*(F_end-F_start)/(nb_freq-1)+F_start);
                    elseif Polarisation==0
        			F_pal=(0:nb_freq-1)*(F_end-F_start)/(nb_freq-1)+F_start;
                    end
                elseif F_end==0
                    if Polarisation==0
                    AGAZ_nl=tsmat(numchoc,'AGAZ;Nl;VALEURS');
    				maxnlprog=max(AGAZ_nl(:,2))*1e19;
    				F_max=floor(1.05*sqrt(0.90*.95*maxnlprog*qe^2/me/epsi0)/2/pi/1e9);
    				F_max=max(F_max,60);F_max=min(F_max,75);
    				F_pal=(0:nb_freq-1)*(F_max-F_start)/(nb_freq-1)+F_start;
                    elseif Polarisation==1
                    F_max=106/6*4;
                    F_pal=3/2*((0:nb_freq-1)*(F_max-F_start)/(nb_freq-1)+F_start);
                    elseif Polarisation==2 %numchoc>41200
    				AGAZ_nl=tsmat(numchoc,'AGAZ;Nl;VALEURS');
    				maxnlprog=max(AGAZ_nl(:,2))*1e19;
    				F_max=floor(1.05*sqrt(0.90*.95*maxnlprog*qe^2/me/epsi0)/2/pi/1e9);
    				F_max=max(F_max,60);F_max=min(F_max,18.75*4);
                    F_pal_O=(0:nb_freq-1)*(F_end-F_start)/(nb_freq-1)+F_start;
                    F_pal_X=3/2*((0:nb_freq-1)*(F_end-F_start)/(nb_freq-1)+F_start);
                    end
                end
        elseif mod_inegaux
            %nb_freq=length(F_pal_O);
                if Polarisation==2 %numchoc>41200
                F_pal_O=F_pal(1:nb_freq);;
                F_pal_X=3/2*F_pal(1:nb_freq);
                end
        end
    end %elseif numchoc>40051

mdB=31;%10*log10(max(mssignal));
Mdb=60;%10*log10(max(diff(v1-mean(v1)).^2));
    
    if mdB<30,ymax=40;
    elseif mdB>40,ymax=65;
    else ymax=55;
    end
    
    if taille_seg<100000
	nfft=1024;
    elseif taille_seg<=200000
	nfft=4096*2;
    elseif	taille_seg>200000
	nfft=4096;
    end
    
fmax=2000;fmax_X=3000;
f=(-nfft/2:nfft/2-1)/nfft/dtDifdop;
kfC=find(f<100e3 &(f> -100e3));

%if (iter_nb==1) iter_nb=2; end 
    for iter=1:(nb_iter*length(TDIFDOP))
    % for iter=1:nb_iter
  		
        
         %wf=gdp_wavelet_reduce(sfreq,1200,5,16,2);
          % wf=wf2(1:length(sfreq));
          
            %Angle=mean(wf(tlent>=t0FA(iter) & tlent<t0FA(iter)+.10));
    t0FA(iter)=Ftime(kDecl(iter))-ignitron;
    wf=smooth(sfreq,128);
    Angle=mean(wf(tlent>=t0FA(iter) & tlent<t0FA(iter)+.10));
  	difdop_anglepol(iter)=Angle;%Angle;

  		if ~isempty(tnl);
   		nld=mean(nl(tnl>=t0FA(iter) & tnl<t0FA(iter)+.030));
  		else
   		nld=[];
        end
        
  	ne1d=mean(GNE(tgne>=t0FA(iter) & tgne<t0FA(iter)+.030),1);
  	ne2d=mean(GNE(tgne>=t0FA(iter) & tgne<t0FA(iter)+.030),2);
  	fp=sqrt(ne1d*qe^2/me/epsi0)/2/pi/1e9;
  	if ~isempty(ne1d);nec(iter)=ne1d;else,nec(iter)=NaN;end

  	m=2;n=ceil(nb_freq/2);
  	Sk=[];Sk2=[];
  	figure(10+iter);clf;
    title(sprintf('TS# %g, d. %g, t=%5.1fs, \\theta= %5.1f, ',...
                    numchoc,iter,t0FA(iter),Angle),'fontsize',14);
    		
    [difdop_phase,difdop_amplitude,difdop_angletor]=difdop_get_cor(numchoc);
  		
        for ifreq=1:nb_freq
        indseg=((1:taille_seg)+(ifreq-1)*taille_seg)+(iter-1)*nb_freq*taille_seg;
            
            % correction desequilibre et dephasage des voies
            %av1=mean(v1(indseg).^2);av2=mean(v2(indseg).^2);
            %if numchoc<40051
               % cor=exp(-i*7/180*pi)*sqrt(av1/av2);%.8; %*10^{-3};
           % else
              %  cor=exp(-i*3/180*pi)*sqrt(av1/av2);
            %end
            %x+i*y -> x*exp(i*phi)+ i*a*y
         window_pwelch=hanning(nfft);
           if max(indseg)<=length(v1)
           %[P,F]=spectrum(v1(indseg)+i*cor*v2(indseg),nfft,[],[],1/dtDifdop);
           %[P2,F2]=spectrum(v1(indseg).*exp(i*difdop_phase(1)/180*pi)+i*difdop_amplitude(1)*v2(indseg),nfft,[],[],1/dtDifdop);
           %[P,F]=pwelch(v1(indseg).*exp(i*difdop_phase(1)/180*pi)+i*difdop_amplitude(1)*v2(indseg),nfft,[],[],1/dtDifdop);
           [P,F]=pwelch(v1(indseg).*exp(i*difdop_phase(1)/180*pi)+i*difdop_amplitude(1)*v2(indseg),window_pwelch,[],nfft,1/dtDifdop);
           else
           indmax=fix((length(v1)-min(indseg)+1)/nfft);
           indsegt=min(indseg)-1+(1:indmax*nfft);
           %[P,F]=spectrum(v1(indsegt)+i*cor*v2(indsegt),nfft,[],[],1/dtDifdop);
           %[P,F]=pwelch(v1(indsegt).*exp(i*difdop_phase(1)/180*pi)+i*difdop_amplitude(1)*v2(indsegt),nfft,[],[],1/dtDifdop);
           [P,F]=pwelch(v1(indsegt).*exp(i*difdop_phase(1)/180*pi)+i*difdop_amplitude(1)*v2(indsegt),window_pwelch,[],nfft,1/dtDifdop);
           end
            
         Sk=[Sk fftshift(P(:,1))]; IC(iter,ifreq)=sum(Sk(kfC,ifreq))/nfft;
         %Sk2=[Sk2 fftshift(P2(:,1))];
         subplot(n,m,ifreq);hold on;
         plot(f'/1e3,10*log10(Sk(:,ifreq)),f(kfC)'/1e3,10*log10(Sk(kfC,ifreq)),':r')
        %set(gca,'xlim',[-fmax fmax],'ylim',[0 40]);
         set(gca,'xlim',[-fmax fmax]);
         grid on
        	% 	set(gca,'ylim',[-5 ymax]);grid on
         xlabel(' F [kHz]');ylabel('dB');
                
			if numchoc>41200
        		%	text(.65*fmax,.8*ymax,sprintf('F=%4.1f',F_pal_O(ifreq)))
            title(sprintf('F=%4.1f',F_pal_O(ifreq)))
			else
        		%text(.65*fmax,.8*ymax,sprintf('F=%4.1f',F_pal(ifreq)))
            title(sprintf('F=%4.1f',F_pal(ifreq)))
            end
        		
     
        end %for ifreq=1:nb_freq
        
    Skw{iter}=Sk;Skw2{iter}=Sk2;
    subplot(n,m,1);
        if numchoc>41200
        title(sprintf('TS# %g, d. %g, t=%5.1fs, \\theta= %5.1f, F = %5.1fGHz',...
                    numchoc,iter,t0FA(iter),Angle,F_pal_O(1)),'fontsize',14);
        else
        title(sprintf('TS# %g, d. %g, t=%5.1fs, \\theta= %5.1f, F = ',...
                    numchoc,iter,t0FA(iter),Angle,F_pal(1)),'fontsize',14);
        end
        
            
    		%subplot(n,m,2);
    		%title(sprintf(' nl=%5.1f, nec=%5.1f, fp=%5.1f GHZ',...
           %         nld,ne1d/1e19,fp),'fontsize',14);

  

 		if numchoc>41200
  		Sk_X=[];
  		figure(100+iter);clf;
        title(sprintf('TS# %g, d. %g, t=%5.1fs, \\theta= %5.1f?, ',...
                    numchoc,iter,t0FA(iter),Angle),'fontsize',14);
    	
            for ifreq=1:nb_freq
     		indseg=((1:taille_seg_W)+(ifreq-1)*taille_seg_W)+(iter-1)*nb_freq*taille_seg_W;
          			% correction desequilibre et dephasage des voies
          			%  av3=mean(v3(indseg).^2);av4=mean(v4(indseg).^2);
           			% cor=exp(i*5/180*pi)*sqrt(av3/av4);
                 if max(indseg)<=length(v3)
                 %[P,Fi]=spectrum(v3(indseg)+i*cor*v4(indseg),nfft,[],[],1/dtDifdop_W);
                 [P,Fi]=pwelch(v3(indseg).*exp(i*difdop_phase(2)/180*pi)+i*difdop_amplitude(2)*v4(indseg),window_pwelch,[],nfft,1/dtDifdop_W);
                 else
                 indmax=fix((length(v3)-min(indseg)+1)/nfft);
                 indsegt=min(indseg)-1+(1:indmax*nfft);
                 %[P,F]=spectrum(v3(indsegt)+i*cor*v4(indsegt),nfft,[],[],1/dtDifdop_W);
                 [P,Fi]=pwelch(v3(indsegt).*exp(i*difdop_phase(2)/180*pi)+i*difdop_amplitude(2)*v4(indsegt),window_pwelch,[],nfft,1/dtDifdop_W);
                 end
                 
        	Sk_X=[Sk_X fftshift(P(:,1))]; IC(iter,ifreq)=sum(Sk_X(kfC,ifreq))/nfft;  
        	subplot(n,m,ifreq);hold on;
        	plot(f'/1e3,10*log10(Sk_X(:,ifreq)),f(kfC)'/1e3,10*log10(Sk_X(kfC,ifreq)),':r')
            %	set(gca,'xlim',[-fmax fmax],'ylim',[0 40]);
            %set(gca,'xlim',[-fmax_X fmax_X]);
            grid on
        	%set(gca,'ylim',[-5 ymax]);grid on
        	xlabel(' F [kHz]');ylabel('dB');
        	%text(.65*fmax,.6*ymax,sprintf('F=%4.1f',F_pal_X(ifreq)))
            title(sprintf('F=%4.1f',F_pal_X(ifreq)))
     
            end %for ifreq=1:nb_freq
            
        Skw_X{iter}=Sk_X;
  		subplot(n,m,1);
    	title(sprintf('TS# %g, d. %g, t=%5.1fs, \\theta= %5.1f, F =%5.1fGHz ',...
                    		numchoc,iter,t0FA(iter),Angle,F_pal_X(1)),'fontsize',14);
    			

    	orient tall
    	print('-dpsc2','-append',fileprintS)
        end %if numchoc>41200
  
  
    end %for iter=1:(nb_iter*length(TDIFDOP))
Data=['./data/'];
	if exist(Data) ~=7
	Data='/Home/LV219680/DIFDOP_LV/data/';
    end
    
filemat=[Data,'Sdifdop',int2str(numchoc)]
disp(['sauvegarde des spectres calcules dans ',filemat])
save(filemat,'numchoc','username','nfft','dtDifdop*','t0FA','t0Fp','taille_seg*','Polarisation','difdop_anglepol','nec','F_pal*','Skw*','nb_freq');
%save(filemat,'numchoc','username','nfft','dtDifdop','t0FA','NbCoups_Rab','Polarisation','difdop_anglepol','nec','F_pal','Skw');


%elseif Mode_acqui==2  % continu
 %t=Ftime(1)+(0:length(v1)-1)/4e6;
 %figure
 %plot(t(1:32768)-ignitron,v1(1:32768))
            
 %nfft=4096;Sk=[];
 %[P,F]=spectrum(v1+i*v2,nfft,[],[],4e6);
 %Sk=[Sk fftshift(P(:,1))];   

% fmax=2000;dtDifdop=.25e-6;
% f=(-nfft/2:nfft/2-1)/nfft/dtDifdop;
% kfC=find(f<100e3 &(f> -100e3));
 %       figure;
% plot(f'/1e3,10*log10(Sk))
%        set(gca,'xlim',[-fmax fmax],'ylim',[0 40]);
 %       set(gca,'xlim',[-fmax fmax]);
 %       set(gca,'ylim',[0 35]);grid on
 %       xlabel(' F [kHz]');ylabel('dB');
 %title(sprintf('TS# %g, t=%5.1fs, \\theta= %5.1f?, ',...
 %                   numchoc,t(1)-ignitron,Angle_Ini),'fontsize',14);

 %[B,F,T]=specgram(v1+i*v2,nfft,4e6);
 %S_plot=20*log10(abs(fftshift(B)));
 %figure
 %imagesc(T+t(1),f/1e3,S_plot);
 %axis xy
 %colormap(jet)

%end % mod_acqui
     
end % exist fast acq
