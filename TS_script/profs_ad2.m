
close all
clear

workdir=('/Home/LV219680/DIFDOP_LV/')




numc=input(' shot number : ');
ad=1
addpath /Home/LV219680/prog_TS/equilibre/
if ~exist('username'),username=getenv('LOGNAME');end;
if ~exist('userdir'),userdir=getenv('HOME');end;
FigureDir=['/Home/LV219680/DIFDOP_LV/DopExp/profiles'];

%load plasma parameters

duree_plasma=tsbase(numc, 'RDUREE');
B0=tsbase(numc, 'RBTOR');
[Ip, tps_Ip] = tsbase(numc, 'SIMAG');	
[zeff,tps_zeff]=tsbase(numc,'szfbrm'); 
[nl,tps_nl]=tsbase(numc,'gnl%4'); 
[R0,tps_R0] = tsbase(numc, 'SRMAJ');
[a, tps_a] = tsbase(numc, 'SAEQA');
[am, tps_am] = tsbase(numc, 'SAMIN');
[Z0, tps_Z0] = tsbase(numc, 'SZPOS');
[d0, tps_Z0] = tsbase(numc, 'SD0MAG');
[piqd, tps_Z0] = tsbase(numc, 'SSMAG');


[q_efit,tps_efit] = tsbase(numc, 'GEFQ');
[rho_efit] = tsbase(numc, 'SEFPSIN');
%load power







[Pohm,tps_ohm]=tsbase(numc,'spohm');
%[Picrh, tps_icrh] = tsbase(numc, 'GBILAN%3');
[Picrh, tps_icrh] = tsbase(numc, 'GPUIFCI%4');
if isempty(Picrh)==1
  Picrh = zeros(length(tps_icrh),1);
end
	
%reference time (can be easly changed)
if ~isempty(tps_Ip)
tref=tps_Ip;
 else
   disp('No Ip data')
return
end
l=length(tref);
disp ('reference time => tps_Ip')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(zeff)
  zeff=interp1q(tps_zeff,zeff,tref);
else
  disp('No data for Zeff') 
   zeff=zeros(1,l);
end


Ip=interp1q(tps_Ip,Ip,tref);

if ~isempty(nl)
  nl=interp1q(tps_nl,nl,tref);
 else 
disp('No data for nl') 
 nl=zeros(1,l);
end

if ~isempty(Pohm)
Pohm=interp1q(tps_ohm,Pohm,tref);
else 
disp('No data for Pohm') 
 Pohm=zeros(1,l);
end

if ~isempty(Picrh)
Picrh=interp1q(tps_icrh,Picrh,tref);
else 
disp('No data for Picrh') 
  Picrh=zeros(l,1);
end	
   

%%%%%% DIFDOP data %%%%%%%%%%%%%%%%%%%%%%%

[svideo,tlent,srms,tlent,sposverin,tlent,sfreq,tlent] = tsmat(numc,...
     'DIFDOP-SRMSVIDEO','DIFDOP-SRMSCOS','DIFDOP-SPOSVERIN','DIFDOP-SFREQ');

% en volts
sposverin=sposverin*10/2048;
srms=srms*10/2048;
svideo=svideo*.014-8.3;  % x 10 dBm calibration AT
% le signal sfreq est depuis 2006 devolu a l'inclinometre
sfreq=asin(sfreq/2048)*180/pi-7.5;% inclinom?tre asin(SV/10)
%%%%  calibration signal verin en angle
if numc > 39100 %campagne 2007 (calibr CH inclinom?tre)
  sangle=sposverin.^2*0.00785+sposverin*1.09 + 4.64; % 2007
  if numc > 40050, sangle=sangle+7;end;
elseif numc > 34300 %campagne 2005 (mesure 30/05/2005 0.01966 1.18296 7.42788)
  sangle=sposverin.^2*0.0115+sposverin*1.1156 +6.7398;
elseif numc > 32040 %campagne 2003, novembre
  sangle=(sposverin+2.3427)/.3365 -7.5;
else
  sangle=[];%pas de calib ou de signal moteur
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
h1=figure(1);
coord=get(h1,'Position');
set(h1,'Position',[650,500,coord(3),coord(4)])
set(gca,'FontSize',16)
title(['# ', num2str(numc) 'B(T) = ', num2str(B0)]) 
%plot(tref, Ip,'k',tref,smooth(Pohm,10),'b',tref,Picrh, 'r',tref,Picrh+smooth(Pohm,10),'m',tref,nl./1e19,'g',tref,zeff,'c', 'LineWidth',2)
plot(tref, Ip,'k',tref,smooth(Pohm,10),'b',tref,nl./1e19,'g',tref,zeff,'c', 'LineWidth',2)
hold on
plot(tref,Picrh, 'r',tref,Picrh+smooth(Pohm,10))
axis([0 duree_plasma+2 0 1.2*max(nl./1e19)])
grid on
legend('I_p','P_{OH}','nl','zeff','P_{ICRH}','P_{tot}')
xlabel('time is s')
title(['# ', num2str(numc),'B(T) = ', num2str(B0)])


h2=figure(2);
coord=get(h2,'Position');
set(h2,'Position',[650,200,coord(3),coord(4)./2])
set(gca,'FontSize',16)
plot(tlent,svideo,tlent,srms,tlent,sfreq,tlent,sangle)
grid on
axis([0 duree_plasma+2 1.2*min(svideo) 1.2*max(sangle)])
title(['# ', num2str(numc) 'DIFDOP acquisition']) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('load ECE data')
%ECE
[te_ece, tps_ece, poubelle, poubelle, r_ece] = tsbase(numc,'GSHTE', 'GSHR');
if isempty (te_ece)
        [te_ece, tps_ece, poubelle, poubelle, r_ece] = tsbase(numc,'GSHTENV', 'GSHR');
end
 

disp('load Thomson data')
%Thomson
       [te_thom, tps_thom, z_thom,cthom] = tsbase(numc,'GTETHOM'); 
        [ne_thom, tps_thom, z_thom,cthom] = tsbase(numc,'GNETHOM'); 
        %tps_thom=[];
%tps_thom=[];
% if ~isempty(tps_thom)
%       [tps_thom,rhothom,tethom,rayonlaser,zspec,certif_thom,certif_mag,mode_sortie_th]=tsrhothom(numc,1)
        
% else
 %  disp('No data from Thomson')
%end

disp('load Interfero data') 
%Interfero
[gne, tps_gne, rho_gne] = tsbase(numc, 'GNE');



disp('load DREFLEC profiles')
%dreflec
[ne_dreflec,tps_dreflec] = tsbase(numc, 'GREFNEX');
[r_dreflec,tps_dreflec] = tsbase(numc, 'GREFRX');
if isempty(tps_dreflec)
  disp('no profiles from DREFLEC in the database')
%check=input('Do you want to check in the local directory /reflec/TREFLEX_profils/ ? (1 or 0) :   ')
  check=0;
if check==1
filemat=['/home/sccp/gttm/reflec/TREFLEX_profils/choc' int2str(numc) '_prof.mat'];
load(filemat)
  r_dreflec=RX;
ne_dreflec=NEX;
tps_dreflec=tps;
end
end 

disp('load DREFLUC profiles')
%drefluc
[ne_drefluc_ref,tps_drefluc_ref,ii] = tsbase(numc, 'GNEFLUCREF');
[r_drefluc_ref,tps_drefluc_ref,ii] = tsbase(numc, 'GRFLUCREF');

[ne_drefluc_gne,tps_drefluc_gne,ii] = tsbase(numc, 'GNEFLUCGNE');
[r_drefluc_gne,tps_drefluc_gne,ii] = tsbase(numc, 'GRFLUCGNE');


if isempty(tps_drefluc_ref)&isempty(tps_drefluc_gne)
  disp('no profiles from DREFLUC in the database')
  %check=input('Do you want to check in the local directory /drefluc/Outils_drefluc/ ? (1 or 0) :   ')
  check=0;
  if check==1
       addpath /home/sccp/gttm/drefluc/Outils_drefluc/
       NE=NE_drefluc(numc);


       if ~isempty(NE.t_gne)
             r_drefluc_gne=NE.R_gne;
             ne_drefluc_gne=NE.ne_gne;
             tps_drefluc_gne=NE.t_gne;
       end

       if ~isempty(NE.t_ref)
            ne_refluc_ref=NE.R_ref;
            r_drefluc_ref=NE.ne_ref;
            tps_drefluc_ref=NE.t_ref;
      end

 end 

end

	%Choice of the time period
	t1=input(' Time period [t1,t2](s) : t1 = ');
	t2=input(' and t2 =  ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mean between t1 and t2

i1=iround(tps_a,t1);
i2=iround(tps_a,t2);

a_mean=mean(a(i1:i2));
Z0_mean=mean(Z0(i1:i2));
R0_mean=mean(R0(i1:i2));
d0_mean=mean(d0(i1:i2));
piqd_mean=mean(piqd(i1:i2));
zeff_mean=mean(zeff(i1:i2));

ie1=iround(tps_efit,t1);
ie2=iround(tps_efit,t2);

q_mean=mean(q_efit(ie1:ie2,:));





%%%%%%% ECE %%%%%%%%%%%%
cor=input('correction of some channels ? : ');

if ~isempty(te_ece)
	ind_ece_1=iround(tps_ece,t1);
	ind_ece_2=iround(tps_ece,t2);
%if (cor==1)
%    te_ece_cor(ind_ece_1:ind_ece_2,:)=te_ece(te_ece(ind_ece_1:ind_ece_2,:)~=0);
%end

if (cor==1)
    if (numc==43076 | numc==43085 | numc==43082 | numc==43083 | numc==43080 | numc==43087)
te_ece_cor(ind_ece_1:ind_ece_2,:)=[te_ece(ind_ece_1:ind_ece_2,1:27) te_ece(ind_ece_1:ind_ece_2,29)];%43085
r_ece_cor(ind_ece_1:ind_ece_2,:)=[r_ece(ind_ece_1:ind_ece_2,1:27) r_ece(ind_ece_1:ind_ece_2,29) ]; 
    end
%te_ece_cor(ind_ece_1:ind_ece_2,:)=[te_ece(ind_ece_1:ind_ece_2,17:32)];
%r_ece_cor(ind_ece_1:ind_ece_2,:)=[r_ece(ind_ece_1:ind_ece_2,17:32)]; 
if (numc>45475 & numc<45513)
te_ece_cor(ind_ece_1:ind_ece_2,:)=[te_ece(ind_ece_1:ind_ece_2,1:29) ];%43089
r_ece_cor(ind_ece_1:ind_ece_2,:)=[r_ece(ind_ece_1:ind_ece_2,1:29)];
end

if (numc==43089 | numc==43090 | numc==43100 | numc==43093 | numc==43096| numc==43097)
te_ece_cor(ind_ece_1:ind_ece_2,:)=[te_ece(ind_ece_1:ind_ece_2,1:27) te_ece(ind_ece_1:ind_ece_2,29) te_ece(ind_ece_1:ind_ece_2,31:32) ];%43089
r_ece_cor(ind_ece_1:ind_ece_2,:)=[r_ece(ind_ece_1:ind_ece_2,1:27) r_ece(ind_ece_1:ind_ece_2,29) r_ece(ind_ece_1:ind_ece_2,31:32)];
end


if (numc==43081 |numc==43072 | numc==43094| numc==43069 )
te_ece_cor(ind_ece_1:ind_ece_2,:)=[te_ece(ind_ece_1:ind_ece_2,1:27) te_ece(ind_ece_1:ind_ece_2,29)];%43081
r_ece_cor(ind_ece_1:ind_ece_2,:)=[r_ece(ind_ece_1:ind_ece_2,1:27) r_ece(ind_ece_1:ind_ece_2,29)];
end




if (numc==43379)
  te_ece_cor(ind_ece_1:ind_ece_2,:)=[te_ece(ind_ece_1:ind_ece_2,1:5) te_ece(ind_ece_1:ind_ece_2,7:32)];%43379
r_ece_cor(ind_ece_1:ind_ece_2,:)=[r_ece(ind_ece_1:ind_ece_2,1:5) r_ece(ind_ece_1:ind_ece_2,7:32)];
end




%te_ece_cor(ind_ece_1:ind_ece_2,:)=[te_ece(ind_ece_1:ind_ece_2,8:20)];
%r_ece_cor(ind_ece_1:ind_ece_2,:)=[r_ece(ind_ece_1:ind_ece_2,8:20)];

te_ece_mean=mean(te_ece_cor(ind_ece_1:ind_ece_2,:));
r_ece_mean=mean(r_ece_cor(ind_ece_1:ind_ece_2,:));
else
        te_ece_mean=mean(te_ece(ind_ece_1:ind_ece_2,:));
        r_ece_mean=mean(r_ece(ind_ece_1:ind_ece_2,:));
end

else
        disp('no profiles from ECE')
        te_ece_mean = zeros(1,100);
        r_ece_mean = zeros(1,100);
end

z_ece=zeros(1,length(r_ece_mean));
%if (B0<3)
 % rep_ece=find(r_ece_mean<R0_mean);
%else
 %rep_ece=find(r_ece_mean>R0_mean);
%end

%te_ece_red=te_ece(:,rep_ece);
%te_ece_mean=te_ece_mean(rep_ece);
%z_ece=z_ece(rep_ece);
%r_ece_mean=r_ece_mean(rep_ece);
r_ece_mean=r_ece_mean;

%equilibrum from ??
[rho_ece,erreur]=rz2rho(r_ece_mean,z_ece,a_mean,R0_mean,Z0_mean,d0_mean,piqd_mean);

%equilibrium from EFIT
%tse=tsefitgeo(numc,linspace(t1,t2,200));
tse=[];
if ~isempty(tse)
  [rho_ece_efit,theta_ece_efit]=rztorho_refluc(r_ece_mean,z_ece,tse);
end

%equilibrium from TMAG
tmag= tmageo(numc,linspace(t1,t2,200));

if ~isempty(tmag)
  [rho_ece_tmag,theta_ece_tmag]=rztorho_refluc(r_ece_mean,z_ece,tmag);
end


%check_ece1=input('Do you want to check all Te profiles from ECE ? (1 or 0) :' );
check_ece1=0;
if check_ece1==1
	  h3=figure(3);
          coord=get(h3,'Position');
name=get(h3,'Name');
set(h3,'Name','Te profiles from ECE')
          set(h3,'Position',[650,500,coord(3),coord(4)])
          set(gca,'FontSize',16)
	    if ~isempty(rho_ece)
	      for i=ind_ece_1:ind_ece_2
	             plot(rho_ece,te_ece(i,:))
		     pause
	      end
	    end
end

%check_ece2=input('Do you want to check mean of Te profiles from ECE ? (1 or 0) : ' );
check_ece2=1;
if check_ece2==1	  
		h4=figure(4);
          coord=get(h4,'Position');
          set(h4,'Position',[650,500,coord(3),coord(4)])
          set(gca,'FontSize',16)
	 if ~isempty(rho_ece)
	   
	            
%	   plot(rho_ece,te_ece_red(ind_ece_1:ind_ece_2,:),'k')
	   hold on
	 plot(rho_ece,te_ece_mean,'r','LineWidth',2)
	      end
title(['# ', num2str(numc) 'TE Profiles fron ECE']) 
end
%%%%%%% THOMSON %%%%%%%%%%%%
%if ~isempty(tps_thom)
%    ind_thom_1=iround(tps_thom,t1);
	%ind_thom_2=iround(tps_thom,t2);
  %  for i=ind_thom_1:ind_thom_2
  %      [rho_thom_s(i,:),ind_s]=sort(rhothom(i,:));
   %     te_thom_s(i,:)=te_thom(i,ind_s);
  %  end
   
  %  rho_thom_mean=mean(rho_thom_s(ind_thom_1:ind_thom_2,:));
  %  te_thom_mean=mean(te_thom_s(ind_thom_1:ind_thom_2,:));
%end
 


%%%%%%% TI CXRS %%%%%%%%%%%%
[Ti,tps_cxs,RTi] = tsbase(numc, 'GTICXS');
%[RTi] = tsbase(numc, 'SRAYCXS');
[Wphi,tps_cxs,RWphi] = tsbase(numc, 'GWPHICXS');

if ~isempty(tps_cxs)
    ZTi=zeros(length(Ti),1);
   % tmag= tmageo(numc,linspace(tps_cxs-0.1,tps_cxs+0.1,10));

%tmag= tmageo(numc,tps_cxs);
tmag= tmageo(numc,7);

        if ~isempty(tmag)
            [rho_cxs_tmag,theta_cxs_tmag]=rztorho_refluc(RTi,ZTi',tmag);
            rho_cxs_tmag_m=mean(rho_cxs_tmag);
        end
end


%%%%%%%% INTERFEREO %%%%%%%%%%%%%%%%
if ~isempty(tps_gne)
    ind_gne_1=iround(tps_gne,t1);
	ind_gne_2=iround(tps_gne,t2);
    gne_mean=mean(gne(ind_gne_1:ind_gne_2,:));
end
%%%%%%%% DREFLEC %%%%%%%%%%%%%%%%

if ~isempty(tps_dreflec)	
    ind_dreflec_1=iround(tps_dreflec,t1);
	ind_dreflec_2=iround(tps_dreflec,t2);
%[r_dreflec_mean,ne_dreflec_mean]=mean_dreflec(r_dreflec,ne_dreflec,ind_dreflec_1,ind_dreflec_2);
if ind_dreflec_1==ind_dreflec_2
    r_dreflec_mean=r_dreflec(ind_dreflec_1,:); 
    ne_dreflec_mean=ne_dreflec(ind_dreflec_1,:);
    
else
r_dreflec_mean=mean(r_dreflec(ind_dreflec_1:ind_dreflec_2,:));
ne_dreflec_mean=mean(ne_dreflec(ind_dreflec_1:ind_dreflec_2,:));
end

z_dreflec=zeros(1,length(r_dreflec_mean));
z_dreflec_t=zeros(ind_dreflec_2,length(r_dreflec_mean));
rho_test=zeros(ind_dreflec_2,length(r_dreflec_mean));
theta_test=zeros(ind_dreflec_2,length(r_dreflec_mean));
%equilibrum from ??
[rho_dreflec,erreur]=rz2rho(r_dreflec_mean,z_dreflec,a_mean,R0_mean,Z0_mean,d0_mean,piqd_mean);

%equilibrium from EFIT
%tse=tsefitgeo(numc,linspace(t1,t2,200));
%if ~isempty(tse)
%  [rho_dreflec_efit,theta_dreflec_efit]=rztorho_refluc(r_dreflec_mean,z_dreflec,tse);
%end

%equilibrium from TMAG
tmag= tmageo(numc,tps_dreflec(ind_dreflec_1:ind_dreflec_2));

if ~isempty(tmag)
  [rho_dreflec_tmag,theta_dreflec_tmag]=rztorho_refluc(r_dreflec_mean,z_dreflec,tmag);
end


else
disp('no profiles from DREFLEC')
ne_dreflec_mean = zeros(1,100);
r_dreflec_mean = zeros(1,100);
end	   


%check_dreflec1=input('Do you want to check all Ne profiles from DREFLEC ? (1 or 0) :' );
check_dreflec1=0;
if check_dreflec1==1
	  h3=figure(3);
          coord=get(h3,'Position');
          set(h3,'Position',[650,500,coord(3),coord(4)])
          set(gca,'FontSize',16)
	    if ~isempty(r_dreflec)
	      for i=ind_dreflec_1:ind_dreflec_2
	             %i
                 kk=i-ind_dreflec_1+1;
                 tmagi.shot=tmag.shot;
                 tmagi.time=tmag.time(kk);
                 tmagi.piqE=tmag.piqE(kk);
                 tmagi.ellip=tmag.ellip(kk);
                 tmagi.rmag=tmag.rmag(kk);
                 tmagi.rc=tmag.rc(kk);
                 tmagi.zc=tmag.zc(kk);
                 tmagi.d=tmag.d(kk);
                 tmagi.piqD=tmag.piqD(kk);
                 tmagi.rminor=tmag.rminor(kk);
                 
	             plot(r_dreflec(i,:),ne_dreflec(i,:))
                % [rho_test(i,:),theta_test(i,:)]=rztorho_refluc(r_dreflec(i,:),z_dreflec_t(i,:),tmagi);
                 
		     %pause
	      end
	    end
end

	    %check_dreflec2=input('Do you want to check mean of Ne profiles from DREFLEC ? (1 or 0) :' ); 
check_dreflec2=1;
if check_dreflec2==1
h5=figure(5);
          coord=get(h5,'Position');
          set(h5,'Position',[650,500,coord(3),coord(4)])
          set(gca,'FontSize',16)
	 if ~isempty(r_dreflec)
	   
	            
	   plot(r_dreflec(ind_dreflec_1:ind_dreflec_2,:)',ne_dreflec(ind_dreflec_1:ind_dreflec_2,:)','k')
	   hold on
	 plot(r_dreflec_mean,ne_dreflec_mean,'r','LineWidth',2)
	      end
end
title(['# ', num2str(numc) '  NE Profiles fron DREFLEC']) 


%%%%%%%% DREFLUC %%%%%%%%%%%%%%%%

if ~isempty(tps_drefluc_ref)	
        ind_drefluc_ref1=iround(tps_drefluc_ref,t1);
	ind_drefluc_ref2=iround(tps_drefluc_ref,t2);
if ind_drefluc_ref1==ind_drefluc_ref2
   r_drefluc_ref_mean= r_drefluc_ref(ind_drefluc_ref1,:);
   ne_drefluc_ref_mean=ne_drefluc_ref(ind_drefluc_ref1,:);
else
r_drefluc_ref_mean=mean(r_drefluc_ref(ind_drefluc_ref1:ind_drefluc_ref2,:));
ne_drefluc_ref_mean=mean(ne_drefluc_ref(ind_drefluc_ref1:ind_drefluc_ref2,:));
%r_drefluc_ref_mean=r_drefluc_ref(ind_drefluc_ref1:ind_drefluc_ref2,:);
%ne_drefluc_ref_mean=ne_drefluc_ref(ind_drefluc_ref1:ind_drefluc_ref2,:);
end


z_drefluc_ref=zeros(1,length(r_drefluc_ref_mean));


%equilibrum from ??
%if (B0>3.2) 
 % rep=(1:1:length(r_drefluc_ref_mean));
%end

%if (B0<3.2) 
  rep=find(r_drefluc_ref_mean>R0_mean);
%end

[rho_drefluc_ref,erreur]=rz2rho(r_drefluc_ref_mean(rep),z_drefluc_ref(rep),a_mean,R0_mean,Z0_mean,d0_mean,piqd_mean);

%equilibrium from EFIT
%tse=tsefitgeo(numc,linspace(t1,t2,100));
%if ~isempty(tse)
%  [rho_drefluc_ref_efit,theta_drefluc_ref_efit]=rztorho_refluc(r_drefluc_ref_mean(rep),z_drefluc_ref(rep),tse);
%end

%equilibrium from TMAG
tmag= tmageo(numc,linspace(t1,t2,200));

if ~isempty(tmag)
  [rho_drefluc_ref_tmag,theta_drefluc_ref_tmag]=rztorho_refluc(r_drefluc_ref_mean(rep),z_drefluc_ref(rep),tmag);
end


else
disp('no profiles from DREFLUC WITH REF')
ne_drefluc_ref_mean = zeros(1,100);
r_drefluc_ref_mean = zeros(1,100);


end	   


if ~isempty(tps_drefluc_gne)	
        ind_drefluc_gne1=iround(tps_drefluc_gne,t1);
	ind_drefluc_gne2=iround(tps_drefluc_gne,t2);
r_drefluc_gne_mean=mean(r_drefluc_gne(ind_drefluc_1:ind_drefluc_2,:));
ne_drefluc_gne_mean=mean(ne_drefluc_gne(ind_drefluc_1:ind_drefluc_2,:));
z_drefluc_gne=zeros(1,length(r_drefluc_gne_mean));

%equilibrum from ??

[rho_drefluc_gne,erreur]=rz2rho(r_drefluc_gne_mean,z_drefluc_gne,a_mean,R0_mean,Z0_mean,d0_mean,piqd_mean);

%equilibrium from TMAG
tmag= tmageo(numc,linspace(t1,t2,200));
if ~isempty(tmag)
  [rho_drefluc_gne_tmag,theta_drefluc_gne_tmag]=rztorho_refluc(r_drefluc_gne_mean,z_drefluc_gne,tmag);
end

else
disp('no profiles from DREFLUC WITH GNE')
ne_drefluc_gne_mean = zeros(1,100);
r_drefluc_gne_mean = zeros(1,100);


end	   
	    %check_drefluc=input('Do you want to check mean of Ne profiles from DREFLUC ? (1 or 0) :' );
check_drefluc=0; 

%if check_drefluc==1

h6=figure(6);
          coord=get(h6,'Position');
          set(h6,'Position',[650,500,coord(3),coord(4)])
          set(gca,'FontSize',16)
	    if ~isempty(r_drefluc_ref)	            
	      plot(r_drefluc_ref(ind_drefluc_ref1:ind_drefluc_ref2,:)',ne_drefluc_ref(ind_drefluc_ref1:ind_drefluc_ref2,:)','k')
	      hold on
	      plot(r_drefluc_ref_mean,ne_drefluc_ref_mean,'r','LineWidth',2)
	    end

title(['# ', num2str(numc) 'NE Profiles fron DREFLUC REF']) 
%end

 

%	 fig=input('Do you want to save figures ? (1 or 0) :  ');
fig=1

if fig==1  

  
filegraph1=[FigureDir,'/',int2str(numc),'PROF','-',int2str(t1),'-',int2str(t2),'-all.ps'];

figure(1)  
%print('-dpsc',filegraph1)
figure(2)  
%print('-dpsc',filegraph1,'-append')
figure(4)  
%print('-dpsc',filegraph1,'-append')
figure(5)  
%print('-dpsc',filegraph1,'-append')
figure(6)  
%print('-dpsc',filegraph1,'-append')


end




  disp('Chosen equilibrium --> from TMAG')
 


rho=linspace(0,1.2,100);
rho_ece_tmag_mean=mean(rho_ece_tmag);
rma=max(rho_ece_tmag_mean);
rmi=min(rho_ece_tmag_mean);
rep_hf=find(gradient(rho_ece_tmag_mean)<0);
rep_t=find(gradient(rho_ece_tmag_mean)>0);
if isempty(rep_t)
        Te=interp1(rho_ece_tmag_mean,te_ece_mean,rho,'cubic').*1e3;
     prof.Te=smooth(Te',5)';
else
    Te=interp1(rho_ece_tmag_mean(rep_t),te_ece_mean(rep_t),rho,'cubic').*1e3;
   % prof.Te=smooth(Te',5)';
   prof.Te=Te;
   
end
rr=find(rho<rmi | rho>rma);

prof.Te(rr)=0;
%if ~isempty(tps_thom)
%rmat=max(rhothom_mean);
%rmit=min(rhothom_mean);

%rhothom_mean=mean(rhothom(ind_thom_1:ind_thom_2,1:6));
%tethom_mean=mean(te_thom(ind_thom_1:ind_thom_2,1:6));
%prof.Te2=interp1(rhothom_mean,tethom_mean,rho,'cubic').*1e3;
%prof.Te3=interp1(rho_thom_mean,te_thom_mean,rho,'cubic').*1e3;
%end
if ~isempty(tmag)
    if ~isempty(tps_dreflec)
        %prof.Ne=interp1(rho_dreflec_tmag,ne_dreflec_mean,rho);
        prof.Ne=interp1(mean(rho_dreflec_tmag),ne_dreflec_mean,rho);
        rmae=max(mean(rho_dreflec_tmag));
        rmie=min(mean(rho_dreflec_tmag));
        re=find(rho<rmie | rho>rmae);
        prof.Ne(re)=0;
    else
        prof.Ne=zeros(1,length(rho));
    end
    if ~isempty(tps_drefluc_ref)
        prof.Ne2=interp1(mean(rho_drefluc_ref_tmag),ne_drefluc_ref_mean(rep),rho);
        rmau=max(mean(rho_drefluc_ref_tmag));
        rmiu=min(mean(rho_drefluc_ref_tmag));
        ru=find(rho<rmiu | rho>rmau);
        prof.Ne2(ru)=0;
    else
       prof.Ne2=zeros(1,length(rho));
    end
end
prof.Ne3=interp1(rho_gne,gne_mean,rho);
rmai=max(rho_gne);
rmii=min(rho_gne);
ri=find(rho<rmii | rho>rmai);
prof.Ne3(ri)=0;
prof.rho=rho;

if ~isempty(tps_drefluc_ref)

%if (B0>3.2) 
%prof.Ne2=interp1(mean(rho_drefluc_ref_tmag),ne_drefluc_ref_mean,rho);
%prof.Ne2=interp1(rho_drefluc_ref_tmag,ne_drefluc_ref_mean,rho);
%end
%if (B0<3.2) 
prof.Ne2=interp1(mean(rho_drefluc_ref_tmag),ne_drefluc_ref_mean(rep),rho);
%end
 %else 
%prof.Ne2=0.*prof.Ne;
end
prof.q=interp1(rho_efit,q_mean,rho);
%tex
%prof.B0
%prof.Ti=prof.Te;
if ~isempty(tps_cxs)
prof.Ti=interp1(rho_cxs_tmag,Ti(2,:),rho);
%prof.Ti=interp1(rho_cxs_tmag,Ti,rho);
else
   prof.Ti=prof.Te; 
end
%R-->R0_mean
rr=0.8
%rr=input('valeurs locales, quel rayon (r/a) ? : ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  calculation of density gradient length  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rr1=rr-0.05;
rr2=rr+0.05;
ind=iround(prof.rho,rr);
ind1=iround(prof.rho,rr1);
ind2=iround(prof.rho,rr2);
[p,s]=polyfit(prof.rho(ind1:ind2),prof.Ne(ind1:ind2),1);
[p2,s2]=polyfit(prof.rho(ind1:ind2),prof.Ne2(ind1:ind2),1);
[p3,s3]=polyfit(prof.rho(ind1:ind2),prof.Ne3(ind1:ind2),1);


 
LN=abs(prof.Ne./(gradient(prof.Ne)./gradient(prof.rho)));
LNN=abs(prof.Ne./(gradient(prof.Ne,prof.rho)));
RLN=(a_mean+R0_mean)./LN;
RLNN=abs(gradient(prof.Ne,prof.rho.*(a_mean+R0_mean))./prof.Ne);

LN2=abs(prof.Ne2./(gradient(prof.Ne2)./gradient(prof.rho)));
RLN2=(a_mean+R0_mean)./LN2;

LN3=abs(prof.Ne3./(gradient(prof.Ne3)./gradient(prof.rho)));
RLN3=(a_mean+R0_mean)./LN3;


RLNI=abs((a_mean+R0_mean)*p(1)./mean(prof.Ne(ind1:ind2)));
RLNI2=abs((a_mean+R0_mean)*p2(1)./mean(prof.Ne2(ind1:ind2)));
RLNI3=abs((a_mean+R0_mean)*p3(1)./mean(prof.Ne3(ind1:ind2)));


figure
plot(prof.rho,prof.Ne,'b')
hold on
plot(prof.rho,prof.Ne2,'r')
plot(prof.rho,prof.Ne3,'k')
%plot([rr1 rr1],[0 1e20],'--')
%plot([rr2 rr2],[0 1e20],'--')
xlabel('r/a')
ylabel('Ne')
legend(sprintf('DREFLEC : p = %2.4g ',p(1)),sprintf('DREFLUC : p = %2.4g',p2(1)),sprintf('INTERF : p = %2.4g',p3(1)))
grid on

figure
plot(prof.rho,RLN,'b')
hold on
plot(prof.rho,RLN2,'r')
plot(prof.rho,RLN3,'k')
plot([rr1 rr1],[0 15],'--')
plot([rr2 rr2],[0 15],'--')
xlabel('r/a')
ylabel('RLN')
legend(sprintf('DREFLEC : p = %2.4g ',RLNI),sprintf('DREFLUC : p = %2.4g',RLNI2),sprintf('INTERF : p = %2.4g',RLNI3))
axis([0 1.4 0 15])
grid on







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  calculation of electron temperature gradient length  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LTE=abs(prof.Te./(gradient(prof.Te)./gradient(prof.rho)));
RLTE=(a_mean+R0_mean)./LTE;
[p,s]=polyfit(prof.rho(ind1:ind2),prof.Te(ind1:ind2),1);
RLTEI=abs((a_mean+R0_mean)*p(1)./mean(prof.Te(ind1:ind2)));



figure
plot(prof.rho,prof.Te,'b')
hold on
plot([rr1 rr1],[0 max(prof.Te)],'--')
plot([rr2 rr2],[0 max(prof.Te)],'--')
xlabel('r/a')
ylabel('Te (eV)')
legend(sprintf('ECE : p = %2.4g ',p(1)))
grid on



figure
plot(prof.rho,RLTE,'b')
hold on
plot([rr1 rr1],[0 15],'--')
plot([rr2 rr2],[0 15],'--')
xlabel('r/a')
ylabel('RLTE')
legend(sprintf('ECE : p = %2.4g ',RLTEI))
axis([0 1.4 0 15])
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  calculation of ion temperature gradient length  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(tps_cxs)
LTI=abs(prof.Ti./(gradient(prof.Ti)./gradient(prof.rho)));
RLTI=(a_mean+R0_mean)./LTI;
[p,s]=polyfit(prof.rho(ind1:ind2),prof.Ti(ind1:ind2),1);
RLTII=abs((a_mean+R0_mean)*p(1)./mean(prof.Ti(ind1:ind2)));



figure
plot(prof.rho,prof.Ti,'b')
hold on
plot([rr1 rr1],[0 max(prof.Ti)],'--')
plot([rr2 rr2],[0 max(prof.Ti)],'--')
xlabel('r/a')
ylabel('Ti (eV)')
legend(sprintf('CXRS : p = %2.4g ',p(1)))
grid on


figure
plot(prof.rho,RLTI,'b')
hold on
plot([rr1 rr1],[0 15],'--')
plot([rr2 rr2],[0 15],'--')
xlabel('r/a')
ylabel('RLTi')
legend(sprintf('CXRS : p = %2.4g ',RLTII))
axis([0 1.4 0 15])
grid on

end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  calculation of beta, nustar et rhostar  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%constantes
c=3e8;
qe=1.602e-19;
ps0=8.8542e-12;
depi=2*pi;
me=9.1096e-31;
mp=1.672e-27;
mi=2*mp;
muo=4*pi*1e-7;
eps0=8.8542e-12;

%Blocal=interp1(rhoB,Bvide,rho(ifr));

wce=qe*B0./me; wci=qe*B0./(2*mp); 
epsilon=a_mean./R0_mean;
    
    if isfield(prof,'Te')==1 
        vteprof=((2*qe*prof.Te)/me).^(.5);
        vtiprof=((2*qe*prof.Ti)/mi).^(.5);
        csprof=((2*qe*prof.Te)./mi).^(.5);
        rhoe=vteprof./wce;
        prof.rhostare=rhoe./a_mean;
   




%lnc=15.2-0.5.*log(nec./1e20)+log(Tec);% CRONOS
%tei=3*eps0*eps0*depi^1.5*sqrt(me)/qe^4*((qe*Tec.*1e3).^(3/2))./lnc./nec%./zeffc;% CRONOS
%nustar=Rmaj*ones(size(xc)).*qx./vte./tei./(xc.*amin./(Rmaj*ones(size(xc)))).^(1.5);% CRONOS

        if isfield(prof,'Ne')==1
            prof.beta=2*muo*(prof.Ne.*prof.Te.*qe)./(B0.^2).*200;
            lncprof=15.2-0.5.*log(prof.Ne./1e20)+log(prof.Te./1e3); 
            lncprofi=17.3-0.5.*log(prof.Ne./1e20)+3/2*log(prof.Ti./1e3);
            lncprofei=17.3-0.5.*log(prof.Ne./1e20)+3/2*log(prof.Te./1e3);
            teiprof=3*eps0*eps0*depi^1.5*sqrt(me)/qe^4*((qe*prof.Te).^(3/2))./lncprof./prof.Ne;
            
            tiiprof=12*eps0*eps0*pi^1.5*sqrt(2*mp)/qe^4*((qe*prof.Ti).^(3/2))./lncprofi./prof.Ne;
            
            tiieprof=12*eps0*eps0*pi^1.5*sqrt(2*mp)/qe^4*((qe*prof.Te).^(3/2))./lncprofei./prof.Ne; %compatible avec papier XG
            
            prof.nustar=R0_mean.*prof.q./vteprof./teiprof./epsilon.^(1.5)*sqrt(2);
            
            prof.nustari=R0_mean.*prof.q./vtiprof./tiiprof./epsilon.^(1.5)*sqrt(2);
            
            %prof.nustari_loc=R0_mean.^(3/2).*prof.q./vtiprof./tiiprof./epsilon.*sqrt(2./(a_mean.*rho));
            prof.nustari_loc=R0_mean.*prof.q./vtiprof./tiiprof./(a_mean.*rho./R0_mean).^(1.5); %%% le bon à utiliser, celui passer à GD le 09/02/2017 !!il n'y a pas de sqrt(2) car dans wb c'est la vitesse perp
            prof.nustar_loc=R0_mean.*prof.q./csprof./tiieprof./(a_mean.*rho./R0_mean).^(1.5)
            %prof.nustarie_loc=R0_mean.^(3/2).*prof.q./csprof./tiieprof./epsilon.*sqrt(2./(a_mean.*rho));
            
            prof.colli=2.3031*1e-5*lncprofi.*(prof.Ne./1e19).*a_mean./((prof.Ti./1000).^2); %coll for GENE and Lref = a 
            prof.coll=2.3031*1e-5*lncprofei.*(prof.Ne./1e19).*a_mean./((prof.Te./1000).^2);
            prof.coll_XG=4*sqrt(pi)./3./(4*pi*eps0).^2*qe.^4./mi^2.*lncprofei.*(prof.Ne)./vtiprof.^3; %! dans cette definition Vt=sqrt(T/m) donc prendre 1/tiiprof
            
            wb0=(a_mean./R0_mean).^(1.5)*sqrt(qe.*prof.Ti/mi)./prof.q./R0_mean;
            Nbob=18; delta_rip=0.01
            qr1=1./(Nbob*prof.q).^(1.5);
            qr2=Nbob.*prof.q*(delta_rip./epsilon).^2;
            
            ZZ=linspace(Z0_mean-a_mean,Z0_mean+a_mean,100);
            RR=linspace(R0_mean,R0_mean+a_mean,100);
            pphi=linspace(0,2*pi,200);
            theta=linspace(0,2*pi,200);
           
            
            %%CALCUL DU RIPPLE
            b1=R0_mean-2.04;
            b2=0.26;
            a1=2.2e-4;
            a2=5;
            a3=1.6;
            rho1=linspace(0,1,100);
            r=rho1.*a_mean;
  
            for ir=1:100
                for it=1:200
                    E=1/2+b2*(r(ir)*cos(theta(it))+b1);
                    F=r(ir).^2+2*b1*r(ir)*cos(theta(it))+b1.^2;
                    rprim_bis=sqrt((1/b2.^2)*(E-sqrt(E.^2-b2.^2*F)));
                    delta_des_bis(ir,it)=a1*exp(a2*rprim_bis+a3*rprim_bis.^2);
            
                    rprim=sqrt((0.5+r(ir).*b2.*cos(theta(it))+b1.*b2-((0.5+r(ir).*b2.*cos(theta(it))+b1.*b2).^2-b2.^2.*(r(ir).^2+2*r(ir).*b1.*cos(theta(it))+b1.^2)).^(1/2)))./b2;
                    delta_des(ir,it)=a1*exp(a2*rprim+a3*rprim.^2);
                       
                    
                    R(ir,it)=R0_mean+rho1(ir).*a_mean*cos(theta(it));
                    Z(ir,it)=Z0_mean+rho1(ir).*a_mean*sin(theta(it));
            
                end
            end
            %plot(rho1,delta_des(:,1))
            %contour(R,Z,delta_des,100)
            %contour(R,Z,delta_des,50,'ShowText','on')


            
            qq=interp1(rho,prof.q,rho1);

            qr2_prof=Nbob.*qq'.*(delta_des(:,1)./epsilon).^2;  
            
            
            %%%%  test regime ripple-plateau %%%%%
            figure
            plot(rho1,qr2_prof) %%si nustar<<qr2_prof => regime ripple-plateau collisional
            
            hold on
            plot(rho,prof.nustari_loc,'r')
            
           
           
            %for iR=1:100
               %for iZ=1:100
                    %for ip=1:200
                        %[bt(ip),br(ip),bz(ip)] = bmex(RR(iR),ZZ(iZ),pphi(ip),924);
                        %delta(iR,iZ)=(max(bt)-min(bt))./(max(bt)+min(bt));
                   % end
               % end
           % end
            
            %for iR=1:100
               % for iZ=1:100
                 %   delta(iR,iZ)=(max(bt(iR,iZ,:))-min(bt(iR,iZ,:)))./(max(bt(iR,iZ,:))+min(bt(iR,iZ,:)));
                %end
           % end
            
        
        if isfield(prof,'Ne2')==1 
            prof.beta2=2*muo*(prof.Ne2.*prof.Te.*qe)./(B0.^2).*200;
            lncprof2=15.2-0.5.*log(prof.Ne2./1e20)+log(prof.Te./1e3); 
            lncprofi2=17.3-0.5.*log(prof.Ne2./1e20)+3/2*log(prof.Ti./1e3);
            lncprofei2=17.3-0.5.*log(prof.Ne2./1e20)+3/2*log(prof.Te./1e3);
            teiprof2=3*eps0*eps0*depi^1.5*sqrt(me)/qe^4*((qe*prof.Te).^(3/2))./lncprof2./prof.Ne2;
            tiiprof2=12*eps0*eps0*pi^1.5*sqrt(2*mp)/qe^4*((qe*prof.Ti).^(3/2))./lncprofi2./prof.Ne2;
            tiieprof2=12*eps0*eps0*pi^1.5*sqrt(2*mp)/qe^4*((qe*prof.Te).^(3/2))./lncprofei2./prof.Ne2;
            
            prof.nustar2=R0_mean.*prof.q./vteprof./teiprof2./epsilon.^(1.5)*sqrt(2);
           
            prof.nustari2=R0_mean.*prof.q./vtiprof./tiiprof2./epsilon.^(1.5)*sqrt(2);
            
            prof.nustari2_loc=R0_mean.^(3/2).*prof.q./vtiprof./tiiprof2./epsilon.*sqrt(2./(a_mean.*rho));
            
            prof.nustarie2_loc=R0_mean.^(3/2).*prof.q./csprof./tiieprof2./epsilon.*sqrt(2./(a_mean.*rho));
             
            prof.colli2=2.3031*1e-5*lncprofi2.*(prof.Ne2./1e19).*a_mean./((prof.Ti./1000).^2); %coll for GENE and Lref = a 
            prof.coll2=2.3031*1e-5*lncprofei2.*(prof.Ne2./1e19).*a_mean./((prof.Te./1000).^2);
            
%prof.nustar2i=R0_mean.*prof.q./vtiprof./teiprof2./epsilon.^(1.5)*sqrt(2);
        end

        if isfield(prof,'Ne3')==1
           prof.beta3=2*muo*(prof.Ne3.*prof.Te.*qe)./(B0.^2).*200;
            lncprof3=15.2-0.5.*log(prof.Ne3/1e20)+log(prof.Te./1e3);
            teiprof3=3*eps0*eps0*depi^1.5*sqrt(me)/qe^4*((qe*prof.Te).^(3/2))./lncprof3./prof.Ne3;
            lncprofi3=17.3-0.5.*log(prof.Ne3/1e20)+3/2*log(prof.Te./1e3);
            tiiprof3=12*eps0*eps0*pi^1.5*sqrt(2*mp)/qe^4*((qe*prof.Te).^(3/2))./lncprofi3./prof.Ne3;
            
            prof.nustar3=R0_mean.*prof.q./vteprof./teiprof3./epsilon.^(1.5)*sqrt(2);
            prof.nustari3=R0_mean.*prof.q./vtiprof./tiiprof3./epsilon.^(1.5)*sqrt(2);
            
        end
    end
    
    if isfield(prof,'Te2')==1
        vteprof2=((2*qe*prof.Te2)/me).^(.5);
        rhoe2=vteprof2./wce;
        prof.rhostare2=rhoe2./a_mean;
        
        if isfield(prof,'Ne3')==1 
            prof.beta4=2*muo*(prof.Ne3.*prof.Te2.*qe)./(B0.^2).*200;
            lncprof4=15.2-0.5.*log(prof.Ne3/1e20)+log(prof.Te2./1e3);
            teiprof4=3*eps0*eps0*depi^1.5*sqrt(me)/qe^4*((qe*prof.Te2).^(3/2))./lncprof4./prof.Ne3;
            prof.nustar4=R0_mean.*prof.q./vteprof2./teiprof4./epsilon.^(1.5)*sqrt(2);
        end
    end
    
	
    if isfield(prof,'Ti')==1
        vteprof3=((2*qe*prof.Ti)/(2*mp)).^(.5);
        rhoe3=vteprof3./wci;
        prof.rhostari=rhoe3./a_mean;
     prof.beta5=2*muo*(prof.Ne3.*(prof.Te+prof.Ti)*qe)./(B0.^2).*100;
    end
	
	

    
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('%%%%%% Valeurs locales  %%%%%')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    prof.local.r=rr;
    
    prof.local.nus=prof.nustar(ind); 
    prof.local.nus2=prof.nustar2(ind);  
    prof.local.nus3=prof.nustar3(ind);
    
    prof.local.rhose=prof.rhostare(ind);
    if ~isempty(tps_cxs)
    prof.local.rhosi=prof.rhostari(ind); 
    end
    
    prof.local.nus_m=mean(prof.nustar(ind1:ind2)); 
    prof.local.nus2_m=mean(prof.nustar2(ind1:ind2)) ; 
    prof.local.nus3_m=mean(prof.nustar3(ind1:ind2));
    
    prof.local.rhose_m=mean(prof.rhostare(ind1:ind2));
    
    if ~isempty(tps_cxs)
    prof.local.rhosi_m=mean(prof.rhostari(ind1:ind2));
    end
    
    disp(sprintf('nus2 = %2.4g',prof.local.nus2))
    disp(sprintf('nus2_m = %2.4g',prof.local.nus2_m))
    disp(sprintf('rhose = %2.4g',prof.local.rhose))
    disp(sprintf('rhose_m = %2.4g',prof.local.rhose_m))
  
    if ~isempty(tps_cxs)
    disp(sprintf('rhosi = %2.4g',prof.local.rhosi))
    disp(sprintf('rhosi_m = %2.4g',prof.local.rhosi_m))
    end
    
	% wesson page 69 and 661
	%lnc=15.2-0.5.*log(nex./1e20)+log(tex./1e3);
	
        %lncprof2=15.2-0.5.*log(prof.Ne2/1e20)+log(prof.Te);
	%tei=3*eps0*eps0*depi^1.5*sqrt(me)/qe^4*((qe*tex).^(3/2))./lnc./nex;
	%tei_wess=1.09.*10e16.*(tex./1000).^1.5./nex./lnc;
	
        %teiprof2=3*eps0*eps0*depi^1.5*sqrt(me)/qe^4*((qe*prof.Te).^(3/2))./lncprof2./prof.Ne2;




%nustar2=R0_mean.*prof.q./vteprof./teiprof2./epsilon.^(1.5)*sqrt(2);


	%nustarprof=Rprof.*qprof./vteprof./teiprof./epsilonprof.^(1.5)/sqrt(2);
	%nustarprof(:,1)=nustarprof(:,2);




h7=figure(15);
coord=get(h7,'Position');
name=get(h7,'Name');
set(h7,'Name','Dimensionless profiles ')
          set(h7,'Position',[650,500,coord(3),coord(4)])
          set(gca,'FontSize',16)
  subplot(221)
  plot(rho,prof.beta)
hold on
%plot(rho,betaprof2)
  grid on
%  axis([0 1 0 1.2*max(prof.beta)])
  xlabel('\rho')
  ylabel('\beta_e (%)')
  subplot(222)
  plot(rho,prof.nustar)
hold on
  %plot(rho,nustar2)
  grid on
  xlabel('\rho')
  ylabel('\nu*')
  axis([0 1 0 2])
  subplot(223)  
  plot(rho,prof.rhostare) 
  grid on
%  axis([0 1 0 1.2*rhostare(1)])
  xlabel('\rho')
  ylabel('\rho*_e')
title('Dimensionless profiles ')

if fig==1  

 
figure(7)  
print('-dpsc',filegraph1,'-append')


end
end

%file=input('Do you want to save files ? (1 or 0) :  ');
file=1
savemat=[FigureDir,'/',int2str(numc),'PROF','-',int2str(t1),'-',int2str(t2),'.mat'];
%save(savemat,'r_dreflec_mean', 'ne_dreflec_mean','rho_ece','te_ece_mean','r_drefluc_ref','ne_drefluc_ref','rho','rhostare','nustar','betaprof','r_drefluc_ref_mean', 'ne_drefluc_ref_mean','r_dreflec', 'ne_dreflec','prof')
figure(7)  
print('-dpsc',filegraph1,'-append')
profq=interp1(rho_efit,q_mean,prof.rho);
RR=interp1(rho_ece_tmag_mean,r_ece_mean,rho);

%f_gam0=(csprof./RR)./(2*pi);
w_gam0=(csprof./RR);
%f_gam1=(sqrt(2+1./(profq.^2)).*csprof./RR)./(2*pi);
%f_gam2=(sqrt(2+1./(profq.^2)+7/2).*csprof./RR)./(2*pi); %%% formule Yanick (ref Watari) avec Te=Ti
%f_gam3=(sqrt(2*prof.Te./prof.Ti+1./(profq.^2)+7/2).*vtiprof./RR)./(2*pi); %avec ti
w_gam=sqrt(11)./2.*vtiprof./RR.*(1+(66)./(121*prof.q.^2)).^(1/2);
taue=0.8
w_gam_taue2=((7+4*taue)./4).*vtiprof.^2./RR.^2.*(1+2*(23+16*taue+4*taue.^2)./((7+4*taue).^2*prof.q.^2));

w_gam_gao2=(7/4+taue).*vtiprof.^2./RR.^2;
gamma_gao=-3./(14+8*taue).*RR;%.*prof.q./tiiprof./vtiprof;

gamma_landau=sqrt(2*pi).*prof.q.^2.*(vtiprof./(prof.q.*RR)).*exp(-prof.q.^2.*RR.^2.*w_gam0.^2./vtiprof.^2).*(prof.q.^4.*RR.^4.*w_gam0.^4./vtiprof.^4);

save(savemat)
r=rho*a_mean;
Bprof=ones(1,length(r))*R0_mean.*(ones(1,length(r)).*B0)./(ones(1,length(r)).*R0_mean+r);

vidia=prof.Ti./Bprof.*gradient(prof.Ne)./gradient(r)./prof.Ne;
vidia2=prof.Ti./Bprof.*gradient(prof.Ne2)./gradient(r)./prof.Ne2;

vedia=prof.Te./Bprof.*gradient(prof.Ne)./gradient(r)./prof.Ne;
vedia2=prof.Te./Bprof.*gradient(prof.Ne2)./gradient(r)./prof.Ne2;

return
%if ripple prediction %load /Home/LV219680/CXRS-DIFDOP/CXS_47181
%tmag= tmageo(numc,tps_cxs(1));
%[rho_cxs_tmag,theta_cxs_tmag]=rztorho_refluc(pos',zeros(length(Ticxs(:,1)),1)',tmag);
%prof.Ti=interp1(rho_cxs_tmag,Ticxs(:,1),rho);

ATi=1./prof.Ti.*gradient(prof.Ti)./gradient(rho*a_mean);
ATe=1./prof.Te.*gradient(prof.Te)./gradient(rho*a_mean);
Ane=1./prof.Ne.*gradient(prof.Ne)./gradient(rho*a_mean);


Vripple_rp=prof.Ti./Bprof.*(Ane+1.5*ATi);
Vripple_c=prof.Ti./Bprof.*(Ane+3.37*ATi);

Vripple_rp_e=prof.Te./Btot_exp'.*(Ane+1.5*ATe);
Vripple_c_e=prof.Te./Btot_exp'.*(Ane+3.37*ATe);

ATicxs=1./Ticxs(:,1).*gradient(Ticxs(:,1))./gradient(rho_cxs_tmag*a_mean);
gne=1./prof.Ne.*gradient(prof.Ne)./gradient(rho*a_mean);
%prof.necxs=interp1(rho,prof.Ne,rho_cxs_tmag);
Anecxs=interp1(rho,gne,rho_cxs_tmag);
%Anecxs=1./prof.necxs.*gradient(prof.necxs)./gradient(rho_cxs_tmag*a_mean);
B_cxs=interp1(rho,Bprof,rho_cxs_tmag);
Vripple_rp_cxs=Ticxs(:,1)./B_cxs'.*(Anecxs+1.5*ATicxs);
