%	RZ2RHO	effectue le changement de coordonees (R,Z) -> Rho
%---------------------------------------------------------------------------------------------
%
% fonction rz2rho.m
% 
% Cette fonction sert a effectuer le changement de coordonees (R,Z) -> Rho
%
% FONCTIONS UTILISEES : 
%
%	rz2rhomex.mex*
%
% SYNTAXE : 
%
%	[rho,erreur]=rz2rho(r,z,a,r0,z0,d0,{piqd,{e1,ep1}});
%
%	 {} denote les parametres optionnels
%
% ENTREES :
%
% Les coordonnees a fournir :
%
%	r		->		coordonnee 'R' de chaque point 
%	z		->		coordonnee 'Z' de chaque point 
%
%	remarque : les deux coordonnees doivent avoir memes dimensions
%
% Les parametres de la geometrie du plasma :
%
%	a		-> 		petit rayon du plasma (SAMIN de TEMPETE)
%	r0		-> 		grand rayon du plasma (SRMAJ de TEMPETE)
%	z0		-> 		decalage vertical du plasma (SZPOS de TEMPETE)
%
% Le profil de decentrement de Schafranov :
%
%	1 - profil parabolique :
%
%			d(rho)=d0*(1-rho^2)
%
%			d0		-> 		decentrement de shafranov en rho=0 (SD0MAG de TEQUILA)
%			piqd	-> 		[]
%		
%			dans ce cas l'ellipticite n'est pas prise en compte !0
%
%	2 - profil calculer sur donnees les magnetiques seules :
%
%			d(rho)=d0*(1-rho^piqd)
%
%			d0		-> 		decentrement de shafranov en rho=0 (SD0MAG de TEQUILA)
%			piqd	-> 		piquage du profil de decentrement (SSMAG de TEQUILA)
%
%	2 - profil calculer sur donnees les magnetiques et les donnees de la polarimetrie :
%
%			d(rho)=d0*(1-rho^piqd)
%
%			d0		-> 		decentrement de shafranov en rho=0 (SD0POL de TEQUILA)
%			piqd	-> 		piquage du profil de decentrement (SSPOL de TEQUILA)
%
%
% L'ellipticite du plasma :
%
%	1 - pas prise en compte  :
%
%		e1		-> 		[]
%		ep1		->  	[]
%
%	2 - prise en compte 
%
%		e1		-> 		elipticite au bord (SELLIP de TEQUILA)
%		ep1		->  	derive de l'elipticite au bord (SELLIPIQ de TEQUILA)
%
% 		le profil d'ellipticite est donnee par :
%
%			e(rho)=ep1*(1-rho^2)+e1
%
% SORTIE :
%
% 	rho		-> 		coordonnee 'Rho' de chaque point 
%	erreur	->		erreur de la recherche de zero (utilise pour les tests)
%
% fonction ecrite par J-F Artaud, poste 46-78 
% version 2, derniere mise a jour le 08/06/95
%---------------------------------------------------------------------------------------------
%
function [rho,erreur]=rz2rho(r,z,a,r0,z0,d0,piqd,e1,ep1)
		%
		% test des arguments d'entree
		%
		if nargin < 6
			error('nombre d''arguments dentree insufisant');
		elseif nargin <7
			piqd=[];
			e1=[];
			ep1=[];
			mode2=1;
			elip=0;
		elseif nargin <9
			e1=[];
			ep1=[];
			elip=0;
			mode2=0;
		elseif isempty(piqd),
			elip=0;
			mode2=1;
		elseif isempty(e1)|isempty(ep1)
			elip=0;
			mode2=0;
		else
			mode2=0;
			elip=1;
		end
		%
		% taille des entrees
		%
		if (~all(size(a)==size(r0)))| ...
			(~all(size(a)==size(z0)))| ...
			(~all(size(a)==size(d0)))| ...
			((mode2==0)&(~all(size(a)==size(piqd))))| ...
			((elip==1)&(~all(size(a)==size(e1))))| ...
			((elip==1)&(~all(size(a)==size(ep1)))),
			
			error('Les donnees de la geometrie du plasma doivent toutes avoir memes dimensions !')
		end
		%
		% les donnees 
		%
		if (~all(size(r)==size(z)))
			error('Les donnees r et z doivent toutes avoir memes dimensions !')			
		end
		%
		% selon la dimension de rho
		%
		if all(size(r)==1)
			r=r*ones(size(a));
			z=z*ones(size(a));
		elseif size(r,1)==1
			r=ones(size(a,1),1)*r;
			z=ones(size(a,1),1)*z;
		elseif size(r,1)~=size(a,1),
			error('La dimension 1 de r (ou z) doit etre la meme que celles de a, r0 ou z0 !');
		end
		%
		% selon la dimension d'espace de a , ...
		%
		if size(a,2)==1,
			a=a*ones(1,size(r,2));
			r0=r0*ones(1,size(r,2));
			z0=z0*ones(1,size(r,2));
			d0=d0*ones(1,size(r,2));
			if ~isempty(piqd),
				piqd=piqd*ones(1,size(r,2));
			end
			if ~isempty(e1),
				e1=e1*ones(1,size(r,2));		
				ep1=ep1*ones(1,size(r,2));
			end
		elseif size(a,2)~=size(r,2),
			error('La dimension 2 de a (ou r0, z0, d0, ...) doit etre 1 ou celle de r !')
		end
		%
		% gestion des modes
		%
		if (isempty(e1))|(elip==0)|(mode2==1),
			ep1=zeros(size(a));
			e1=ones(size(a));		
		end
		%
		if (isempty(piqd))|(mode2==1),
			piqd=2*ones(size(a));		
		end					
		%
		% debut
		%
		rho=nan*ones(size(r));
		erreur=rho;
		test=~(isnan(r)|isnan(z));
		if any(any(test)),
			[rho(test),erreur(test)]=rz2rhomex(r(test),z(test),a(test),r0(test),z0(test), ...
									d0(test),piqd(test),e1(test),ep1(test));
		end
		%
		% remplace rho <0 par nan
		%
		matnan=	(rho<0);
		cpt=nan*ones(sum(sum(matnan)),1);
		rho(matnan)=cpt;
		%
		% fin 
		%
end