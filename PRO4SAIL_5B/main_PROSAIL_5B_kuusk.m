% %*************************************************************************
% %*                                                                       *
% 	main_PROSAIL
%
% 09 01 2011
% This program allows modeling reflectance data from canopy
% - modeling leaf optical properties with PROSPECT-5 (feret et al. 2008)
% - modeling leaf inclination distribution function with the subroutine campbell
% (Ellipsoidal distribution function caracterised by the average leaf
% inclination angle in degree), or dladgen (2 parameters LIDF)
% - modeling canopy reflectance with 4SAIL (Verhoef et al., 2007)

% This version has been implemented by Jean-Baptiste Feret
% Jean-Baptiste Feret takes the entire responsibility for this version
% All comments, changes or questions should be sent to:
% jbferet@stanford.edu

% References:
% 	Verhoef et al. (2007) Unified Optical-Thermal Four-Stream Radiative
% 	Transfer Theory for Homogeneous Vegetation Canopies, IEEE TRANSACTIONS
% 	ON GEOSCIENCE AND REMOTE SENSING, VOL. 45, NO. 6, JUNE 2007
% 	F閞et et al. (2008), PROSPECT-4 and 5: Advances in the Leaf Optical
% 	Properties Model Separating Photosynthetic Pigments, REMOTE SENSING OF
% 	ENVIRONMENT
% The specific absorption coefficient corresponding to brown pigment is
% provided by Frederic Baret (EMMAH, INRA Avignon, baret@avignon.inra.fr)
% and used with his autorization.
% the model PRO4SAIL is based on a version provided by
%	Wout Verhoef
%	NLR
%	April/May 2003

% The original 2-parameter LIDF model is developed by and described in:
% 	W. Verhoef, 1998, "Theory of radiative transfer models applied in
%	optical remote sensing of vegetation canopies", Wageningen Agricultural
%	University,	The Netherlands, 310 pp. (Ph. D. thesis)
% the Ellipsoidal LIDF is taken from:
%   Campbell (1990), Derivtion of an angle density function for canopies
%   with ellipsoidal leaf angle distribution, Agricultural and Forest
%   Meteorology, 49 173-176
%*                                                                       *
%*************************************************************************
clc

TypeLidf=1;
% if 2-parameters LIDF: TypeLidf=1
if (TypeLidf==1)
    % LIDFa LIDF parameter a, which controls the average leaf slope
    % LIDFb LIDF parameter b, which controls the distribution's bimodality
    %	LIDF type 		a 		 b
    %	Planophile 		1		 0
    %	Erectophile    -1	 	 0
    %	Plagiophile 	0		-1
    %	Extremophile 	0		 1
    %	Spherical 	   -0.35 	-0.15
    %	Uniform 0 0
    % 	requirement: |LIDFa| + |LIDFb| < 1
    LIDFa	=	-0.35;
    LIDFb	=	-0.15;
    % if ellipsoidal LIDF: TypeLidf=2
elseif (TypeLidf==2)
    % 	LIDFa	= average leaf angle (degrees) 0 = planophile	/	90 = erectophile
    % 	LIDFb = 0
    LIDFa	=	30;
    LIDFb	=	0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LEAF CHEM & STR PROPERTIES	%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cab		=	40;		% chlorophyll content (礸.cm-2)
Car		=	8;		% carotenoid content (礸.cm-2)
Cbrown	=	0.0;	% brown pigment content (arbitrary units)
Cw		=	0.01;	% EWT (cm)
Cm		=	0.009;	% LMA (g.cm-2)
N		=	1.5;	% structure coefficient

data=dataSpec_P5B;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Soil Reflectance Properties	%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rsoil1 = dry soil
% rsoil2 = wet soil
Rsoil1=data(:,10);Rsoil2=data(:,11);
psoil	=	1;		% soil factor (psoil=0: wet soil / psoil=1: dry soil)
rsoil0=psoil*Rsoil1+(1-psoil)*Rsoil2;

rsoil0(648-400+1,1)=0.0741;
rsoil0(858-400+1,1)=0.3667;

%angles
va_red=zeros(722,4);
load Users_VA_Kuusk_SAILH_RED.dat
va_red(:,1:2)=Users_VA_Kuusk_SAILH_RED;
clear Users_VA_Kuusk_SAILH_RED

va_nir=zeros(815,4);
load Users_VA_Kuusk_SAILH_NIR.dat
va_nir(:,1:2)=Users_VA_Kuusk_SAILH_NIR;
clear Users_VA_Kuusk_SAILH_NIR

% red band
for t=1:722
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%	4SAIL canopy structure parm	%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LAI		=	1.86;     	% leaf area index (m^2/m^2)
    hspot	=	0.25;       % hot spot
    tts		=	52.;		% solar zenith angle (?
    tto		=	va_red(t,1);		% observer zenith angle (?
    
    psi		=	va_red(t,2);         % azimuth (?
    if psi>180
        psi=psi-360;
    end
    psi=abs(psi);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%          CALL PRO4SAIL       %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rdot: hemispherical-directional reflectance factor in viewing direction
    % rsot: bi-directional reflectance factor
    % rsdt: directional-hemispherical reflectance factor for solar incident flux
    % rddt: bi-hemispherical reflectance factor
    [rdot,rsot,rddt,rsdt]=PRO4SAIL_kuusk(N,Cab,Car,Cbrown,Cw,Cm,LIDFa,LIDFb,TypeLidf,LAI,hspot,tts,tto,psi,rsoil0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	direct / diffuse light	%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the direct and diffuse light are taken into account as proposed by:
    % Francois et al. (2002) Conversion of 400?100 nm vegetation albedo
    % measurements into total shortwave broadband albedo using a canopy
    % radiative transfer model, Agronomie
    % Es = direct
    % Ed = diffuse
    Es=data(:,8);Ed=data(:,9);
    rd=pi/180;
    % skyl	=	0.847- 1.61*sin((90-tts)*rd)+ 1.04*sin((90-tts)*rd)*sin((90-tts)*rd); % % diffuse radiation
    % PARdiro	=	(1-skyl)*Es;
    % PARdifo	=	(skyl)*Ed;
    %
    % % resv : directional reflectance
    % resv	= (rdot.*PARdifo+rsot.*PARdiro)./(PARdiro+PARdifo);
    
    
    
    skyl	=	0.180; % % diffuse radiation
    PARdiro	=	1-skyl;
    PARdifo	=	skyl;
    
    % resv : directional reflectance
    resv	= (rdot.*PARdifo+rsot.*PARdiro)./(PARdiro+PARdifo);
    
    
    
    va_red(t,3)=resv(648-400+1);
    va_red(t,4)=resv(858-400+1);
    
end

% NIR band
for t=1:815
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%	4SAIL canopy structure parm	%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LAI		=	1.86;     	% leaf area index (m^2/m^2)
    hspot	=	0.25/4.2;       % hot spot
    tts		=	45.;		% solar zenith angle (?
    tto		=	va_nir(t,1);		% observer zenith angle (?
    
    psi		=	va_nir(t,2);         % azimuth (?
    if psi>180
        psi=psi-360;
    end
    psi=abs(psi);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%          CALL PRO4SAIL       %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rdot: hemispherical-directional reflectance factor in viewing direction
    % rsot: bi-directional reflectance factor
    % rsdt: directional-hemispherical reflectance factor for solar incident flux
    % rddt: bi-hemispherical reflectance factor
    [rdot,rsot,rddt,rsdt]=PRO4SAIL_kuusk(N,Cab,Car,Cbrown,Cw,Cm,LIDFa,LIDFb,TypeLidf,LAI,hspot,tts,tto,psi,rsoil0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	direct / diffuse light	%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the direct and diffuse light are taken into account as proposed by:
    % Francois et al. (2002) Conversion of 400?100 nm vegetation albedo
    % measurements into total shortwave broadband albedo using a canopy
    % radiative transfer model, Agronomie
    % Es = direct
    % Ed = diffuse
    Es=data(:,8);Ed=data(:,9);
    rd=pi/180;
    % skyl	=	0.847- 1.61*sin((90-tts)*rd)+ 1.04*sin((90-tts)*rd)*sin((90-tts)*rd); % % diffuse radiation
    % PARdiro	=	(1-skyl)*Es;
    % PARdifo	=	(skyl)*Ed;
    %
    % % resv : directional reflectance
    % resv	= (rdot.*PARdifo+rsot.*PARdiro)./(PARdiro+PARdifo);
    
    
    
    skyl	=	0.081; % % diffuse radiation
    PARdiro	=	1-skyl;
    PARdifo	=	skyl;
    
    % resv : directional reflectance
    resv	= (rdot.*PARdifo+rsot.*PARdiro)./(PARdiro+PARdifo);
      
    va_nir(t,3)=resv(648-400+1);
    va_nir(t,4)=resv(858-400+1);
    
end

dlmwrite('E:\非均质混合像元\RAMI\KUUSK\pine\SAILNEW_RED.dat',va_red(:,3),'delimiter','\t','precision',3);
dlmwrite('E:\非均质混合像元\RAMI\KUUSK\pine\SAILNEW_NIR.dat',va_nir(:,4),'delimiter','\t','precision',3);

% x=zeros(29,1);
% for t=1:15
%     x(t,1)=-5*(15-t);
% end
% 
% for t=1:14
%     x(t+15,1)=5*t;
% end
% 
% subplot(2,2,1),plot(x,va(1:29,3),'r','Linewidth',2);
% subplot(2,2,2),plot(x,va(29+1:29+29,3),'k','Linewidth',2);
% subplot(2,2,3),plot(x,va(1:29,4),'r','Linewidth',2);
% subplot(2,2,4),plot(x,va(29+1:29+29,4),'k','Linewidth',2);
% dlmwrite('Refl_CAN.txt',[data(:,1),resv],'delimiter','\t','precision',5)
