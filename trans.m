function trans=trans(taucloud, Recloud, albsurf, solarzenith); 
% computing Arctic low level cloud transmittance using the algorithms of
% Fitzpatrick et al. (2004)
%

a1=0.58;
b1=0.74;
b2=-0.1612;
b3=-0.8343;
k1=1.9785;
k2=0.2828;
k3=2.3042;
aa=a1+(1-a1)*exp(-k1*taucloud);
bb = b1*(1+b2.*exp(-k2*taucloud)+b3*(exp(-k3*taucloud))); 

trans=(aa+bb*cosd(solarzenith)) ./ (1+(cc-dd.*albsurf).*taucloud) ;
end