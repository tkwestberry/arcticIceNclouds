% compute surface SW transmittance
%
% The observational variables are
%
%   CALIPSO cloud layer products, e.g., depol (volume depolarization ratios
%   of cloud layer), bk532 and bk1064 (layer integrated backscatter by clouds at 532
%   nm and 1064 nm), atm532 and atm1064 (column integrated atmospheric backscatter above clouds)
%   CALIPSO cloud backscatter profiles, e.g., profile532 and profile 1064
%   (cloud backscatter profile at 532 nm and 1064 nm),
%   profiledep(depolarization profiles of cloud backscatter)
%   
%   Cloud optical depth (taucloud) and droplet size (Recloud) derived from CALIPSO using the the machine learning
%   trained with colloated MODIS optical depth and droplet size (Hu and Lu,
%   2024)
% 
%   Cloud microphysical property, such as extinction coefficient (cloudext) and droplet number concentration (cloudNd)
%   derived from CALIPSO cloud backscatter (Hu et al., 2021)
%
%   The CALIPSO cloud product is publically available on google drive (we
%   will send the link upon request)
%  
%
% load the CALIPSO cloud database
for kyear=2008:2017
    for kmonth=1:12
        kkk=kyear*100+kmonth;
        fnmonth=num2str(kkk);
load(['calipso_cloud_profile_' fnmonth '.mat']);

% computing cloud optical depth and effective radius using the Hu and Lu (2024) algorithm
load taunets.mat;
NNoutput = taunets(depol,bk532,bk1064,atm532,atm1064);
taucloud(kyear-2007,kmonth,:) = NNoutput(1,:);
Recloud(kyear-2007,kmonth,:) = NNoutput(2,:);

% computing cloud extinction coeff and droplet number concentration using Hu et al. (2021) 
[cldext cldNd] = microphysics(profile532,profiledp,profile1064,Recloud); 

% estimating surface albedo using the monthly snow/ice map
albsurf = albedo(timecalipso,lat,lon);
% estimating solar zenith angle 
sza = solarzenith(timecalipso,lat,lon);

% computing surface net downward solar transmittance using Fitzpatrick et
% al. (2004)
transmit(kyear,kmonth,:) = trans(taucloud, Recloud, albsurf, sza); 

    end
end

save(['surfnet_' fnmonth '.mat'],'lat','lon','taucloud','Recloud','cldext','cldNd');

