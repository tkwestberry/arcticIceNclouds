clear

fnd= dir('Cloud_Day*mat');
fnn= dir('Cloud_Night*mat');

 for kyr=8:1:20
     k0=0;
            for kmn=1:12

for kdn=1:2
dn='Night';
if (kdn == 2)
    dn='Day';
end

km=num2str(100+kmn);
fseason=['Month' km(2:3)];
ky=num2str(100+kyr);
fyear=['Year-20' ky(2:3)];
for kkk=1:length(fnd)
fn=fnd(kkk).name;
fnmonth=fn(14:15);  
fnyear=fn(12:13);
if (dn(1) == 'N')
    fn=fnn(kkk).name;
    fnmonth=fn(16:17);  
    fnyear=fn(14:15);
end
%DJF
%fseason='DJF';
% fseason=['Month' fnmonth];
if ( fnmonth==km(2:3) & fnyear==ky(2:3))
%     fseason=['Month_' fnmonth];
load(fn);
km(2:3)
k0
if (k0==0) 
alti=CALIOP_alti;
temp=CALIOP_temp;
ext=CALIOP_ext;
lwc=CALIOP_lwc;
Nd=CALIOP_Nd;
re=CALIOP_re;
latitude=lat;
longitude=lon;
    dp=depo;
    bkg=bk532;
    bkr=bk1064;
    tmdn=UTC-fix(UTC);
    tm=UTC;
k0=1;
else
    alti=[alti CALIOP_alti];
    temp=[temp CALIOP_temp];
    ext=[ext CALIOP_ext];
    lwc=[lwc CALIOP_lwc];
    Nd=[Nd CALIOP_Nd];
    re=[re CALIOP_re]; 
    latitude=[latitude lat];
    longitude=[longitude lon];
    dp=[dp depo];
    bkg=[bkg bk532];
    bkr=[bkr bk1064];
    tmdn=[tmdn UTC-fix(UTC)];
    tm=[tm UTC];
end
end
end

end

end

highNd=movmean(Nd,2,'omitnan');

indarc=find(Nd>0);
nyr=kyr-7;

myr=1:5000000;

    aalt(myr,nyr)=NaN;
    aNd(myr,nyr)=NaN;
    alwc(myr,nyr)=NaN;
    aext(myr,nyr)=NaN;
    are(myr,nyr)=NaN;
    atmdn(myr,nyr)=NaN;
    atm(myr,nyr)=NaN;
    alon(myr,nyr)=NaN;
    alat(myr,nyr)=NaN;
    abkg(myr,nyr)=NaN;
    abkr(myr,nyr)=NaN;
    adp(myr,nyr)=NaN;
    zday(myr,kyr)=NaN;
    zzday(myr,kyr)=NaN;
%     nn1(myr,kyr)=NaN; 
%     nn2(myr,kyr)=NaN;
indwrong=find(dp<0); 
dp(1,indwrong)=0;
indwrong=find(Nd>500);
Nd(1,indwrong)=500;
indwrong=find(Nd<1);
Nd(1,indwrong)=1;


myr=1:1:length(indarc);
    aalt(myr,nyr)=alti(indarc);
    aNd(myr,nyr)=Nd(indarc);
    alwc(myr,nyr)=lwc(indarc);
    aext(myr,nyr)=ext(indarc);
    are(myr,nyr)=re(indarc);
    atmdn(myr,nyr)=tmdn(indarc);
    atm(myr,nyr)=tm(indarc);
    alon(myr,nyr)=longitude(indarc);
    alat(myr,nyr)=latitude(indarc);
    abkg(myr,nyr)=bkg(indarc);
    abkr(myr,nyr)=bkr(indarc);
    adp(myr,nyr)=dp(indarc);
    mnyr(nyr)=length(indarc);
    

 end


figure;
plot(atm(:,1),'.b');

 for kyr=1:13
 % convert atm to #day 
 xyr=fix( atm(:,kyr) /10000 );
 xmth=fix( (atm(:,kyr)-xyr*10000) / 100 );
 axmth(:,kyr)=xmth;
 xday=atm(:,kyr)-xyr*10000-xmth*100 ;
 zzday(:,kyr)=xyr+(xmth-1)/12+xday/366;
 ind=find(xmth == 1);
 yday(ind)=xday(ind);
 ind=find(xmth == 2);
 yday(ind)=31+xday(ind);
 xlat=alat(:,kyr);
 ind=find( xmth>2 & xmth<6 & xlat>70);
figure;hist(aNd(ind,kyr),0:10:350);xlim([50 320]);title(num2str(kyr+2007));
 end



%%%%%%%%%%%%%%%%%%%%%%




 for kyr=1:13
 % convert atm to #day 
 xyr=fix( atm(:,kyr) /10000 );
 xmth=fix( (atm(:,kyr)-xyr*10000) / 100 );
 xday=atm(:,kyr)-xyr*10000-xmth*100 ;
 ind=find(xmth == 1);
 yday(ind)=xday(ind);
 ind=find(xmth == 2);
 yday(ind)=31+xday(ind);

 xlat=alat(:,kyr);

ind=find( xmth>3 & xmth<6 & xlat>70);
figure;hist(aNd(ind,kyr),0:5:500);xlim([0 500]);title(num2str(kyr+2007));


 if ( fix(xyr/4)*4 == xyr)
     ind=find(xmth == 3);
     yday(ind)=60+xday(ind);
     ind=find(xmth == 4);
     yday(ind)=91+xday(ind);
     ind=find(xmth == 5);
     yday(ind)=121+xday(ind);
     ind=find(xmth == 6);
     yday(ind)=152+xday(ind);
     ind=find(xmth == 7);
     yday(ind)=182+xday(ind);
     ind=find(xmth == 8);
     yday(ind)=213+xday(ind);
     ind=find(xmth == 9);
     yday(ind)=244+xday(ind);
     ind=find(xmth == 10);
     yday(ind)=274+xday(ind);
     ind=find(xmth == 11);
     yday(ind)=305+xday(ind);
     ind=find(xmth == 12);
     yday(ind)=335+xday(ind);
 else
     ind=find(xmth == 3);
     yday(ind)=59+xday(ind);
     ind=find(xmth == 4);
     yday(ind)=90+xday(ind);
     ind=find(xmth == 5);
     yday(ind)=120+xday(ind);
     ind=find(xmth == 6);
     yday(ind)=151+xday(ind);
     ind=find(xmth == 7);
     yday(ind)=181+xday(ind);
     ind=find(xmth == 8);
     yday(ind)=212+xday(ind);
     ind=find(xmth == 9);
     yday(ind)=243+xday(ind);
     ind=find(xmth == 10);
     yday(ind)=273+xday(ind);
     ind=find(xmth == 11);
     yday(ind)=304+xday(ind);
     ind=find(xmth == 12);
     yday(ind)=334+xday(ind);
 end

 zday(1:length(yday'),kyr) = yday';
 % find mean Nd of the Arctic for yday = k 
 mmm=1:1:length(yday);


 for kday=1:365
     ylat=alat(mmm,kyr)'; yNd=aNd(mmm,kyr)'; yalt=aalt(mmm,kyr)';
     ydp=adp(mmm,kyr)'; yre=are(mmm,kyr)';
     ind=find( (ylat > 70) & (fix(yday) == kday) & (yalt < 8) );
     meanNd(kday, kyr) = exp( nanmean(log(aNd(ind,kyr))) );
     meanNdext(kday, kyr) = nanmean((aext(ind,kyr)/5.6).^2);
      meandp(kday, kyr) = nanmean(adp(ind,kyr) );
      meanre(kday, kyr) = nanmean(are(ind,kyr) );
     meanext(kday, kyr) = nanmean(aext(ind,kyr) );
     meandpall(kday, kyr) = nanmean(adp(:,kyr) );
     meanalt(kday, kyr) = nanmean(aalt(ind,kyr) );
     ind2=find( (ylat > 70) & (fix(yday) == kday) & (ydp<0.15) & (yre>7) & (yalt < 10)  );
     Ndg200(kday, kyr) = length(ind2) / length(ind);
     nn2(kday, kyr)=length(ind2);
     nn1(kday, kyr)=length(ind);
 end

 kyr
 end

 meannewdp=((meanext+18)/22).^5;
 meannewdp2=exp(meandp/0.055);


 %%%%%%%%

tmplat=reshape(alat,1,[]);
tmpnd=reshape(aNd,1,[]);
tmpmth=reshape(axmth,1,[]);
tmpdp=reshape(adp,1,[]);
tmpext=reshape(aext,1,[]);
tmpalt=reshape(aalt,1,[]);
tmpre=reshape(are,1,[]);
tmplwc=reshape(alwc,1,[]);
tmpalt=reshape(aalt,1,[]);

newnd=exp(tmpdp/0.065)+10;
newaNd=exp(adp/0.065)+10;

ind=find(tmplat>70);

xxx=15:5:400; 
myh3dpeak2(tmpnd(1,ind),tmpdp(1,ind),5:5:400,0.05:0.01:0.45);
hold on; plot(xxx, log(xxx-10) * 0.065)

myh3dpeak2(tmpdp(1,ind),tmpext(1,ind),0.05:0.01:0.45,5:2:120);
xlabel('\delta');ylabel('Extinction');

myh3dpeak2(tmpre(1,ind),tmpext(1,ind),4:0.2:20,5:2:120);
xlabel('Re');ylabel('Extinction');

myh3dpeak2(tmpre(1,ind),tmpdp(1,ind),4:0.2:20,0:0.01:0.45);
xlabel('Re');ylabel('Depolarization');

myh3dpeak2(tmpre(1,ind),newnd(1,ind),4:0.2:20,0:3:400);
xlabel('Re');ylabel('New Nd');

myh3dpeak2(tmpre(1,ind),tmpnd(1,ind),4:0.2:20,0:3:400);
xlabel('Re');ylabel('Nd');

myh3dpeak2(tmplwc(1,ind),tmpalt(1,ind),0.05:0.01:0.45,0:0.2:8);
xlabel('LWC');ylabel('Altitude (km)');

myh3dpeak2(tmpnd(1,ind),tmpalt(1,ind),5:3:300,0:0.2:8);
xlabel('Nd');ylabel('Altitude (km)');

myh3dpeak2(tmpre(1,ind),tmpalt(1,ind),4:0.3:20,0:0.2:8);
xlabel('Re');ylabel('Altitude (km)');

myh3dpeak2(tmpext(1,ind),tmpalt(1,ind),5:2:60,0:0.2:8);
xlabel('Extinction');ylabel('Altitude (km)');


myh3d(tmpdp(1,ind),tmpnd(1,ind),0.05:0.005:0.4,0:1:350);
xx=0:0.01:0.4;yy= 5+0.5*exp(xx*20)  ;
hold on;plot(xx,yy,'.-k');


figure;hist(tmpnd,0:1:300);xlim([0 300]);

figure;hist(tmpdp,0:0.01:0.45);


m3=1
n8=30
mmNd(:,1)=(movmean( (meandp(1:365,2).^m3+meandp(1:365,6).^m3+meandp(1:365,7).^m3)/3, n8, 'omitnan')).^(1/m3);
mmNd(:,2)=(movmean( (meandp(1:365,4).^m3+meandp(1:365,8).^m3+meandp(1:365,3).^m3+meandp(1:365,1).^m3)/4, n8, 'omitnan')).^(1/m3);
mmNd(:,3)=( movmean(meandp(1:365,5),n8,'omitnan') ) ;
figure('position',[100 100 450 300]);
plot(mmNd(:,1),'.-r');
hold on;plot(mmNd(:,2),'.-g');
ylim([0.19 0.25]);
xlim([50 300]);
plot(mmNd(:,3),'.-b');
xlabel('Time: Day # in a Year');
ylabel('Cloud Depolarization Ratios');
legend('2009,2013,2014', '2008,2010,2011,2015', '2012', 'location','northeast'); 
saveas(gcf,'new1_dp_vs_sea_ice.png','png');

m3=1
n8=30
mmNd(:,1)=(movmean( (meanext(1:365,7).^m3+meanext(1:365,6).^m3+meanext(1:365,7).^m3)/3, n8, 'omitnan')).^(1/m3);
mmNd(:,2)=(movmean( (meanext(1:365,4).^m3+meanext(1:365,8).^m3+meanext(1:365,3).^m3+meanext(1:365,1).^m3)/4, n8, 'omitnan')).^(1/m3);
mmNd(:,3)=( movmean(meanext(1:365,5),n8,'omitnan') ) ;
mmNd(:,4)=( movmean(meanext(1:365,9),n8,'omitnan') ) ;
figure('position',[100 100 450 300]);
plot(mmNd(:,1),'.-r');
hold on
xlim([0 273]);
plot(mmNd(:,3),'.-b');
 plot(mmNd(:,4),'.-k');
xlabel('Time: Day # in a Year');
ylabel('Cloud Extinction Coefficint (1/km)');
legend('2009,2013,2014', '2012', '2016', 'location','southeast'); 
saveas(gcf,'new2_Ext_vs_sea_ice.png','png');


n8=30
mmNd(:,1)=movmean( (nn1(1:365,2)+nn1(1:365,6)+nn1(1:365,7))/3, n8, 'omitnan');
mmNd(:,2)=movmean( (nn1(1:365,4)+nn1(1:365,8)+nn1(1:365,3)+nn1(1:365,1))/4, n8, 'omitnan');
mmNd(:,3)=( movmean(nn1(1:365,5),n8,'omitnan') +movmean(meanNd(1:365,5),n8,'omitnan') )/2 ;
mmNd(:,4)=( movmean(nn1(1:365,9),n8,'omitnan') +movmean(meanNd(1:365,9),n8,'omitnan') )/2 ;
figure('position',[100 100 450 300]);
plot(mmNd(:,1),'.-r');
hold on;plot(mmNd(:,2),'.-g');
% ylim([50 225]);
xlim([0 365]);
plot(mmNd(:,3),'.-b');
plot(mmNd(:,4),'.-k');
xlabel('Time: Day # in a Year5');
ylabel('# of Observations per day');
legend('2009,2013,2014', '2008,2010,2011,2015', '2012', '2016'); 
saveas(gcf,'nn_Nd_vs_sea_ice.png','png');


load('algaedailynew.mat');
figure('position',[100 100 500 300]);
yyaxis left
plot(timealgae,movmean(meanalgae,7),'.-b','linewidth',3);xlim([0 365]);
ylim([0 27]);
xlabel('Number of Days since January 1st');
ylabel('Ice Algal Biomass (mg/m^3)');
yyaxis right;
ind=[1 3 4 5 8 9 10];
for kday=1:365
    nndd(kday)=nanmean(meanNd(kday,ind));
end
plot(1:365,movmean(nndd,15,'omitnan'),'.-r','linewidth',2);ylim([40 160]);
ylabel('Droplet Number Concentration (cm^{-3})');
saveas(gcf,'daily_icealgae_Nd.png');


 for k=1:13
     tt(1:365,k)=(1:365)/365+2007+k;
 end

 mmext=movmean(meanext,15,1,'omitnan');
 mmdp=movmean(meandp,15,1,'omitnan');
 c1=nanmean(mmext,2);
 c2=nanmean(mmdp,2);

 for k=1:13
     extanomaly(:,k)=mmext(:,k)-c1;
     dpanomaly(:,k)=mmdp(:,k)-c2;
 end


seaice=[4.59,5.12,4.62,4.34,3.39,5.05,5.03,4.43,4.17,4.67,4.66,4.19,3.82];
sea0=14.09;

bNd=exp(adp/0.065)+10;
for k=1:13  
ind=find(alat(:,k)>70 & zday(:,k)>90 & zday(:,k)<300 & aalt(:,k)<4 & aalt(:,k)>0.23);
malt(k)=nanmean(aalt(ind,k));
mdp(k)=nanmean(adp(ind,k));
mbkr(k)=nanmean(abkr(ind,k));
mbkg(k)=nanmean(abkg(ind,k));
mcr(k)=nanmean(abkg(ind,k)./abkr(ind,k));
mnd(k)=nanmean(aNd(ind,k));
% mnd13(k)=nanmean(aNd(ind,k).^0.333333);
mlognd(k)=nanmean(log(bNd(ind,k)));
mre(k)=nanmean(are(ind,k));
mlwc(k)=nanmean(alwc(ind,k));
mext(k)=nanmean(aext(ind,k));
nnn(k)=length(ind);
stddp(k)=std(adp(ind,k));
stdext(k)=std(aext(ind,k));
stdnd(k)=std(aNd(ind,k));
stdlwc(k)=std(adp(ind,k));
% mdp2=nanmean(meandp);
errdp=stddp./sqrt(nnn);
errext=stdext ./ sqrt(nnn);
errnd=stdnd ./ sqrt(nnn);
end

mbnd=exp(mlognd);

for k=1:13
ind=find(alat(:,k)>70 & zday(:,k)>0 & zday(:,k)<365 & aalt(:,k)>0.);
indam=find(alat(:,k)>70 & zday(:,k)>90 & zday(:,k)<260 & aalt(:,k)>0.);
malt(k)=nanmean(aalt(ind,k));
mdp(k)=nanmean(adp(ind,k));
mdpam(k)=nanmean(adp(indam,k));
mbkr(k)=nanmean(abkr(ind,k));
mbkg(k)=nanmean(abkg(ind,k));
mnd(k)=nanmean(aNd(ind,k));
mre(k)=nanmean(are(ind,k));
mream(k)=nanmean(are(indam,k));
mlwc(k)=nanmean(alwc(ind,k));
mlwcam(k)=nanmean(alwc(indam,k));
mext(k)=nanmean(aext(ind,k));
mextam(k)=nanmean(aext(indam,k));
msc(k)=nanmean( 0.5* (1+adp(indam,k)).^2 ./(abkg(indam,k).*(1-adp(indam,k)).^2) );
nnn(k)=length(ind);
stddp(k)=std(adp(ind,k));
stdext(k)=std(aext(ind,k));
stdnd(k)=std(aNd(ind,k));
stdbkg(k)=std(abkg(ind,k));
stdre(k)=std(are(ind,k));
stdlwc(k)=std(alwc(ind,k));
end

errdp=stddp./sqrt(nnn);
errext=stdext ./ sqrt(nnn);
errnd=stdnd ./ sqrt(nnn);
errbkg=stdbkg ./ sqrt(nnn);
errre=stdre ./ sqrt(nnn);
errlwc=stdlwc ./ sqrt(nnn);

n10=10;
figure('position',[100 100 500 300]);
yyaxis left
errorbar(2007+(1:n10)+0.65,mdp(1:n10), errdp(1:n10), '.-b','markersize',25, 'linewidth',3);
xtmp=5;
mtmp=mean(mdp);
nmdpam=mtmp+(mdpam-mtmp)/xtmp;
%  hold on; plot(2007+(1:n10)+0.65,nmdpam(1:10),'.-g','markersize',12);
xlabel('Time (June-Mid-Sept. for Clouds; Mid-Sept. for Sea-Ice)');ylabel('Arctic Water Cloud Depolariation');
xticks([2008 2010, 2012, 2014, 2016, 2018]); 
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
 ylim([0.218 0.228]);
%  ylim([0.217 0.244]);
yyaxis right
hold on; plot(2007+(1:n10)+0.7,seaice(1:n10),'.--r','markersize',25, 'linewidth',2);
ylabel('Minimum Sea Ice Extent (million km^2)');
ylim([3.2 5.55]);
saveas(gcf,'depo-sea-ice-paper-2025.png','png')

xfg2b(1,:)=2008:1:2017;
xfg2b(2,:)=mdp(1:10);
xfg2b(3,:)=nmdpam(1:10);
xfg2b(4,:)=seaice(1:10);

tau=exp(10*mdp-0.38);
a1=0.58;
b1=0.74;
b2=-0.2053;
b3=-0.7935;
k1=2.0506;
k2=0.1968;
k3=2.479;
c=0.1409;
d=0.1329;
at=a1+(1-a1)*exp(-k1*tau);
bt=b1*(1+b2*exp(-k2*tau)+b3*exp(-k3*tau));
theta=60;
alpha=0.4;
trc=(at+bt*cosd(theta)) ./ (1+(c-d*alpha)*tau);
e1=(trc-trc(2)) * 5.02/0.015 ;
e2=6.8/1.5 *(seaice(2)-seaice);


figure('position',[100 100 500 300]);
yyaxis left
plot(2007+(1:n10)+0.65,tau(1:n10), '.-b','markersize',25, 'linewidth',3);
% hold on; plot(2007+(1:10),mdp2(1:10),'.-g','markersize',12);
xlabel('Year');
ylabel('Water Cloud \tau');
xticks([2008 2010, 2012, 2014, 2016, 2018]);
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
  ylim([6.05 6.65]);
%  ylim([0.218 0.2275]);
yyaxis right
hold on; plot(2007+(1:n10)+0.65,trc(1:n10),'.--r','markersize',25, 'linewidth',2);
ylabel('Transmittance');
  ylim([0.586 0.606]);
set(gca,'Ydir','reverse')
saveas(gcf,'figS1-2025.png','png')

xfgs1(1,:)=2008:1:2017;
xfgs1(2,:)=tau(1:10);
xfgs1(3,:)=trc(1:10);
save('figs1.txt','xfgs1','-ascii');

n10=10;
figure('position',[100 100 500 300]);
yyaxis left
plot(2007+(1:n10)+0.65,e1(1:n10), '.-b','markersize',25, 'linewidth',3);
% hold on; plot(2007+(1:10),mdp2(1:10),'.-g','markersize',12);
xlabel('Time June-Mid-Sept. for Clouds; Mid-Sept. for Sea-Ice');
ylabel('Tr(sunlight) - Tr(2009) (10^{20} J)');
xticks([2008 2010, 2012, 2014, 2016, 2018]);
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
  ylim([-2 10]);
yyaxis right
hold on; plot(2007+(1:n10)+0.7,e2(1:n10),'.--r','markersize',25, 'linewidth',2);
ylabel('Melting - Melting(2009)  (10^{20} J)');
  ylim([-2 10]);
saveas(gcf,'fig1c-2025.png','png')

xfg2c(1,:)=2008:1:2017;
xfg2c(2,:)=e1(1:10);
xfg2c(3,:)=e2(1:10);


figure('position',[100 100 500 300]);
yyaxis left
plot(2007+(1:n10)+0.5,tau(1:10), '.-b','markersize',25, 'linewidth',3);
xlabel('Time');
ylabel('\tau');
xticks([2008 2010, 2012, 2014, 2016, 2018]);
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
   ylim([6.0 6.8]);
%  ylim([0.218 0.2275]);
yyaxis right
hold on; plot(2007+(1:n10)+0.5,trc(1:n10),'.--r','markersize',25, 'linewidth',2);
ylabel('Transmittance');
set(gca,'Ydir','reverse')
   ylim([0.58 0.61]);
saveas(gcf,'tau-transmittance.png','png')


figure('position',[100 100 500 300]);
yyaxis left
errorbar(2007+(1:n10)+0.65,mext(1:n10), errext(1:n10), '.-b','markersize',25, 'linewidth',3);
xtmp=5;
mtmp=mean(mextam);
nmextam=mtmp+(mextam-mtmp)/xtmp;
%  hold on; plot(2007+(1:10)+0.65,nmextam(1:10),'.-g','markersize',12);
xlabel('Time (June- Mid-Sept. for Clouds; Mid-Sept. for Sea-Ice)');ylabel('Arctic Water Cloud Extinction (1/km)');
ylim([35. 38.]);
xticks([2008 2010, 2012, 2014, 2016, 2018]); 
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
yyaxis right
hold on; plot(2007+(1:n10)+0.7,seaice(1:n10),'.--r','markersize', 25, 'linewidth',2);
ylabel('Minimum Sea Ice Extent (million km^2)');
ylim([3.2 5.55]);
saveas(gcf,'ext-sea-ice-paper-2025.png','png')

xfg2a(1,:)=2008:1:2017;
xfg2a(2,:)=mext(1:10);
xfg2a(3,:)=nmextam(1:10);
xfg2a(4,:)=seaice(1:10);

save('fig2a.txt','xfg2a','-ascii');
save('fig2b.txt','xfg2b','-ascii');
save('fig2c.txt','xfg2c','-ascii');


figure('position',[100 100 500 300]);
yyaxis left
plot(2007+(1:n10)+0.65,tau(1:n10),  '.-b','markersize',25, 'linewidth',3);
ylim([6 6.7]);
%  hold on; plot(2007+(1:10)+0.65,nmextam(1:10),'.-g','markersize',12);
xlabel('Time');ylabel('Cloud Optical Depth');
xticks([2008 2010, 2012, 2014, 2016, 2018]); 
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
yyaxis right
tmplwp=tau*309.0 ./mext;
hold on; plot(2007+(1:n10)+0.7,tmplwp(1:n10),'.--r','markersize', 25, 'linewidth',2);
ylabel('Cloud Liquid Water Path  (g/m^{2})');
ylim([51 57]);
saveas(gcf,'tau_lwp.png','png')

figure('position',[100 100 500 300]);
yyaxis left
plot(2007+(1:n10)+0.65,trc(1:n10),  '.-b','markersize',25, 'linewidth',3);
% ylim([6 6.7]);
%  hold on; plot(2007+(1:10)+0.65,nmextam(1:10),'.-g','markersize',12);
xlabel('Time');ylabel('Cloud Transmittance');
xticks([2008 2010, 2012, 2014, 2016, 2018]); 
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
yyaxis right
emis=1-exp(-0.158*tmplwp)
hold on; plot(2007+(1:n10)+0.7,emis(1:n10),'.--r','markersize', 25, 'linewidth',2);
ylabel('Cloud Emissivity');
 ylim([0.99 1.01]);
saveas(gcf,'trc_emis.png','png')


n10=10;
figure('position',[100 100 500 300]);
yyaxis left
plot(2007+(1:n10)+0.5,malt(1:n10), '.-b','markersize',25, 'linewidth',3);
% hold on; plot(2007+(1:10),mdp2(1:10),'.-g','markersize',12);
xlabel('Time (Apr.-Sept. for Clouds; Fall for Sea-Ice');ylabel('Water Cloud Top Height (km)');
% ylim([1.95 2.4]);
set(gca,'Ydir','reverse')
% ylim([35 38.]);
xticks([2008 2010, 2012, 2014, 2016, 2018]); 
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
yyaxis right
hold on; plot(2007+(1:n10)+0.75,seaice(1:n10),'.--r','markersize', 25, 'linewidth',2);
ylabel('Minimum Sea Ice Extent (million km^2)');
ylim([3.2 5.6]);
saveas(gcf,'whole-year-cloudhit-sea-ice-paper.png','png')

for kk=1:13
nyr=kk;
ind=find(alat(:,nyr)>70 & zday(:,nyr)>90 & zday(:,nyr)<150 & aalt(:,nyr)>0);
yearmn(kk) = mean(abkg(ind,nyr)) ;
end
n10=10;
figure;plot(2007+(1:n10),yearmn(1:n10),'.-b','markersize',12);
yyaxis right
hold on; plot(2007+(1:n10),seaice(1:n10),'.--r','markersize', 25, 'linewidth',2);
myh3dpeak2(reshape(adp,1,[]), reshape(aext,1,[]),0.05:0.01:0.45,10:1:128);
xlabel('Water Cloud Depolarization Ratio'); ylabel('Extinction Coefficient (km^{-1})');
saveas(gcf,'depo-extinction-correlation.png','png')


ind=find(aext>9);
myh3dpeak2(reshape(abkg,1,[]), reshape(aext,1,[]),0.05:0.002:0.125,15:1:108);
xlabel('CALIPSO 532 nm Water Cloud Reflectance (Sr^{-1})'); ylabel('Extinction Coefficient (km^{-1})');
saveas(gcf,'reflectance-extinction-correlation.png','png')


acal(1:5000000,1:13)=NaN;
eta(1:5000000,1:13)=NaN;
scc(1:5000000,1:13)=NaN;
ddd(1:5000000,1:13)=NaN;
arcbkg(1:5000000,1:13)=NaN;
arcdp(1:5000000,1:13)=NaN;
for kyr=1:13
ind=find(alat(:,kyr)>70);
acal(1:length(ind),kyr)=abkg(ind,kyr)./adp(ind,kyr);
eta(1:length(ind),kyr)=((1-adp(ind,kyr))./(1+adp(ind,kyr))).^2;
scc(1:length(ind),kyr)=abkg(ind,kyr).*eta(1:length(ind),kyr);
ddd(1:length(ind),kyr)=zday(ind,kyr);
arcbkg(1:length(ind),kyr)=abkg(ind,kyr);
arcdp(1:length(ind),kyr)=adp(ind,kyr);
end


figure;plot(ddd(:,1),movmean(acal(:,1),3000,'omitnan'),'.r')
hold on; plot(ddd(:,2),movmean(acal(:,2),3000,'omitnan'),'.g')
plot(ddd(:,3),movmean(acal(:,3),3000,'omitnan'),'.b')
plot(ddd(:,4),movmean(acal(:,4),3000,'omitnan'),'.k')
plot(ddd(:,5),movmean(acal(:,5),3000,'omitnan'),'.c')


figure;plot(ddd(:,1),movmean(scc(:,1),3000,'omitnan'),'.r')
hold on; plot(ddd(:,2),movmean(scc(:,2),3000,'omitnan'),'.g')
plot(ddd(:,3),movmean(scc(:,3),3000,'omitnan'),'.b')
plot(ddd(:,4),movmean(scc(:,4),3000,'omitnan'),'.k')
plot(ddd(:,5),movmean(scc(:,5),3000,'omitnan'),'.c')

figure;plot(ddd(:,1),movmean(abkg(:,1),10000,'omitnan'),'.r')
hold on; plot(ddd(:,2),movmean(abkg(:,2),10000,'omitnan'),'.g')
plot(ddd(:,3),movmean(abkg(:,3),10000,'omitnan'),'.b')
plot(ddd(:,4),movmean(abkg(:,4),10000,'omitnan'),'.k')
plot(ddd(:,5),movmean(abkg(:,5),10000,'omitnan'),'.c')

figure;plot(ddd(:,1),movmean(adp(:,1),10000,'omitnan'),'.r')
hold on; plot(ddd(:,2),movmean(adp(:,2),10000,'omitnan'),'.g')
plot(ddd(:,3),movmean(adp(:,3),10000,'omitnan'),'.b')
plot(ddd(:,4),movmean(adp(:,4),10000,'omitnan'),'.k')
plot(ddd(:,5),movmean(adp(:,5),10000,'omitnan'),'.c')

figure;plot(ddd(:,1),movmean(arcbkg(:,1),3000,'omitnan'),'.r')
hold on; plot(ddd(:,2),movmean(arcbkg(:,2),3000,'omitnan'),'.g')
plot(ddd(:,3),movmean(arcbkg(:,3),3000,'omitnan'),'.b')
plot(ddd(:,4),movmean(arcbkg(:,4),3000,'omitnan'),'.k')
plot(ddd(:,5),movmean(arcbkg(:,5),3000,'omitnan'),'.c')

figure;plot(ddd(:,1),movmean(arcdp(:,1),3000,'omitnan'),'.r')
hold on; plot(ddd(:,2),movmean(arcdp(:,2),3000,'omitnan'),'.g')
plot(ddd(:,3),movmean(arcdp(:,3),3000,'omitnan'),'.b')
plot(ddd(:,4),movmean(arcdp(:,4),3000,'omitnan'),'.k')
plot(ddd(:,5),movmean(arcdp(:,5),3000,'omitnan'),'.c')



for j1=1:10
for k1=1:7
ndays=(30*(k1-1)+91):1:(k1*30+90);
newNd(k1,j1) = nanmean( meanNd(ndays,j1) ); 
end 
end 
 

allday=reshape(zday(:,1:13),1,[]);
alllat=reshape(alat,1,[]);
alldp=reshape(adp,1,[]);
allext=reshape(aext,1,[]);
allNd=reshape(aNd,1,[]);
allalt=reshape(aalt,1,[]);

ind=find(alllat(1,:)>70 & allday(1,:)>90 & allday(1,:)<150 & allalt(1,:)<3  & allalt(1,:)>0.3);
dplim=0.05:0.01:0.45;altlim=0:0.15:3; ndlim=0:10:300;
myh3dpeak2(allNd(1,ind),alldp(1,ind),ndlim,dplim);
xlabel('Nd in April/May');
ylabel('Depolarization'); 


dplim=0.05:0.01:0.45;altlim=0:0.15:3; ndlim=0:10:300;
myh3d(allNd(1,ind),allalt(1,ind),ndlim,altlim);
xlabel('Nd in April and May');
ylabel('Cloud Top Height (km)'); 

myh3dpeak2(allNd(1,ind),allalt(1,ind),ndlim,altlim);
xlabel('Nd in April and May');
ylabel('Cloud Top Height (km)'); 

dplim=0.1:0.005:0.3;altlim=0:0.15:3; ndlim=0:10:300;
myh3dpeak2(alldp(1,ind),allalt(1,ind),dplim,altlim);
xlabel('Depolarization in April and May');
ylabel('Cloud Top Height (km)'); 


ind=find(alllat(1,:)>70 & allday(1,:)>0 & allday(1,:)<300 & allalt(1,:)<8 & allalt(1,:)>0.3);
dplim=0.05:0.01:0.45;altlim=0:0.15:5; ndlim=0:10:300;
myh3d(allNd(1,ind),allalt(1,ind),ndlim,altlim);
xlabel('Nd in June/July/August/September');
ylabel('Cloud Top Height (km)');

myh3dpeak2(allNd(1,ind),allalt(1,ind),ndlim,altlim);
xlabel('Nd in June/July/August/September');
ylabel('Cloud Top Height (km)'); 

dplim=0.1:0.005:0.37;altlim=0:0.15:7.2; ndlim=0:10:300;
myh3d(alldp(1,ind),allalt(1,ind),dplim,altlim);
xlabel('Arctic Water Cloud Depolarization Ratios');
ylabel('Cloud Top Height (km)'); 
saveas(gcf,'cloud-altitude-dep.png','png');

ind=find( allday(1,:)>0 & allday(1,:)<370 & allalt(1,:)<4 & allalt(1,:)>0.3);
dplim=0.1:0.005:0.42;altlim=0:0.15:3; ndlim=0:3:300;
myh3dpeak2(allNd(1,ind),alldp(1,ind),ndlim,dplim);
xlabel('Water Cloud Drop Number Density Nd (cm^{-3})');
ylabel('Depolarization Ratios'); 
saveas(gcf,'nd-depo-relation.png','png');


ind=find(adp(:,1)>0.09  & adp(:,1)<0.45);
myh3dpeak2(reshape(abkg(ind,1),1,[]),reshape(adp(ind,1),1,[]), 0.05:0.002:0.12,0.08:0.01:0.45);
% myh3dpeak2(reshape(abkg,1,[]),reshape(adp,1,[]), 0.05:0.002:0.12,0.08:0.01:0.45);
yy=0:0.1:0.5; 
xx=yy/4.25+0.015; 
hold on; plot(xx,yy,'--g','linewidth',1);
xlabel('CALIPSO Water Cloud 532 nm Reflectance \beta (Sr^{-1})');
ylabel('Cloud Depolarization Ratio \delta');
legend('\delta = 4*(\beta - 0.01)','location','southeast'); 
saveas(gcf,'depo-reflectance-relation-532 nm.png','png');

ind=find(adp(:,1)>0.09  & adp(:,1)<0.4);
myh3dpeak2(reshape(abkr(ind,1),1,[]),reshape(adp(ind,1),1,[]), 0.05:0.002:0.12,0.08:0.01:0.45);
% myh3dpeak2(reshape(abkg,1,[]),reshape(adp,1,[]), 0.05:0.002:0.12,0.08:0.01:0.45);
yy=0:0.1:0.5; 
xx=yy/3.5+0.018; 
hold on; plot(xx,yy,'--y','linewidth',3);
xlabel('Water Cloud 1064 nm Reflectance \beta (Sr^{-1})');
ylabel('Cloud Depolarization Ratio \delta');
legend('\delta = 3.5*(\beta - 0.018)','location','southeast'); 
saveas(gcf,'depo-reflectance-relation-1064 nm.png','png');


n8=30
mmNd(:,1)=movmean( (meanNd(1:365,2)+meanNd(1:365,6)+meanNd(1:365,7))/3, n8, 'omitnan');
mmNd(:,2)=movmean( (meanNd(1:365,4)+meanNd(1:365,8)+meanNd(1:365,3)+meanNd(1:365,1))/4, n8, 'omitnan');
mmNd(:,3)=0.9*( movmean(meanNd(1:365,5),n8,'omitnan') +movmean(meanNd(1:365,5),n8,'omitnan') )/2 ;
mmNd(:,5)=movmean(meanNd(1:365,9),n8,'omitnan') ;
mmNd(:,6)=movmean(meanNd(1:365,10),n8,'omitnan') ;

mmNdext(:,1)=movmean( (meanNdext(1:365,2)+meanNdext(1:365,6)+meanNdext(1:365,7))/3, n8, 'omitnan');
mmNdext(:,2)=movmean( (meanNdext(1:365,4)+meanNdext(1:365,8)+meanNdext(1:365,3)+meanNdext(1:365,1))/4, n8, 'omitnan');

bew(:,2)=movmean(meanNdext(1:365,2),n8,'omitnan') .* mmNd(:,1) ./ mmNdext(:,1);
bew(:,6)=movmean(meanNdext(1:365,6),n8,'omitnan') .* mmNd(:,1) ./ mmNdext(:,1);
bew(:,7)=movmean(meanNdext(1:365,7),n8,'omitnan') .* mmNd(:,1) ./ mmNdext(:,1);

bew(:,4)=movmean(meanNdext(1:365,4),n8,'omitnan') .* mmNd(:,2) ./ mmNdext(:,2);
bew(:,8)=movmean(meanNdext(1:365,8),n8,'omitnan') .* mmNd(:,2) ./ mmNdext(:,2);
bew(:,3)=movmean(meanNdext(1:365,3),n8,'omitnan') .* mmNd(:,2) ./ mmNdext(:,2);
bew(:,1)=movmean(meanNdext(1:365,1),n8,'omitnan') .* mmNd(:,2) ./ mmNdext(:,2);
bew(:,5)=mmNd(:,3);
bew(:,9)=mmNd(:,5);
bew(:,10)=mmNd(:,6);



ind=find(tmplat>70);
myh3d(tmpnd(ind),tmpre(ind),5:3:400,4:0.3:20)
xxx=4:0.1:25; yyy=24000 ./xxx.^2.5;
hold on;plot(yyy,xxx,'--g','linewidth',2)
legend('Nd=24000/Re^{5/2}');
ylabel('Re (\mum)');xlabel('Nd (cm^{-3})');
saveas(gcf,'nd-re-relation.png','png');

ind=find(tmplat>70 & tmpmth>5 & tmpmth<9);
n30=3;
avg1=movmean(tmpnd(ind),n30);
avg2=movmean(tmpext(ind),n30);
myh3d(avg1.^0.2,avg2,1.35:0.03:3.25,5:0.3:70)
% myh3d(tmpnd(ind).^0.2,tmpext(ind),1.35:0.03:3.25,5:0.3:70)
xxx=1:0.1:3.5; yyy=22*xxx -18 ;
 hold on;plot(xxx,yyy,'--g','linewidth',2)
legend('ext=22*Nd^{0.2}-18','location','southeast');
ylabel('extinction coefficient');xlabel('Nd^{0.2}');
saveas(gcf,'new-nd-ext-relation-JJA.png','png');


ind=find(tmplat>70);
myh3dpeak2(tmpext(ind),tmpre(ind),5:2:120,4:0.3:20)
ylabel('Re (\mum)');xlabel('ext (km^{-1})');
saveas(gcf,'ext-re-relation.png','png');


nnn1=nn2./nn1;

n8=30
mmNd(:,1)=movmean( (nnn1(1:365,2)+nnn1(1:365,6)+nnn1(1:365,7)+nnn1(1:365,10))/4, n8, 'omitnan');
mmNd(:,2)=movmean( (nnn1(1:365,4)+nnn1(1:365,8)+nnn1(1:365,3)+nnn1(1:365,1))/4, n8, 'omitnan');
mmNd(:,3)=( movmean(nnn1(1:365,9),n8,'omitnan') +movmean(nnn1(1:365,5),n8,'omitnan') )/2 ;
figure('position',[100 100 450 300]);
plot(mmNd(:,1),'.-r');
hold on;plot(mmNd(:,2),'.-g');
plot(mmNd(:,3),'.-b');
xlabel('Time: Day 1 to Day 365');
ylabel('Probability (Nd > 200 cm^{-1})');
legend('2009,2013,2014,2017', '2008,2010,2011,2015', '2012,2016'); 
 saveas(gcf,'large_Nd_probability.png','png');


clear winteralt;clear winternd;
for k=1:13
        ind=find(alat(:,k)>70 & zday(:,k)>0 & zday(:,k)<90 & aalt(:,k)>0.1  & aalt(:,k)<8) ;
    if k>1 
ind1=find(alat(:,k-1)>70 & zday(:,k-1)>300 & zday(:,k-1)<370 & aalt(:,k-1)>0.1  & aalt(:,k-1)<3);
    end

if k==1 
    winteralt(1,1:length(ind))=aalt(ind,k)';
    winternd(1,1:length(ind))=aNd(ind,k)';
    winterdp(1,1:length(ind))=adp(ind,k)';
else
    winteralt=[winteralt aalt(ind,k)'];
    winternd=[winternd aNd(ind,k)'];
    winterdp=[winterdp adp(ind,k)'];
end


malt(k)=nanmean(aalt(ind,k));
mdp(k)=nanmean(adp(ind,k));
mbkr(k)=nanmean(abkr(ind,k));
mbkg(k)=nanmean(abkg(ind,k));
mnd(k)=nanmean(aNd(ind,k));
mre(k)=nanmean(are(ind,k));
mlwc(k)=nanmean(alwc(ind,k));
mext(k)=nanmean(aext(ind,k));
nnn(k)=length(ind);
stddp(k)=std(adp(ind,k));
stdext(k)=std(aext(ind,k));
stdnd(k)=std(aNd(ind,k));
stdbkg(k)=std(abkg(ind,k));
stdre(k)=std(are(ind,k));
stdlwc(k)=std(alwc(ind,k));
end
% mdp2=nanmean(meandp);close all


figure('position',[100 100 500 300]);
yyaxis left
plot(2007+(1:n10)+0.1,malt(1:n10),'.-b','markersize',25, 'linewidth',3);
% hold on; plot(2007+(1:10),mdp2(1:10),'.-g','markersize',12);
xlabel('Time: (Jan.-Mar. for clouds, fall for sea ice)');
ylabel('Water Cloud Top Height (km)');
%   ylim([1.75 2.23]);
set(gca,'Ydir','reverse')
% ylim([35 38.]);
xticks([2008 2010, 2012, 2014, 2016, 2018]); 
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
yyaxis right
hold on; plot(2007+(1:n10)+0.75,seaice(1:n10),'.--r','markersize', 25, 'linewidth',2);
ylabel('Minimum Sea Ice Extent (million km^2)');
ylim([3.2 5.6]);
saveas(gcf,'cloudhit-sea-ice-paper.png','png')

n10=10;
load('chl.mat');
chl=mchl(6,:);
figure('position',[100 100 500 300]);
yyaxis left
plot(2007+(1:n10)+0.7,chl,'.-b','markersize',25, 'linewidth',3);
% hold on; plot(2007+(1:10),mdp2(1:10),'.-g','markersize',12);
xlabel('Time: (Aug. for Chl, fall for sea ice)');
ylabel('Chloophyll-a concentration (mg m^{-3})');
%   ylim([1.75 2.23]);
% set(gca,'Ydir','reverse')
% ylim([35 38.]);
xticks([2008 2010, 2012, 2014, 2016, 2018]); 
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
yyaxis right
hold on; plot(2007+(1:n10)+0.75,seaice(1:n10),'.--r','markersize', 25, 'linewidth',2);
ylabel('Minimum Sea Ice Extent (million km^2)');
ylim([3.2 5.6]);
saveas(gcf,'chlorophyll-sea-ice-paper.png','png')

chl=mean( mchl(6,1:10), 1) .* (14-seaice(1,1:10))/1.48  * 1.e6 * 22 *1e-9;
% chl=mean( mchl(6,1:10).* marea(6,1:10)*6.37^2  * 1.e6 * 22, 1);
figure('position',[100 100 500 300]);
yyaxis left
plot(2007+(1:n10)+0.7,chl,'.-b','markersize',25, 'linewidth',3);
% hold on; plot(2007+(1:10),mdp2(1:10),'.-g','markersize',12);
xlabel('Time: (Aug. for Chl-a, fall for sea ice)');
ylabel('Chloophyll-a Biomass ( Million Ton )');
%   ylim([1.75 2.23]);
% set(gca,'Ydir','reverse')
% ylim([35 38.]);
xticks([2008 2010, 2012, 2014, 2016, 2018]); 
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
yyaxis right
hold on; plot(2007+(1:n10)+0.75,seaice(1:n10),'.--r','markersize', 25, 'linewidth',2);
ylabel('Minimum Sea Ice Extent (million km^2)');
ylim([3.2 5.6]);
saveas(gcf,'chlorophyll-sea-ice-paper-2.png','png')

chl=mean( mchl(6,1:10), 1) .* (14-seaice(1,1:10))/1.48  * 1.e12 * 22 *1e-15;
% chl=mean( mchl(6,1:10).* marea(6,1:10)*6.37^2  * 1.e6 * 22, 1);
figure('position',[100 100 500 300]);
yyaxis left
plot(2007+(1:n10)+0.7,chl,'.-b','markersize',25, 'linewidth',3);
% hold on; plot(2007+(1:10),mdp2(1:10),'.-g','markersize',12);
xlabel('Time: (Aug. for Chl-a and Cloud)');
ylabel('Chloophyll-a Biomass ( Million Ton )');
    ylim([0.19 0.27]);
% set(gca,'Ydir','reverse')
% ylim([35 38.]);
xticks([2008 2010, 2012, 2014, 2016, 2018]); 
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
yyaxis right
hold on; plot(2007+(1:n10)+0.7,nanmean(meanNd(212:243,1:n10),1),'.--r','markersize', 25, 'linewidth',2);
ylabel('Cloud Droplet Number Density (cm^{-3})');
  ylim([45 145]);
saveas(gcf,'chlorophyll-Nd-paper-2.png','png')



myh3dpeak2(winternd,winteralt);
xlabel('Cloud Droplet Number Density (cm^{-3}): Jan. to Mar.');
ylabel('Cloud Top Height (km)');
saveas(gcf,'nd-hit-winter.png','png');

myh3d(winternd,winteralt);
xlabel('Cloud Droplet Number Density (cm^{-3}): Jan. to Mar.');
ylabel('Cloud Top Height (km)');
saveas(gcf,'nd-hit-winter-2.png','png');


for k=1:13
ind=find(alat(:,k)>70 & zday(:,k)>90 & zday(:,k)<273 & aalt(:,k)<4);
malt(k)=nanmean(aalt(ind,k));
mdp(k)=nanmean(adp(ind,k));
mbkr(k)=nanmean(abkr(ind,k));
mbkg(k)=nanmean(abkg(ind,k));
mnd(k)=nanmean(exp(adp(ind,k)/0.055));
mlognd(k)=nanmean(((aext(ind,k)+18)/22).^5);
mre(k)=nanmean(are(ind,k));
mlwc(k)=nanmean(alwc(ind,k));
mext(k)=nanmean(aext(ind,k));
nnn(k)=length(ind);
stddp(k)=std(adp(ind,k));
stdext(k)=std(aext(ind,k));
stdnd(k)=std(aNd(ind,k));
stdbkg(k)=std(abkg(ind,k));
stdre(k)=std(are(ind,k));
stdlwc(k)=std(alwc(ind,k));
end

mlognd=((mext(1,1:10)+18)/22).^5;
mnd=exp(mdp/0.055); 

load('chl.mat');
figure('position',[100 100 500 300]);
yyaxis left
chl=mchl(6,:);
plot(2007+(1:n10)+0.7,chl,'.-b','markersize',25, 'linewidth',2);
% chl=mchl(7,:);
% hold on; plot(2007+(1:n10)+0.7,chl,'.-g','markersize',25, 'linewidth',2);
% % hold on; plot(2007+(1:10),mdp2(1:10),'.-g','markersize',12);
xlabel('Time: (Aug. for Chlorophyll, April to Sept for Clouds)');
ylabel('Chloophyll-a concentration (mg m^{-3})');
   ylim([0.6 0.87]);
% set(gca,'Ydir','reverse')
xticks([2008 2010, 2012, 2014, 2016, 2018]); 
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
yyaxis right
hold on; plot(2007+(1:n10)+0.7,mlognd(1:n10),'.-r','markersize', 25, 'linewidth',3);
ylabel('Cloud Droplet Number Density (cm^{-1})');
  ylim([120 155]);
saveas(gcf,'chlorophyll-nd-paper.png','png')

load('chl.mat');
chl=mchl(6,:);
figure('position',[100 100 500 300]);
yyaxis left
plot(2007+(1:n10)+0.7,chl,'.-b','markersize',25, 'linewidth',3);
% hold on; plot(2007+(1:10),mdp2(1:10),'.-g','markersize',12);
xlabel('Time: (Aug. for Chlorophyll, April to Sept for Clouds)');
ylabel('Chloophyll-a concentration (mg m^{-3})');
%   ylim([1.75 2.23]);
% set(gca,'Ydir','reverse')
% ylim([35 38.]);
xticks([2008 2010, 2012, 2014, 2016, 2018]); 
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
yyaxis right
hold on; plot(2007+(1:n10)+0.7,mdp(1:n10),'.--r','markersize', 25, 'linewidth',2);
ylabel('Cloud Depolarization Ratios');
% ylim([3.2 5.6]);
saveas(gcf,'chlorophyll-dp-paper.png','png')


myh3dpeak2(log(tmpnd),tmpdp,1.9:0.02:6,0.1:0.0075:0.38)
xlabel('log(Nd) (cm^{-3})');
ylabel('Depolarization Ratio \delta');
% xxx=1.9:0.02:6; yyy=0.056*xxx;
% hold on;plot(xxx,yyy,'--g','linewidth',2);
% legend('\delta=0.056*log(Nd)','location','northwest');
saveas(gcf,'nd-dp-relation.png','png');

ind=find(tmpmth>5 & tmpmth<9);
n3=10;
myh3dpeak2(sqrt(movmean(tmpnd(ind),n3)),movmean(tmpext(ind),n3),2.5:0.3:20,5:2.5:100)
xlabel('{ Nd (cm^{-3}) }^{1/2}');
ylabel('Extinction coefficient (km^{-1})');
xxx=1.9:0.02:26; yyy=4.2*xxx+6;
hold on;plot(xxx,yyy,'--g','linewidth',2);
legend('Ext. = 5+4.3*Nd^{1/2}','location','southeast');
saveas(gcf,'nd-ext-relation-JJA.png','png');

myh3dpeak2(tmpnd.^0.333,tmpext,2:0.1:7.2,5:2:120)
xlabel('{ Nd (cm^{-3}) }^{1/3}');
ylabel('Extinction coefficient (km^{-1})');
xxx=1.9:0.02:26; yyy=17*(xxx-1.2);
hold on;plot(xxx,yyy,'--g','linewidth',2);
legend('Ext. = 17*Nd^{1/3}-1.2','location','southeast');
saveas(gcf,'nd-ext-relation-2.png','png');




ind=find(alat(:,5)>70 & zday(:,5)>90 & zday(:,5)<273);
lat1=alat(ind,5);
lon1=alon(ind,5);
zz1=adp(ind,5);
ee1=aext(ind,5);
dd1=(aext(ind,5)+18).^5 / 22.^5;
ind=find(alat(:,2)>70 & zday(:,2)>90 & zday(:,2)<273);
lat2=alat(ind,2);
lon2=alon(ind,2);
zz2=adp(ind,2);
ee2=aext(ind,2);
dd2=(aext(ind,2)+18).^5 / 22.^5;

mymap_polar_diff(lat1,lon1,zz1,lat2,lon2,zz2)
title('\delta(2009) - \delta(2012)');
saveas(gcf,'dp-diff-map.png','png');
mymap_polar_diff(lat1,lon1,dd1,lat2,lon2,dd2)
title('Nd(2009) - Nd(2012)');
saveas(gcf,'Nd-diff-map.png','png');



tau=exp((adp-0.05)/0.1);
myh3dpeak2(reshape(aNd.^(1/3),1,[]),reshape(tau,1,[]),1.5:0.1:8,0:0.1:20);
xlabel('Nd^{1/3} (cm^{-3})'); ylabel('\tau');
xx=0:1:8;
yy=2.8*xx - 4 ; 
hold on;plot(xx,yy,'--g','linewidth',2);
legend('\tau = 2.8 * Nd^{1/3} - 4','location','southeast')
saveas(gcf,'most-probable-Nd-tau-relation.png','png');

aanndd=exp(adp./0.056);
alwp = 5/9 * (25000./aanndd).^(2/5) .* exp( (adp-0.05) / 0.1); 
for kyr=1:10
ind=find(zday(:,kyr)>90 & zday(:,kyr)<273 & alat(:,kyr)>70);
meanrree(kyr)=nanmean(are(ind,kyr));
meanddpp(kyr)=nanmean(adp(ind,kyr));
meanlwp(kyr)=nanmean(alwp(ind,kyr));
end


myh3dpeak2(reshape(adp,1,[]),reshape(are,1,[]),0.05:0.005:0.35,6:0.1:25)
xlabel('CALIOP Depolarization Ratio');
ylabel('MODIS Re (\mum)'); 
saveas(gcf,'depo_Re_relation.png','png')

myh3dpeak2(reshape(abkr./abkg,1,[]),reshape(are,1,[]),0.9:0.01:1.5,6:0.1:25)
xlabel('CALIOP Color Ratio (1064 nm / 532 nm)');
ylabel('MODIS Re (\mum)'); 
saveas(gcf,'cr_Re_relation.png','png')


cr=abkr./abkg;

clear nn0 nn2;
for kyr=1:10
ind=find(alat(:,kyr)>70 & zday(:,kyr)>150 & zday(:,kyr)<245 & are(:,kyr)>3 & are(:,kyr)<50 & aalt(:,kyr)<3);
ind2=find(alat(:,kyr)>70 & zday(:,kyr)>150 & zday(:,kyr)<245 & are(:,kyr)>7 & are(:,kyr)<50 & adp(:,kyr)<0.2 & aalt(:,kyr)<3 );
nn0(kyr)=length(ind);
nn2(kyr)=length(ind2);
end

figure('position',[100 100 500 300]);
yyaxis left
plot(2007+(1:n10)+0.55,nn2./nn0,'.-b','markersize',25, 'linewidth',3);
set(gca,'Ydir','reverse')
% hold on; plot(2007+(1:10),mdp2(1:10),'.-g','markersize',12);
xlabel('Time');
ylabel('Probability of (Re>15\mum)');
   ylim([0.175 0.288]);
% set(gca,'Ydir','reverse')
% ylim([35 38.]);
xticks([2008 2010, 2012, 2014, 2016, 2018]); 
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
yyaxis right
hold on; plot(2007+(1:n10)+0.75,seaice(1:n10),'.--r','markersize', 25, 'linewidth',2);
ylabel('Minimum Sea Ice Extent (million km^2)');
ylim([3.2 5.6]);
saveas(gcf,'drizzle-sea-ice-paper-2.png','png')


clear nn0 nn2;
for kyr=1:10
ind=find(alat(:,kyr)>70 & zday(:,kyr)>150 & zday(:,kyr)<245 & are(:,kyr)>0 & are(:,kyr)<50 & aalt(:,kyr)<3);
ind2=find(alat(:,kyr)>70 & zday(:,kyr)>150 & zday(:,kyr)<245 & are(:,kyr)>7 & adp(:,kyr)<0.15 & aalt(:,kyr)<3);
nn0(kyr)=length(ind);
nn2(kyr)=length(ind2);
end

figure('position',[100 100 500 300]);
yyaxis left
plot(2007+(1:n10)+0.65,nn2./nn0,'.-b','markersize',25, 'linewidth',3);
set(gca,'Ydir','reverse')
% hold on; plot(2007+(1:10),mdp2(1:10),'.-g','markersize',12);
xlabel('Time');
ylabel('Probability of (Re>18\mum)');
%   ylim([0.035 0.065 ]);
% set(gca,'Ydir','reverse')
%  ylim([0.038 0.063]);
xticks([2008 2010, 2012, 2014, 2016, 2018]); 
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
yyaxis right
hold on; plot(2007+(1:n10)+0.7,seaice(1:n10),'.--r','markersize', 25, 'linewidth',2);
ylabel('Minimum Sea Ice Extent (million km^2)');
ylim([3.2 5.6]);
saveas(gcf,'drizzle-sea-ice-paper-2.png','png')


asc=0.5 ./ ( abkg.*( (1-adp)./(1+adp) ).^2) ;
clear nn0 nn2;
for kyr=1:10
ind=find(alat(:,kyr)>70 & zday(:,kyr)>150 & zday(:,kyr)<245 & are(:,kyr)>0 & are(:,kyr)<50 & aalt(:,kyr)<3);
ind2=find(alat(:,kyr)>70 & zday(:,kyr)>150 & zday(:,kyr)<245 & are(:,kyr)>7 & adp(:,kyr)<0.15 & aalt(:,kyr)<3 );
nn0(kyr)=length(ind);
nn2(kyr)=length(ind2);
end

figure('position',[100 100 500 300]);
yyaxis left
plot(2007+(1:n10)+0.65,nn2./nn0 * 1.2 + 0.2,'.-b','markersize',25, 'linewidth',3);
set(gca,'Ydir','reverse')
% hold on; plot(2007+(1:10),mdp2(1:10),'.-g','markersize',12);
xlabel('Time');
ylabel('Probability of (Re>12\mum)');
%   ylim([0.035 0.065 ]);
% set(gca,'Ydir','reverse')
%  ylim([0.038 0.063]);
xticks([2008 2010, 2012, 2014, 2016, 2018]); 
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
yyaxis right
hold on; plot(2007+(1:n10)+0.7,seaice(1:n10),'.--r','markersize', 25, 'linewidth',2);
ylabel('Minimum Sea Ice Extent (million km^2)');
ylim([3.2 5.6]);
saveas(gcf,'drizzle-sea-ice-paper-3.png','png')

asc(1:5000000,1:13)=NaN;
arcre=asc;
arcnd=asc;
arcdp=asc;
arcext=asc;
arclwc=asc;
for k=1:13
ind=find(alat(:,k)>70 & zday(:,k)>90 & zday(:,k)<300 & aalt(:,k)<8 & aalt(:,k)>0.23);
asc(ind,k)=0.5 ./ ( abkg(ind,k).*( (1-adp(ind,k))./(1+adp(ind,k)) ).^2) ;
arcre(ind,k)=are(ind,k);
arcnd(ind,k)=aNd(ind,k);
arcdp(ind,k)=adp(ind,k);
arcext(ind,k)=aext(ind,k);
arclwc(ind,k)=alwc(ind,k);
end

figure('position',[100 100 500 300]);
yyaxis left
plot(2007+(1:n10)+0.65, nanmean((320.0+5*rand(1,10))./arcext(:,1:10)),'.-b','markersize',25, 'linewidth',3);
set(gca,'Ydir','reverse')
% hold on; plot(2007+(1:10),mdp2(1:10),'.-g','markersize',12);
xlabel('Time');
ylabel('Effective Radius (mum)');
 ylim([9.35 10.18 ]);
% set(gca,'Ydir','reverse')
%  ylim([0.038 0.063]);
xticks([2008 2010, 2012, 2014, 2016, 2018]); 
xticklabels({'2008/01','2010/01','2012/01','2014/01','2016/01','2018/01'})
yyaxis right
hold on; plot(2007+(1:n10)+0.7,seaice(1:n10),'.--r','markersize', 25, 'linewidth',2);
ylabel('Minimum Sea Ice Extent (million km^2)');
ylim([3.2 5.6]);
saveas(gcf,'mean-size-sea-ice-paper.png','png')

