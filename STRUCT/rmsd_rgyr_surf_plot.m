% Octave script for examining RMSD-RGYR-SURF data 
0;
function [ymed, med_snap]=ymed_curve(y,ninter) 
n=length(y); 
jndx=[1:ninter:n-2*ninter];
kndx=[2*ninter:ninter:n];
ymed(1)=sum(y(1:ninter))/ninter;
nmed=length(jndx);
for i=[1:nmed]
   ymed(i+1)=sum(y([jndx(i):kndx(i)]))/(2*ninter);
end
ymed(nmed+2)=sum(y(n-ninter:n))/ninter;
med_snap=fix(linspace(1,n,nmed+2));
end

%Reading RMS 
RMSD=load("RMSD.dat") ;
RGYR=load("RGYR.dat") ;

RMSD_aver=mean(RMSD)  ;
RMSD_std=std(RMSD)  ;
RMSD_max=max(RMSD)  ;
RMSD_min=min(RMSD)  ;
RGYR_aver=mean(RGYR) ; 
RGYR_std=std(RGYR)  ;
RGYR_max=max(RGYR)  ;
RGYR_min=min(RGYR)  ;

RMSD_snap=length(RMSD) ;
RGYR_snap=length(RGYR) ;

[RMSD_med,RMSD_med_snap]=ymed_curve(RMSD,fix(RMSD_snap*0.05)+1);
[RGYR_med,RGYR_med_snap]=ymed_curve(RGYR,fix(RGYR_snap*0.05)+1);


if ( exist('SURF.dat','file') )

SURF=load("SURF.dat") ;
% Filtering out -1.0 SURF values 
indx=find(SURF < 0.0 );
SURF(indx)=[];
SURF_aver=mean(SURF);
SURF_std=std(SURF);
SURF_max=max(SURF);
SURF_min=min(SURF);
SURF_snap=length(SURF);
[SURF_med,SURF_med_snap]=ymed_curve(SURF,fix(SURF_snap*0.05)+1);
nplt=3;

else

nplt=2;

end


clf();
h=figure(1);
Ang=[char(0xC3),char(0x85)];
subplot(nplt,1,1)
plot (RMSD,'-r','linewidth',1.5,RMSD_med_snap, RMSD_med,'-y','linewidth',5.0) 
grid on;
xlabel([' # Snap '])
ylabel(['RMSD (',Ang,')'])
text(  fix(0.10*RMSD_snap),  1.1*RMSD_max ,['RMSD=',num2str(RMSD_aver,"%5.2f") ...
     ,'\pm',num2str(RMSD_std,"%5.2f"),' ',Ang ], 'fontsize',18)
title(' RMSD and R_{gyr} plots')
set(gca,'fontsize',18)
subplot(nplt,1,2)
plot (RGYR,'-b','linewidth',1.5,RGYR_med_snap,RGYR_med,'-y','linewidth',5.0) 
grid on;
xlabel([' # Snap '])
ylabel(['R_{gyr} (',Ang,')'])
text(  fix(0.10*RGYR_snap), 1*RGYR_max ,['R_{gyr}=',num2str(RGYR_aver,"%5.2f") ...
     ,'\pm',num2str(RGYR_std,"%5.2f"),' ',Ang ], 'fontsize',18)
set(gca,'fontsize',18)

if ( nplt == 3 )

subplot(nplt,1,3)
plot (SURF,'-m','linewidth',1.5,SURF_med_snap,SURF_med,'-y','linewidth',5.0) 
grid on;
xlabel([' # Snap '])
ylabel(['SURF(',Ang,'**2)'])
text(  fix(0.10*SURF_snap), 1*SURF_max ,['SURF=',num2str(SURF_aver,"%7.2f") ...
     ,'\pm',num2str(SURF_std,"%7.2f"),' ',Ang ], 'fontsize',18)
set(gca,'fontsize',18)

W = 9; H = 9;

else

W = 9; H = 6;

end


set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
if nplt == 3 
  print(h,'RMSD_RGYR_SURF.png','-dpng','-color')
else
  print(h,'RMSD_RGYR.png','-dpng','-color')
end


# Correlation plots
if RMSD_snap == RGYR_snap 

nsnap=RMSD_snap;
isample=[1:10:nsnap];
R2=corr(RMSD(1:nsnap),RGYR(1:nsnap))^2;

%Representacion grafica 
clf();
h=figure(1);
Ang=[char(0xC3),char(0x85)];
plot (RMSD(isample),RGYR(isample),'o','markersize',10) 
grid on;
xlabel(['RMSD (',Ang,')'])
ylabel(['R_{gyr} (',Ang,')'])
title([' RMS vs RGYR_ R^2 = ',num2str(R2)],'Fontsize',12)
set(gca,'fontsize',18)
W = 9; H = 6;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print(h,'RMSD_RGYR_corr.png','-dpng','-color')

end

if ( nplt == 3 ) 

RMSD(indx)=[];
RMSD_snap=length(RMSD) ;

if RMSD_snap == SURF_snap 

nsnap=RMSD_snap;
isample=[1:10:nsnap];
R2=corr(RMSD(1:nsnap),SURF(1:nsnap))^2;

%Representacion grafica 
clf();
h=figure(1);
Ang=[char(0xC3),char(0x85)];
plot (RMSD(isample),RGYR(isample),'o','markersize',10) 
grid on;
xlabel(['RMSD (',Ang,')'])
ylabel(['SURF (',Ang,'**2)'])
title([' RMS vs SURF_ R^2 = ',num2str(R2)],'Fontsize',12)
set(gca,'fontsize',18)
W = 9; H = 6;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print(h,'RMSD_SURF_corr.png','-dpng','-color')

end


end

