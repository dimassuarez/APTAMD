% Octave script for examining RMSD-RGYR metrics 

%Reading RMS 
RMS=load("DUMMY_RMS") ;
nsnap=length(RMS(:,1));
%Reading SURF
SURF=load("DUMMY_SURF") ;
nsurf=length(SURF(:,1));
if ( nsurf  != nsnap )
   display('Inconsistency in RMSD/SURF data')
   if ( nsurf < nsnap )
          exit
   else 
      display(['Using only ',num2str(nsnap),' data from SURF'])
      SURF(nsnap+1:nsurf,:)=[];
   end
end

%Filtering -1.0 SURF values
indx=find(SURF(:,2) < 0.0 );
RMS(indx,:)=[];
SURF(indx,:)=[];

nsnap=length(RMS(:,1));

isample=[1:10:nsnap];
R2_SURF=corr(RMS(:,2),SURF(:,2))^2;
%Representacion grafica 
clf();
h=figure(1);
Ang=[char(0xC3),char(0x85)];
plot (RMS(isample,2),SURF(isample,2),'o','markersize',10) 
grid on;
xlabel(['RMSD (',Ang,')'])
ylabel(['SURF (',Ang,'**2)'])
title([' RMS vs SURF_ R^2 = ',num2str(R2_SURF)],'Fontsize',12)
set(gca,'fontsize',18)

W = 9; H = 6;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print(h,['DUMMY_PNG','_rmsd_rgyr_surf.png'],'-dpng','-color')

