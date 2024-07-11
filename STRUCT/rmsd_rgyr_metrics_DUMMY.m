% Octave script for examining RMSD-RGYR metrics 

%Reading RMS 
RMS=load("DUMMY_RMS") ;
nsnap=length(RMS(:,1));
%Reading RGYR
RGYR=load("DUMMY_RGYR") ;
if (length(RGYR(:,1)) != nsnap )
   display('Inconsistency in RMSD/RGYR data')
   if (length(RGYR(:,1)) < nsnap )
          exit
   else 
      display(['Using only ',num2str(nsnap),' data from RGYR'])
   end
end

isample=[1:10:nsnap];
R2_RGYR=corr(RMS(1:nsnap,2),RGYR(1:nsnap,2))^2;
%Representacion grafica 
clf();
h=figure(1);
Ang=[char(0xC3),char(0x85)];
plot (RMS(isample,2),RGYR(isample,2),'o','markersize',10) 
grid on;
xlabel(['RMSD (',Ang,')'])
ylabel(['R_{gyr} (',Ang,')'])
title([' RMS vs RGYR_ R^2 = ',num2str(R2_RGYR)],'Fontsize',12)
set(gca,'fontsize',18)

W = 9; H = 6;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print(h,['DUMMY_PNG','_rmsd_surf.png'],'-dpng','-color')

