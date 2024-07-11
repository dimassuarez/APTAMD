% Octave script for examining RMSD-INF metrics 
0;

function ymed=ymed_curve(x,y,ninter) 
nx=length(x); 
n=length(y); 
for ix=[1:nx] 
  i=x(ix);
  if  i > ninter && i < n-ninter 
     ymed(ix)=mean(y(i-ninter:i+ninter)); 
  elseif i <= ninter 
     ymed(ix)=mean(y(1:i+ninter+(ninter-i))); 
  elseif i >= n-ninter 
   ymed(ix)=mean(y(i-ninter:n)); 
  end 
end 
end 

% Reading RMS and INF.dat 
% Pruning out INF=-1 data 

INF=load("INF.dat") ;
INF_snap=length(INF); 
indx=find(INF != -1);
INF=INF(indx);

RMSD=load("RMSD.dat") ;
RMSD_snap=length(RMSD);

if ( INF_snap != RMSD_snap ) 
  disp( ' Possible inconsistency in RMSD / INF data ')
end
RMSD=RMSD(indx);
INF_snap=length(indx); 
RMSD_snap=INF_snap;

INF_aver=mean(INF);
INF_std=std(INF);
INF_max=max(INF);
INF_min=min(INF);

x= round(linspace(1,INF_snap,min([ 100, INF_snap ])));
INF_med=ymed_curve(x,INF,fix(INF_snap*0.05));

clf();
h=figure(1);
title(' INF plot')
set(gca,'fontsize',18)
plot (INF,'ob','markersize',3,x,INF_med,'-r','linewidth',5.0) 
grid on;
xlabel([' # Snap '])
ylabel('INF')
text(  fix(0.10*INF_snap), 0.9*INF_max ,['INF=',num2str(INF_aver,"%5.2f") ...
     ,'\pm',num2str(INF_std,"%5.2f")], 'fontsize',18)
set(gca,'fontsize',18)
W = 9; H = 6;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print(h,'INF.png','-dpng','-color')

# Correlation plot 
if RMSD_snap == INF_snap 

nsnap=RMSD_snap;
isample=[1:1:nsnap];
R2=corr(RMSD(1:nsnap),INF(1:nsnap))^2;
%Representacion grafica 
clf();
h=figure(1);
Ang=[char(0xC3),char(0x85)];
plot (RMSD(isample),INF(isample),'o','markersize',10) 
grid on;
xlabel(['RMSD (',Ang,')'])
ylabel('INF')
title([' RMS vs INF_ R^2 = ',num2str(R2)],'Fontsize',12)
set(gca,'fontsize',18)
W = 9; H = 6;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print(h,'RMSD_INF_corr.png','-dpng','-color')

end
