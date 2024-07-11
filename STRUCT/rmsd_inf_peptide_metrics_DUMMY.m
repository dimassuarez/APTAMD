% Octave script for examining RMSD INF metrics for peptide simulations

%The following arrays contain the number of HB interactions 
hb=load("DUMMY_HB");
nsnap=length(hb)-1;

% reference data
nhb=length(hb(1,:))-1;
hb_ref=hb(1,2:nhb+1);
nocontact_REF=sum(hb_ref);   % A ficticious contact accounts for ZERO interactions
if nocontact_REF == 0 
   hb_ref(nhb+1)=1.0;
else
   hb_ref(nhb+1)=0.0;
end
nhb_ref=sum(hb_ref);

% Computing INF metric index for each snapshot

fod=fopen('DUMMY_FOUT','w+');

fprintf(fod,'# Isnap, n_hb, INF_hb \n') 
fprintf(fod,'# REF , %i, \n',nhb_ref) 

k0_hb=0;
for isnap=[1:nsnap]
  hb_snap=hb(isnap+1,2:nhb+1);
  nocontact=sum(hb_snap);
  if nocontact == 0 
     hb_snap(nhb+1)=1.0;   % Ficticious interaction means NO interaction
  else
     hb_snap(nhb+1)=0.0;   % Ficticious interaction means NO interaction
  end

  SDIF_HB = hb_ref - hb_snap; 
% False negatives: interactions in REF but not in snap
  FN_hb=numel(find(SDIF_HB ==  1.0)); 
% False positives: interactions in snap but not in REF
  FP_hb=numel(find(SDIF_HB == -1.0)); 
%True positives; interactions both in REF and snap
  TP_hb=nhb_ref - FN_hb ;

% MMC correlation coeffs
  PPV=(TP_hb)/((TP_hb)+(FP_hb));
  STY=(TP_hb)/((TP_hb)+(FN_hb));
  INF_hb=sqrt(PPV*STY);

  fprintf(fod,' %i, %i, %f \n', ...
          isnap,sum(hb_snap),INF_hb)

  INF_mat(isnap,1)=INF_hb;

end
fclose(fod);

iplot=DUMMY_DOPLOT;
if iplot == 0
   exit
end

%Reading RMS 
RMS=load("DUMMY_RMSD") ;
if (length(RMS(:,1)) != nsnap )
   display('Inconsistency in RMS data')
   if (length(RMS(:,1)) < nsnap )
          exit
   else 
      display(['Using only ',num2str(nsnap),' data from RMS'])
   end
end

isample=[1:10:nsnap];
R2=corr(RMS(1:nsnap,2), INF_mat(1:nsnap,1))^2;
%Representacion grafica 
clf();
h=figure(1);
Ang=[char(0xC3),char(0x85)];
plot (RMS(isample,2),INF_mat(isample,1),'o','markersize',10)
grid on;
xlabel(['RMSD (',Ang,')'])
ylabel(['INF'])
ylim([0 1])
title([' RMS vs INF R^2 = ',num2str(R2)],'Fontsize',12)
set(gca,'fontsize',18)
W = 9; H = 6;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print(h,['DUMMY_PNG','_rmsd_inf.png'],'-dpng','-color')
