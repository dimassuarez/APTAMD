% Octave script for determining INF metrics for ssDNA simulations

nres=DUMMY_NRES;

%The following arrays contain the number of pair interactions 
% and non-pair interactions in each snapshot
bp_n=load("DUMMY_BP_N");
bp_n ( find ( bp_n == -1 ) ) = [];
np_n=load("DUMMY_NP_N");
np_n ( find ( np_n == -1 ) ) = [];

nsnap_bp=length(bp_n); nsnap_np=length(np_n);
if nsnap_bp != nsnap_np 
   display('Problem while reading BP_N/NP_N arrays')
   exit
else
   nsnap=nsnap_bp;
end

%The following two-column array contains the I-J residue indexes
%for the base pair interactions
BP=load("DUMMY_BP_IJ");
ndat_bp_ij=length(BP(:,1));

if (sum(bp_n) != ndat_bp_ij)
   display('Inconsistency in BP_N and BP_IJ data')
   sum(bp_n)
   ndat_bp_ij
   exit 
end

% The following three-column array contains the I-J residue indexes
% and TYPE for the non-pair interactions. If TYPE=1 then it is stacking
% else it is interbase
NP=load("DUMMY_NP_IJ");
ndat_np_ij=length(NP(:,1));

if (sum(np_n) != ndat_np_ij)
   display('Inconsistency in NP_N and NP_IJ data')
   exit 
end

% Loading reference data
S_BP_REF=zeros(nres);
fid=fopen("DUMMY_REF");
nbp_ref=fscanf(fid,'%i','C');
for k=[1:nbp_ref]
   [i,j]=fscanf(fid,'%i%i','C');
    S_BP_REF(i,j)=1.0;
    S_BP_REF(j,i)=1.0;
end

nnptot_ref=fscanf(fid,'%i','C');
S_STACK_REF=zeros(nres);
S_NP_REF=zeros(nres);
nnp_ref=0;
nstack_ref=0;
for l=[1:nnptot_ref]
   [i,j,k]=fscanf(fid,'%i%i%i','C');
    if k == 0
       S_NP_REF(i,j)=1.0;
       S_NP_REF(j,i)=1.0;
       nnp_ref=nnp_ref+1;
    else
       S_STACK_REF(i,j)=1.0;
       S_STACK_REF(j,i)=1.0;
       nstack_ref=nstack_ref+1;
    end 
end
fclose(fid);

% Computing INF metric index for each snapshot

fod=fopen('DUMMY_FOUT','w+');
fprintf(fod,'# Isnap, n_bp, n_np, n_stack, INF_bp_stack, INF_nb_stack, MC_bp_np_stack \n') 
fprintf(fod,'# REF , %i, %i, %i \n',nbp_ref,nnp_ref,nstack_ref) 

k0_bp=0;
k0_np=0;
for isnap=[1:nsnap]

  if ( bp_n(isnap) == -1 ) && ( np_n(isnap) == -1 )

    nnp=0;nstack=0;
    INF_bp_stack=INF_np_stack=INF_bp_np_stack=-1; 
    fprintf(fod,' %i, %i, %i, %i, %i, %i, %i \n', ...
          isnap,bp_n(isnap),nnp,nstack,INF_bp_stack,INF_np_stack,INF_bp_np_stack)

  else 

  S_BP=zeros(nres);
  S_NP=zeros(nres);
  S_STACK=zeros(nres);
  
  for k=[k0_bp+1:k0_bp+bp_n(isnap)] 
      i=BP(k,1);j=BP(k,2);
      S_BP(i,j)=1.0;
      S_BP(j,i)=1.0;
  end 
  k0_bp=k0_bp+bp_n(isnap);

  nnp=0;nstack=0;
  for k=[k0_np+1:k0_np+np_n(isnap)] 
      i=NP(k,1);j=NP(k,2);itype=NP(k,3);
      if ( itype == 0 ) 
        S_NP(i,j)=1.0;
        S_NP(j,i)=1.0;
        nnp=nnp+1;
      else
        S_STACK(i,j)=1.0;
        S_STACK(j,i)=1.0;
        nstack=nstack+1;
      end
  end 
  k0_np=k0_np+np_n(isnap);

% False negatives: interactions in REF but not in snap
% False postives: interactions in snap but not in REF
  SDIF_BP = S_BP_REF - S_BP;
  FN_bp=numel(find(SDIF_BP ==  1.0))/2;
  FP_bp=numel(find(SDIF_BP == -1.0))/2;
  SDIF_NP = S_NP_REF - S_NP;
  FN_np=numel(find(SDIF_NP ==  1.0))/2;
  FP_np=numel(find(SDIF_NP == -1.0))/2;
  SDIF_STACK = S_STACK_REF - S_STACK;
  FN_stack=numel(find(SDIF_STACK ==  1.0))/2;
  FP_stack=numel(find(SDIF_STACK == -1.0))/2;
%True positives; interactions both in REF and snap
  TP_bp=nbp_ref - FN_bp;
  TP_np=nnp_ref - FN_np;
  TP_stack=nstack_ref - FN_stack;

% MM1GC correlation coeffs

  PPV=(TP_bp+TP_stack)/((TP_bp+TP_stack)+(FP_bp+FP_stack));
  STY=(TP_bp+TP_stack)/((TP_bp+TP_stack)+(FN_bp+FN_stack));
  INF_bp_stack=sqrt(PPV*STY);

  PPV=(TP_np+TP_stack)/((TP_np+TP_stack)+(FP_np+FP_stack));
  STY=(TP_np+TP_stack)/((TP_np+TP_stack)+(FN_np+FN_stack));
  INF_np_stack=sqrt(PPV*STY);

  PPV=(TP_bp+TP_np+TP_stack)/((TP_bp+TP_np+TP_stack)+(FP_bp+FP_np+FP_stack));
  STY=(TP_bp+TP_np+TP_stack)/((TP_bp+TP_np+TP_stack)+(FN_bp+FN_np+FN_stack));
  INF_bp_np_stack=sqrt(PPV*STY);


  fprintf(fod,' %i, %i, %i, %i, %f, %f, %f \n', ...
          isnap,bp_n(isnap),nnp,nstack,INF_bp_stack,INF_np_stack,INF_bp_np_stack)

  end

  INF_mat(isnap,1)=INF_bp_stack;
  INF_mat(isnap,2)=INF_np_stack;
  INF_mat(isnap,3)=INF_bp_np_stack;
 

end
 
iplot=DUMMY_DOPLOT;
if iplot == 0
   exit
end

%Reading RMS 
RMS=load("DUMMY_RMS") ;
if (length(RMS(:,1)) != nsnap )
   display('Inconsistency in RMS data')
   if (length(RMS(:,1)) < nsnap )
          exit
   else 
      display(['Using only ',num2str(nsnap),' data from RMS'])
   end
end

Ang=[char(0xC3),char(0x85)];
for i=[1:3]
  indx=find( INF_mat(:,i) != -1 );
  R2(i)=corr(RMS(indx,2),INF_mat(indx,i))^2;
%Representacion grafica 
  clf();
  h=figure(1);
  isample=round([1:10:max(indx)]);
  plot (RMS(isample,2),INF_mat(isample,i),'o','markersize',10)
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
  print(h,['DUMMY_PNG','_rmsd_inf_',num2str(i),'.png'],'-dpng','-color')
end

