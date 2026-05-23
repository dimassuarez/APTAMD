% Octave script for statistics of BP and NP contacts

pop_thres=0.10;

fileres='reslabel.txt';
if (  exist(fileres,"file") == 0 )
      disp(['reslabel file  ',fileres,' not found'])
      exit
end
fid=fopen(fileres);
XRES=textscan(fid,'%s\n') ;
fclose(fid);
[nres,idummy]=size(XRES{1,1}) ;

BP_data=zeros(nres,nres);
NP_data=zeros(nres,nres);
NP_STACK_data=zeros(nres,nres);
nsnap_tot=0;

filelist='filelist.txt';
if (  exist(filelist,"file") == 0 )
      disp(['file list ',filelist,' not found'])
      exit
end
fid=fopen(filelist);
XFILE=textscan(fid,'%s\n') ;
fclose(fid);
[nfile,idummy]=size(XFILE{1,1}) ;


% Loop over data files

for ifile=[1:nfile]

disp(['file=',num2str(ifile)])
bp_n_file=strcat ( XFILE{1,1}{ifile},'_bp_n.dat');
np_n_file=strcat ( XFILE{1,1}{ifile},'_np_n.dat'); 
bp_ij_file=strcat ( XFILE{1,1}{ifile},'_bp_ij.dat'); 
np_ij_file=strcat ( XFILE{1,1}{ifile},'_np_ij.dat'); 

%The following arrays contain the number of pair interactions 
% and non-pair interactions in each snapshot
bp_n=load(bp_n_file);
bp_n ( find ( bp_n == -1 ) ) = [];
np_n=load(np_n_file);
np_n ( find ( np_n == -1 ) ) = [];

nsnap_bp=length(bp_n); nsnap_np=length(np_n);
if nsnap_bp != nsnap_np 
   display('Problem while reading BP_N/NP_N arrays')
   display(bp_n_file)
   nsnap_bp
   display(np_n_file)
   nsnap_np
   exit
else
   nsnap=nsnap_bp;
end

%The following two-column array contains the I-J residue indexes
%for the base pair interactions
BP=load(bp_ij_file);
ndat_bp_ij=length(BP(:,1));

if (sum(bp_n) != ndat_bp_ij)
   display('Inconsistency in BP_N and BP_IJ data')
   ndat_bp_n=sum(bp_n)
   ndat_bp_ij
   exit 
end


% The following three-column array contains the I-J residue indexes
% and TYPE for the non-pair interactions. If TYPE=1 then it is stacking
% else it is interbase
NP=load(np_ij_file);
ndat_np_ij=length(NP(:,1));

if (sum(np_n) != ndat_np_ij)
   display('Inconsistency in NP_N and NP_IJ data')
   ndat_np_n=sum(np_n)
   ndat_np_ij
   exit 
end

% DO the counting : Smart use of vector operators
ix=BP(:,1); 
iy=BP(:,2);
indx=ix + nres*(iy-1);
[ c , ii ] = hist ( indx , [1:max(indx)] );
BP_data(ii)=BP_data(ii)+c;
clear ix iy c ii indx;

ix=NP(:,1); 
iy=NP(:,2);
iz=NP(:,3);
indx=ix + nres*(iy-1);
jndx= indx .* iz ;
[ c , ii ] = hist ( indx  ( find ( jndx ==  0 )  ) , [1:max(indx)] );
NP_data(ii)=NP_data(ii)+c;
clear c ii ;
[ c , ii ] = hist ( indx  ( find ( jndx >  0 )  ), [1:max(indx)] );
NP_STACK_data(ii)=NP_STACK_data(ii)+c;
clear ix iy iz  c ii indx jndx;

nsnap_tot=nsnap_tot + nsnap 
end

% Getting fractions 
BP_data=(BP_data+BP_data')/nsnap_tot;
NP_data=(NP_data+NP_data')/nsnap_tot;
NP_STACK_data=(NP_STACK_data+NP_STACK_data')/nsnap_tot;
NP_FULL_data=NP_data+NP_STACK_data;
 

%figure settings   

XT=[1:nres];
YT=XT;

clf();
imagesc(XT,YT, BP_data',[0 1])
colorbar();
title('Base Pairing Contacts') 
xlabel('resid')
ylabel('resid')
% Custom gridlines and axis:
M=nres;
axis([1-0.5 M-0.5 1-0.5 M-0.5]);
for i = [1:1:M]; line([i-0.5,i-0.5], [1-0.5, M-0.5] , 'linestyle','-','color','m', 'linewidth',0.25); end;
for i = [1:1:M]; line( [1-0.5, M-0.5], [i-0.5,i-0.5] , 'linestyle','-','color','m', 'linewidth',0.25); end;
print(1,'BP_stat.png','-dpng','-r600')

clf();
imagesc(XT,YT, NP_FULL_data',[0 1])
colorbar();
title('Non Pairing Contacts') 
xlabel('resid')
ylabel('resid')
% Custom gridlines and axis:
M=nres;
axis([1-0.5 M-0.5 1-0.5 M-0.5]);
for i = [1:1:M]; line([i-0.5,i-0.5], [1-0.5, M-0.5] , 'linestyle','-','color','m', 'linewidth',0.25); end;
for i = [1:1:M]; line( [1-0.5, M-0.5], [i-0.5,i-0.5] , 'linestyle','-','color','m', 'linewidth',0.25); end;
print(1,'NP_FULL_stat.png','-dpng','-r600')

clf();
imagesc(XT,YT, NP_FULL_data',[0 1])
colorbar();
title('Non Pairing (Stacking) Contacts') 
xlabel('resid')
ylabel('resid')
% Custom gridlines and axis:
M=nres;
axis([1-0.5 M-0.5 1-0.5 M-0.5]);
for i = [1:1:M]; line([i-0.5,i-0.5], [1-0.5, M-0.5] , 'linestyle','-','color','m', 'linewidth',0.25); end;
for i = [1:1:M]; line( [1-0.5, M-0.5], [i-0.5,i-0.5] , 'linestyle','-','color','m', 'linewidth',0.25); end;
print(1,'NP_STACK_stat.png','-dpng','-r600')
% Statistics


fid=fopen('BP_NP_stat.dat','w');

lxy=(triu(BP_data)>pop_thres);
[ix, iy]= find ( lxy ) ; 
[ pop,  ipop ]=sort ( BP_data(lxy)  )   ;
fprintf(fid,'=============================== \n')
fprintf(fid,'    BASE PAIRING CONTACTS \n')
fprintf(fid,'=============================== \n')
fprintf(fid,'   IRES    JRES   ABUND.(percen) \n')
fprintf(fid,'================================ \n')
for i=[length(ipop):-1:1]
    fprintf(fid,'%6s ... %6s   %8.4f \n', ...
    XRES{1,1}{ix ( ipop(i) ) }, XRES{1,1}{iy ( ipop(i) ) }, pop(i)*100.0)
end
clear lxy ix iy pop ipop ;

lxy=(triu(NP_FULL_data)>pop_thres);
[ix, iy]= find ( lxy ) ; 
[ pop,  ipop ]=sort ( NP_FULL_data(lxy)  )   ;
fprintf(fid,'=============================== \n')
fprintf(fid,'  NON PAIRING (FULL) CONTACTS \n')
fprintf(fid,'=============================== \n')
fprintf(fid,'   IRES    JRES   ABUND.(percen) \n')
fprintf(fid,'================================ \n')
for i=[length(ipop):-1:1]
    fprintf(fid,'%6s ... %6s   %8.4f \n', ...
    XRES{1,1}{ix ( ipop(i) ) }, XRES{1,1}{iy ( ipop(i) ) }, pop(i)*100.0)
end
clear lxy ix iy pop ipop ;

lxy=(triu(NP_STACK_data)>pop_thres);
[ix, iy]= find ( lxy ) ; 
[ pop,  ipop ]=sort ( NP_STACK_data(lxy)  )   ;
fprintf(fid,'=============================== \n')
fprintf(fid,'NON PAIRING (STACKING) CONTACTS \n')
fprintf(fid,'=============================== \n')
fprintf(fid,'   IRES    JRES   ABUND.(percen) \n')
fprintf(fid,'================================ \n')
for i=[length(ipop):-1:1]
    fprintf(fid,'%6s ... %6s   %8.4f \n', ...
    XRES{1,1}{ix ( ipop(i) ) }, XRES{1,1}{iy ( ipop(i) ) }, pop(i)*100.0)
end
clear lxy ix iy pop ipop ;

fclose(fid);

