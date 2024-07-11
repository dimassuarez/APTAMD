% Processing GAMD energy reweight data

disp('Processing Energy Reweigth Data')
% A guess of the maximum energy difference (probably set by inspection)
ETHRES=DUMMY_EMAX;

RC1='DUMMY_RC_1';
RC2='DUMMY_RC_2';

%Second order cumulant expansion data 
A=load('CE2.dat');
indx=find ( A(:,3) < ETHRES ); 
if RC1 == 'INF'
  f1=10.0;
else
  f1=1.0;
end
X=A(indx,1)/f1;
Y=A(indx,2);E=A(indx,3);
ndat=length(X);

ETHRES=round(max(E));

%  Max/min reaction coordinate values
clear A
A=dlmread('DUMMY_COORD.max_min');

xmax=A(1,2)/f1;
xmin=A(2,2)/f1;

fine_grid=DUMMY_GRID;

if fine_grid == 1 
  xgrid=A(3,2)/f1;   % Fine grid
else
  xgrid=A(4,2)/f1;   % Reg grid
end

ymax=A(1,3);
ymin=A(2,3);
if fine_grid == 1 
  ygrid=A(3,3);   % Fine grid
else
  ygrid=A(4,3);   % Reg grid
end

nx=(xmax-xmin)/xgrid+1;
ny=(ymax-ymin)/ygrid+1;
nx=round(nx);ny=round(ny);

% Determining Relative Population of each 2D bin 
POP=zeros(nx,ny);

ix=round (  ((X.-xmin)/xgrid).+1) ;
iy=round (  ((Y.-ymin)/ygrid).+1) ;

% Reading bulk data Reaction coordinate 1 (RC1) and Reaction coordinate (RC2) 

RC=load('DUMMY_COORD.dat');
n=length(RC);

disp('RC1/RC2 data loaded')
IRC(:,1)=round( (RC(:,1)/f1.-xmin)/xgrid.+1 );
IRC(:,2)=round( (RC(:,2).-ymin)/ygrid.+1 );
disp('Building 2D Population map')
for i=1:n
  POP( IRC(i,1) , IRC(i,2) ) =  POP( IRC(i,1) , IRC(i,2) ) + 1;
end

fid=fopen('2D_RW.dat','w');

disp('Population calculated')
disp('# RC1, RC2,  E, Pop')
fprintf(fid,'# RC1, RC2, E, Pop \n');
for i=[1:ndat]
   P(i)=POP(ix(i),iy(i))/n;
   printf(' %f, %f, %f, %15.12f, %i, %i \n',X(i),Y(i),E(i),P(i),ix(i),iy(i))
   fprintf(fid,' %f, %f, %f, %15.12f, %i, %i \n',X(i),Y(i),E(i),P(i),ix(i),iy(i))
end

% POP-Weigthed  E map 
D=-ones(nx,ny);
DP=zeros(nx,ny);

for i=[1:ndat]
    D(ix(i),iy(i))=E(i);
    DP(ix(i),iy(i))=P(i);
end

% Plotting 2D Pop Map
disp('Plotting 2D population map')
Psup=-1.0;
Pmin=-6.0;
for jx=[1:nx]
for jy=[1:ny]
   if D(jx,jy) < 0  || DP(jx,jy) < 1E-7
      DP(jx,jy)=Pmin;
   else 
     DP(jx,jy)=log10(DP(jx,jy));
   end
end
end

for jx=[1:nx]
    XT(jx)=xmin+xgrid*(jx-1);
end
for jy=[1:ny]
    YT(jy)=ymin+ygrid*(jy-1);
end

%Color map (light grey means empty)
clear A;
A=colormap('cubehelix');
% A(64,:)=0.66;  
A(64,:)=1.00;
% Uncomment
colormap(flipud(A));

%figure settings   LATTICE PLOT
clf();
h=figure(1);
FN = findall(h,'-property','FontName');
set(FN,'FontName','/usr/share/fonts/dejavu/DejaVuSerifCondensed.ttf');
FS = findall(h,'-property','FontSize');
set(FS,'FontSize',12);
set (FS,'defaultaxesfontsize', 14) 
W = 5; H = 3.75;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])

imagesc(XT,YT,DP',[Pmin Psup])

colorbar();

if RC1 == 'INF' 
  set(gca,'xtick',round(linspace(0,xmax,5)*100)/100)
end
%set(gca,'ytick',[]);
xlabel(RC1)
ylabel(RC2)

% Custom gridlines and axis:
%axis([1-0.5 M-0.5 1-0.5 M-0.5]);
%for i = [1:1:M]; line([i-0.5,i-0.5], [1-0.5, M-0.5] , 'linestyle','-','color','g', 'linewidth',0.5); end;
%for i = [1:1:M]; line( [1-0.5, M-0.5], [i-0.5,i-0.5] , 'linestyle','-','color','g', 'linewidth',0.5); end;

print(1,'CE2_pop.png','-dpng')


% Plotting 2D Energy Map
disp('Plotting Free Energy 2D map');
Esup=0.0;
Emin=1000.0;

for i=[1:ndat]
    Esup=max(E(i),Esup);
    Emin=min(E(i),Emin);
end
D=D.-Emin;
Esup=Esup-Emin;

for jx=[1:nx]
for jy=[1:ny]
   if D(jx,jy) < 0  
      D(jx,jy)=Esup;
   end
end
end

for jx=[1:nx]
    XT(jx)=xmin+xgrid*(jx-1);
end
for jy=[1:ny]
    YT(jy)=ymin+ygrid*(jy-1);
end

%Color map (light grey means empty)
clear A;
A=colormap('hot');
% A(64,:)=0.66;  
A(64,:)=1.00;
% Uncomment
colormap(A);

%figure settings   LATTICE PLOT
clf();
h=figure(1);
FN = findall(h,'-property','FontName');
set(FN,'FontName','/usr/share/fonts/dejavu/DejaVuSerifCondensed.ttf');
FS = findall(h,'-property','FontSize');
set(FS,'FontSize',12);
set (FS,'defaultaxesfontsize', 14) 
W = 5; H = 3.75;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])

imagesc(XT,YT,D',[0 Esup])

colorbar();


if RC1 == 'INF'
   set(gca,'xtick',round(linspace(xmin,xmax,5)*100)/100)
end
%set(gca,'ytick',[]);
xlabel(RC1)
ylabel(RC2)

% Custom gridlines and axis:
%axis([1-0.5 M-0.5 1-0.5 M-0.5]);
%for i = [1:1:M]; line([i-0.5,i-0.5], [1-0.5, M-0.5] , 'linestyle','-','color','g', 'linewidth',0.5); end;
%for i = [1:1:M]; line( [1-0.5, M-0.5], [i-0.5,i-0.5] , 'linestyle','-','color','g', 'linewidth',0.5); end;

print(1,'CE2.png','-dpng')

%
%  Determining min position
%
disp('Selecting min(E) snapshots') 

[ Esorted, isorted ] = sort (E ); 

disp('Loading ID and GAMD data ') 
ID=load("ID.dat");
GAMD=load("gamd.log");
A=load("SNAP.dat");
nsnap=A(:,2); nstep=A(:,4)./A(:,3);
ncrd=length(nsnap);
clear A ;
printf('Average, STD,  Max and Min  Boost potential: \n')
fprintf(fid,'Average, STD,  Max and Min  Boost potential: \n')
gaver=mean(GAMD(:,7)+GAMD(:,8));
gstd=std(GAMD(:,7)+GAMD(:,8));
gmax=max(GAMD(:,7)+GAMD(:,8));
gmin=min(GAMD(:,7)+GAMD(:,8));
printf("  %f  , %f , %f, %f \n",gaver, gstd, gmax, gmin)
fprintf(fid,"  %f  , %f , %f, %f \n",gaver, gstd, gmax, gmin)

% We checked the five lowest E bins with enough population 
l_selec=0;
clear jx jy;
for l=[ 1 : ndat ]

    k=isorted(l);
    pop_l=POP(ix(k),iy(k))/n;
    E_l=Esorted(l) 
    disp(['Processing  snapshots in bin ',num2str(l)])
    disp ([' k, l =',num2str(k),' ',num2str(l)])
    disp ([' ix, iy =',num2str(ix(k)),' ',num2str(iy(k))])
    if ( pop_l  >  0.01 ) 
 
       jx=find(IRC(:,1) == ix(k) );
       jy=find(IRC(:,2) == iy(k) );
 
       ja=intersect(jx,jy);
       nj=length(ja);
       snaps(1:nj,1)=ja;
       snaps(1:nj,2:3)=RC(ja,1:2);
       snaps(1:nj,4)=ID(ja);
       snaps(1:nj,5)=GAMD(ja,7).+GAMD(ja,8);

       format long
       m=length(snaps);
     
       [x0,ix0]=sort(snaps(:,5));
       printf(" Picking up snapshots bin= %i for  E= %f POP= %5.3f  X= %5.2f  Y= %5.2f \n", ...
                  l, E_l, pop_l, X(k), Y(k) )
       fprintf(fid," Picking up snapshots bin= %i for  E= %f POP= %5.3f  X= %5.2f  Y= %5.2f \n", ...
                  l, E_l, pop_l, X(k), Y(k) )
     
       for jj=[1:5]  
           i=ix0(jj); 
           printf("   Snapshot %i   isnap= %i   ID=%i  RC1= %f RC2= %f Boost= %f \n ", ...
           jj, round(snaps(i,1)), round(snaps(i,4)), snaps(i,2)/f1, snaps(i,3), snaps(i,5) )
           fprintf(fid,"   Snapshot %i   isnap= %i   ID=%i  RC1= %f RC2= %f Boost= %f \n ", ...
           jj, round(snaps(i,1)), round(snaps(i,4)), snaps(i,2)/f1, snaps(i,3), snaps(i,5) )
           isel=round(snaps(i,1));
           nsnap_acc=0;
           for icrd=[1:ncrd]
              nsnap_acc=nsnap_acc+nsnap(icrd);
              if ( nsnap_acc > isel ) 
                 isnap_in_crd=isel-(nsnap_acc-nsnap(icrd));
                 break
              end
           end
           printf("  bin= %i isnap_in_crd= %i  CRD_FILE= %i \n", l,(isnap_in_crd-1)*nstep(icrd)+1, icrd)
           fprintf(fid,"  bin= %i  isnap_in_crd= %i  CRD_FILE= %i \n", l, (isnap_in_crd-1)*nstep(icrd)+1, icrd)
        end
        clear jx jy ja snaps x0, ix0;
        l_selec=l_selec+1;
     else
        printf(" Bin  %i  with E =  %f and POP= %f not considered \n", l,E_l,pop_l)
        fprintf(fid," Bin  %i  with E =  %f and POP= %f not considered \n", l,E_l,pop_l)
     end
     if l_selec == 5 
       break
     end
end

exit

