# Python script for TSNE clustering
#
# Clustering Heterogeneous Conformational Ensembles of
# Intrinsically Disordered Proteins with t-Distributed Stochastic
# Neighbor Embedding
#
# J. Chem. Theory Comput. 2023, 19, 4711-4727
#
# https://doi.org/10.1021/acs.jctc.3c00224
#
# This script is an adaptation of that included in the JCTC paper.
# Some settings are controlled by the do_cluster.sh script
#

import numpy as np
from os.path import join, isfile 
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn import metrics
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Inputs for TSNE: provide distance matrix between all pairs of conformations. 

# Pairwise RMSD obtained from Gromacs in binary vector format.
# data = np.fromfile('rmsd.dat', np.float32)
# rmsd = data.reshape(int(np.sqrt(len(data))), int(np.sqrt(len(data)))) ### To reshape the vector to matrix
inputfile='DUMMY_FILE'
# Perplexity range
# perp0=DMMY_PERP0;perp1=DMMY_PERP1;perpx=DMMY_PERPX;
perp0=DUMMY_PERP0;perp1=DUMMY_PERP1;perpx=DUMMY_PERPX;
perplexityVals = range(perp0, perp1, perpx)
# Kclus number range
# kclus0=DMMY_KCLUS0;kclus1=DMMY_KCLUS1;kclusx=DMMY_KCLUSX
kclus0=DUMMY_KCLUS0;kclus1=DUMMY_KCLUS1;kclusx=DUMMY_KCLUSX
range_n_clusters = range ( kclus0, kclus1, kclusx) 
# Gradient method in tSNE (exact or barnes_hut)
gradient_method='barnes_hut'
lx=ly=0  # Plot width 
bestP=0  # 0 --> Select optimal
bestK=0  # 0 --> Select optimal
silhouette_type=2   # 2 -> SD  3-> HD 4-> SD*HD

do_tsne=1  #  If do_tsne == 0 only plot and analysis 

if ( do_tsne ) :
    # Pairwise RMSD obtained from Cpptraj (rmsd2 command) in standard data type (ASCII) 
    if ( inputfile == 'RMS2D.dat' ):
    	print('Using RMSD2.dat standard data file')
    	data=np.genfromtxt('../RMS2D.dat')
    	rmsd=np.delete(data,0,1)
    	del data
    
    else :
    # Cpptrajpairdist generated by the Cpptraj ckuster command, 
    # Note tha this matrix contains only upper non-diagonal elements
    #
    	print('Using CpptrajPairDist Netcdf file')
    	from netCDF4 import Dataset
    	nc_f = ('../CpptrajPairDist')  # Finename 
    	nc_fid = Dataset(nc_f, 'r')  # Opening 
    	print ("NetCDF dimension information:")
    	nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    	for dim in nc_dims :
    		print ("\tName:", dim )
    		print ("\t\tsize:", len(nc_fid.dimensions[dim]))
    	n = len(nc_fid.dimensions['n_rows']) 
    	x = nc_fid.variables['matrix'][:]  
    	del nc_fid
    	rmsd = np.zeros( (n,n) )
    	iupper =  np.triu_indices( n , 1 )
    	rmsd  [ iupper ] = x  
    	del x
    	rmsd = rmsd + np.transpose(rmsd)    
    
    print ( 'RMSD shape') 
    (nsnap,msnap)=rmsd.shape
    print ( rmsd.shape) 
    
    
    ##To make sure if the matrix is symmetric
    import sklearn.utils.validation as suv
    suv.check_symmetric(rmsd, raise_exception=True)
    
    print('symmetry check completed')
    
    # PCA analysis
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    pca.fit(rmsd)
    print('PCA: explained variance')
    print(pca.explained_variance_ratio_)
    print('PCA: explained variance')
    print(pca.explained_variance_ratio_)
    initial_values_for_tsne=np.transpose(pca.components_)
    del pca
    
    
    # Creating the TSNE object and projection
    
    for perp in perplexityVals:
        tsne_file="tsnep{0}".format(perp)
        seed_tsne_file="seed_tsnep{0}".format(perp)
        if ( perp  > nsnap ) :
            break 
        if ( isfile ( tsne_file ) ) :
            print(" File %s exists. Assuming a previous calc has been done"%tsne_file)
            print(" and thereby we skip perplexity value = %f"%perp)
            tsne=np.loadtxt(tsne_file)
        else :
            if ( isfile ( seed_tsne_file ) ) :
                 print(" Seed file %s exists. Reading data as initial tSNE coords"%seed_tsne_file)
                 initial_values=np.loadtxt(seed_tsne_file)

            print(" Projecting RMSD data on 2D TSNE perp=%f"%perp)
    ### metric is precomputed RMSD distance. if you provide Raw coordinates, 
    #   the TSNE will compute the distance by default with Euclidean metrics
            tsneObject = TSNE(n_components=2, perplexity=perp, early_exaggeration=10.0, learning_rate=100.0,  
            n_iter=3500, n_iter_without_progress=300, min_grad_norm=1e-7, metric="precomputed",  
            init=initial_values_for_tsne, method=gradient_method, angle=0.5, verbose=1) 
            tsne = tsneObject.fit_transform(rmsd)
            np.savetxt(tsne_file, tsne)

        for n_clusters in range_n_clusters:
            kfile_cen='kmeans_'+str(n_clusters)+'clusters_centers_tsnep'+str(perp)
            kfile_tsne='kmeans_'+str(n_clusters)+'clusters_tsnep'+str(perp)+'.dat'
            if ( isfile(kfile_cen)) and  (isfile(kfile_tsne) ):
                 print(" %s and %s already exist!"%(kfile_cen,kfile_tsne))
            else:
                 print("Extracting K-clusters n=%i"%n_clusters)
                 kmeans = KMeans(n_clusters=n_clusters).fit(tsne)
                 np.savetxt(kfile_cen, kmeans.cluster_centers_, fmt='%1.3f')
                 np.savetxt(kfile_tsne, kmeans.labels_, fmt='%1.1d')
    #### Compute silhouette score based on low-dim and high-dim distances        
                 silhouette_ld = silhouette_score(tsne, kmeans.labels_)
                 silhouette_hd = metrics.silhouette_score(rmsd, kmeans.labels_, metric='precomputed')
                 with open('silhouette.txt', 'a') as f:
                      f.write("\n")
                      print(perp, n_clusters, silhouette_ld, silhouette_hd, silhouette_ld*silhouette_hd, file =f)

##### plotting for the best cluster with highest silhouette score######       
s = np.loadtxt('silhouette.txt')

[bestP_tmp,bestK_tmp] = s[np.argmax(s[:,silhouette_type]), 0], s[np.argmax(s[:,silhouette_type]), 1]
if bestP == 0 :
   bestP= bestP_tmp
if bestK == 0 :
   bestK= bestK_tmp
besttsne = np.loadtxt('tsnep'+str(int(bestP)))
bestclust = np.loadtxt('kmeans_'+str(int(bestK))+'clusters_tsnep'+str(int(bestP))+'.dat')

# Abundances 
c_members = {i: np.where(bestclust == i)[0] for i in range(int(bestK))}
nframes_tot=0
nclus=len (c_members)
for i in range ( nclus ) :
     nframes_tot = nframes_tot+ len( c_members[i] )
abun=np.zeros( ( nclus) )
for i in range ( nclus )  :
        nframes = len( c_members[i] )
        abun[i] =  ( nframes / nframes_tot ) * 100

# Loading centers and generating distances
center = np.loadtxt('kmeans_'+str(int(bestK))+'clusters_centers_tsnep'+str(int(bestP))) #kmeans.cluster_centers_
from scipy import spatial
distance,index_dist = spatial.KDTree(besttsne, leafsize=10).query(center, 10)

# Reordering 
# XY coord of clusters
xtsne=np.zeros(int(bestK))
ytsne=np.zeros(int(bestK))
rtsne=np.zeros(int(bestK))
atsne=np.zeros(int(bestK))
qtsne=np.zeros(int(bestK))
q_abun=np.zeros(4)

# getting averaheX/Y and polar coord of each cluster
for i in range(int(bestK)):
    i_members=np.array ( c_members[i] )
    i_members=i_members.flatten()
    xtsne[i]=np.mean( besttsne[ i_members ,0] )
    ytsne[i]=np.mean( besttsne[ i_members ,1] )
    rtsne[i]=np.sqrt ( xtsne[i]**2 + ytsne[i]**2 ) 
    atsne[i]=(180.0/np.pi)*np.arccos(xtsne[i]/rtsne[i])
    if ( ytsne[i] < 0. ):
         atsne[i]=360.-atsne[i]
    if ( atsne[i] > 0.0 ) and ( atsne[i] <= 90.0) :
         qtsne[i]=0
         q_abun[0]=q_abun[0]+abun[i]
    elif ( atsne[i] > 90.0 ) and ( atsne[i] <= 180.0) :
         qtsne[i]=1
         q_abun[1]=q_abun[1]+abun[i]
    elif ( atsne[i] > 180.0 ) and ( atsne[i] <= 270.0) :
         qtsne[i]=2
         q_abun[2]=q_abun[2]+abun[i]
    else: 
         qtsne[i]=3
         q_abun[3]=q_abun[3]+abun[i]
#   print(i,xtsne[i],ytsne[i],rtsne[i],atsne[i],qtsne[i])

# Which quadrant is most populated ?
iqmax = np.argmax ( q_abun ) 
j=0
jndx=np.array([0,0,0,0],dtype=np.int8)
for i in np.flip ( np.argsort(q_abun) ) :
      jndx[i]=j
      j=j+1
q_members = { jndx[i] : np.where(qtsne == i) for i in np.flip ( np.argsort(q_abun) ) }
# Classify clusters in terms of tSNE quadrant
q_members = {i: np.where(qtsne == i) for i in np.argsort(q_abun) }
indx=np.array([],dtype=np.int8) 
nq0=0
for i in range(4):
      q_members_unsorted=np.array ( q_members[i] )
      q_members_unsorted=q_members_unsorted.flatten()
      nq_clus=len( q_members_unsorted )
      if nq_clus >  0 :
            if ( i == 0 ) :
                 rref0=np.mean( rtsne[  q_members_unsorted ] )
            rdist= ( rtsne[  q_members_unsorted ]  - rref0 )**2  
            indx_rdist=np.argsort(rdist)
            q_members_sorted=q_members_unsorted[indx_rdist]
            indx=np.append( indx, q_members_sorted )
            iamax= np.argmax ( atsne [ q_members_unsorted ] )
            rref0= rtsne [ q_members_unsorted [iamax] ]

# Rotating coord so that most populated quadrant is 0 
if ( iqmax == 0 ):
    krot=0
elif  (iqmax == 1 ) :
    krot=3
elif  (iqmax == 2 ) :
    krot=2
elif  (iqmax == 3 ) :
    krot=1

if ( iqmax > 0 ) :
     print('Rotating coords iqmax,krot',iqmax,krot)
     theta = (np.pi/2)*(-krot)
     c,s = np.cos(theta), np.sin(theta)
     R=np.array(((c,-s),(s,c)))
     besttsne = np.matmul (  besttsne , R) 
     xtsne = np.matmul ( np.transpose (np.array ([xtsne,ytsne] ) ), R)[:,0]
     ytsne = np.matmul ( np.transpose (np.array ([xtsne,ytsne] ) ), R)[:,1]

# PLaying the trick: renumbering clusters:
bestclust_sorted=np.zeros_like(bestclust)
for i in range ( len ( bestclust ) ) :
    for j in range ( len (indx ) ):
         if ( indx[j] == bestclust[i] ):
               bestclust_sorted[i]=j
               break


# Plotting tSNE representation

# Initialize the plot
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "Times New Roman"
plt.rcParams["font.size"]=14
plt.rcParams["axes.linewidth"]=3
plt.rc('axes', linewidth=1.5)
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=1)
if ( lx == 0 ) :
    lx = np.max ( np.absolute ( besttsne[:,0] ) )
if ( ly == 0 ) :
    ly = np.max ( np.absolute ( besttsne[:,1] ) )
plt.xlim(-lx, lx)
plt.ylim(-ly,ly)
plt.plot([-lx,lx],[0,0],color='k',linestyle='--',linewidth=1.5)
plt.plot([0,0],[-ly,ly],color='k',linestyle='--',linewidth=1.5)

cmap = cm.get_cmap('jet', bestK)
plt.scatter(besttsne[:,0], besttsne[:,1], c= bestclust_sorted.astype(float), s=50, alpha=0.5, cmap=cmap)

##### Adding center tags
for i in range ( len(center)) :
    xc=center[i][0]
    yc=center[i][1]
    for j in range ( len (indx ) ):
         if ( indx[j] == i ):
               ii=j
               break
    if  ( iqmax > 0 ) :
        [ xc, yc] = np.matmul ( [ xc, yc], R) 
    if  ( ii > len(center)/3 ) :
          fcolor='black'
    else:
          fcolor='yellow'
    plt.text(xc,yc,str(ii),fontsize=14,fontweight='bold',color=fcolor)

plt.gca().tick_params(width=2)

# Save figure 
plt.savefig('plot_tsnep'+str(int(bestP))+'_kmeans'+str(int(bestK))+'.png', dpi=600)    


# Printing summary file with cluster member information abundance
def ismember(B, A):
    x=-1
    for i in range ( len(A)):
        if B == A[i] :
            x=i
            return x
    return x

with open('summary_optimal_TSNE.txt', 'w') as f:
    print('# Optimal tSNE-K-means BestP=%i  BestK=%i '%(int(bestP),int(bestK)),file=f)
    print('# Total number of clusters  ',file=f)
    print('# Frame-ID Cluster-ID  ', file=f)
   
    for i in np.flip(np.argsort( abun )) :
        for k in range ( len (indx ) ):
            if ( indx[k] == i ):
                ii=k
                break
        nframes = len( c_members[i] )
        print('# Cluster %i Nframes %i  Abundance %f t-SNE x/y %f %f '%(ii, nframes, abun[i],xtsne[i],ytsne[i]), file=f)
        for j in range ( len( c_members[i] )) :
            dismin=np.min(distance[i])
            icheck = ismember( c_members[i][j] , index_dist[i] ) 
            if  ( icheck > -1 ) :
                if distance[i][icheck] == dismin :
                    print('   %i   C%s  BEST Dist= %f'%(c_members[i][j]+1,str(ii),distance[i][icheck]),file=f)
                else:
                    print('   %i   C%s  Nearest Dist= %f'%(c_members[i][j]+1,str(ii),distance[i][icheck]),file=f)
            else :
                print('   %i    C%s   '%(c_members[i][j]+1,str(ii)),file=f)
    
# Plottinh silhoutte maps 
A = np.loadtxt('silhouette.txt')
X = A[:, 0]
Y = A[:, 1]
P = np.unique(X)
K = np.unique(Y)

# Get the colormap
cmap = cm.get_cmap('viridis', 256)  
for icol in [2, 3, 4]:  

    S = A[:, icol]
    smin = np.min(S)
    smax = np.max(S)
    is_ = np.round((255 * (S - smin) / (smax - smin))).astype(int)
    ns = len(is_)
    
    ms = [36*36 for _ in range(ns)]
    Scolor = [cmap(is_[i] / 255.0) for i in range(ns)]
    
    ntick = 10
    stick = np.linspace(smin, smax, ntick)

    # Figure settings
    plt.figure(figsize=(9, 6))
    sc = plt.scatter(X, Y, s=ms, c=S, cmap=cmap, marker='o', edgecolor='k')
    plt.grid(True, which='both')
    plt.xlabel('Perplexity', fontweight='bold')
    plt.ylabel('# K clust', fontweight='bold')
    plt.xlim(0, max(P)*1.1)
    plt.ylim(0, max(K)*1.1)

    cb = plt.colorbar(sc)
    cb.set_ticks(np.linspace(smin, smax, ntick))
    cb.set_ticklabels([f'{stick[i]:.2f}' for i in range(ntick)])

    plt.savefig(f'sil_{icol}.png', dpi=300)
    plt.close()

