library(bio3d)
library(ncdf4)



mydcdfile <- "/Users/jansiess/Desktop/main_repository/v/vik_lab/projects/ferredoxin_experiments/mutter_simulations/ann/ann_full/third_run/ann_full.wa.dyn.dcd"
dcd <- read.dcd(mydcdfile)

mypdb <- "/Users/jansiess/Desktop/main_repository/v/vik_lab/projects/ferredoxin_experiments/mutter_simulations/ann/ann_full/third_run/ann_full.wa.pdb"
pdb <- read.pdb(mypdb)

  
ca.inds <- atom.select(pdb, elety="CA")

xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd,
               fixed.inds=ca.inds$xyz,
               mobile.inds=ca.inds$xyz)

rd <- rmsd(xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])
png('ann_full.thirdRun_RMSD.png', pointsize=10, width=3000, height=2300, res=500)
plot(rd, typ="l", ylab="RMSD", xlab="Frame No.", main="ANN (third run)")
points(lowess(rd), typ="l", col="red", lty=2, lwd=2)
dev.off()



# xlab="\u5c"
# xlab="\u212b"
# xlab="Å"



# RMSD Histogram
png('ann_full.thirdRun_RMSD-histo.png', pointsize=10, width=3000, height=2300, res=500)
hist(rd, breaks=40, freq=FALSE, main="RMSD Histogram - ANN (third run)", xlab="RMSD")
lines(density(rd), col="gray", lwd=3)
dev.off()

# RMSF
png('ann_full.thirdRun_RMSF.png', pointsize=10, width=3000, height=2300, res=500)
rf <- rmsf(xyz[,ca.inds$xyz])
plot(rf, ylab="RMSF", xlab="Residue Position", typ="l", main="RMSF - ANN (third run)")
dev.off()


#######
# PCA # 
pc <- pca.xyz(xyz[,ca.inds$xyz])
png('ann_full.thirdRun_PC.png', pointsize=10, width=3000, height=2300, res=500)
plot(pc,, col=bwr.colors(nrow(xyz)) )
dev.off()

# Quick clustering in PC-space to further highlight these distinct conformers
hc <- hclust(dist(pc$z[,1:2]))
grps <- cutree(hc, k=2)
png('ann_full.thirdRun_PC_cluster.png', pointsize=10, width=3000, height=2300, res=500)
plot(pc, col=grps)
dev.off()

# Examine the contribution of each residue to the first two principal components
png('ann_full.thirdRun_pc1_pc2_overlay.png', pointsize=10, width=3000, height=2300, res=500)
plot.bio3d(pc$au[,1], ylab="PC1 (A)", xlab="Residue Position", typ="l",
           main="PC1 PC2 Residue Fluctuation Overlay")
points(pc$au[,2], typ="l", col="blue")
dev.off()

# Interpolating between the most dissimilar structures in the distribution along
# a given principal component
p1 <- mktrj.pca(pc, pc=1, b=pc$au[,1], file="pc1.pdb")
p2 <- mktrj.pca(pc, pc=2, b=pc$au[,2], file="pc2.pdb")



# Normal Mode Calculation
modes <- nma(pdb)
print(modes)
plot(modes, sse=pdb)

# Make a PDB trajectory
mktrj.nma(modes, mode=7)
mktrj.nma(modes, mode=8)




# Writing these trajectory's as AMBER NetCDF
# The below code writes out a visualization of the first principal component of
# the HIV protease. Color scale from blue to red depicts low to high atomic
# displacements.
write.ncdf(p1, "trj_pc1.nc")


##############################
# Cross-Correlation Analysis # 

# Creating a dynamical cross-correlation map (DCCM) 
cij <- dccm(xyz[,ca.inds$xyz])
png('ann_full.thirdRun_dynamic_correlations.png', pointsize=10, width=3000, height=2300, res=500)
#plot(cij)
plot(cij,contour=FALSE, col.regions=bwr.colors(200), at=seq(-1,1,by=0.01))
dev.off()

# Viewing the correlations in pymol
pymol.dccm(cij, pdb, type="launch")
