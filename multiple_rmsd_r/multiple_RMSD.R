library(bio3d)
library(ncdf4)



mydcdfile_1 <- "/Users/jansiess/Desktop/main_repository/v/vik_lab/projects/ferredoxin_experiments/mutter_simulations/snn/single_chain/joint_analysis_in_R/snn_SC.firstRun.wa.dyn.dcd"
dcd_1 <- read.dcd(mydcdfile_1)
mypdb_1 <- "/Users/jansiess/Desktop/main_repository/v/vik_lab/projects/ferredoxin_experiments/mutter_simulations/snn/single_chain/joint_analysis_in_R/snn_SC.firstRun.wa.pdb"
pdb_1 <- read.pdb(mypdb_1)

ca_1.inds <- atom.select(pdb_1, elety="CA")
xyz_1 <- fit.xyz(fixed=pdb_1$xyz, mobile=dcd_1,
              fixed.inds=ca_1.inds$xyz,
              mobile.inds=ca_1.inds$xyz)
rm(pdb_1,dcd_1,mydcdfile_1,mypdb_1)
rd_1 <- rmsd(xyz_1[1,ca_1.inds$xyz], xyz_1[,ca_1.inds$xyz])



mydcdfile_2 <- "/Users/jansiess/Desktop/main_repository/v/vik_lab/projects/ferredoxin_experiments/mutter_simulations/snn/single_chain/joint_analysis_in_R/snn_SC.secondRun.wa.dyn.dcd"
dcd_2 <- read.dcd(mydcdfile_2)
mypdb_2 <- "/Users/jansiess/Desktop/main_repository/v/vik_lab/projects/ferredoxin_experiments/mutter_simulations/snn/single_chain/joint_analysis_in_R/snn_SC.secondRun.wa.pdb"
pdb_2 <- read.pdb(mypdb_2)

ca_2.inds <- atom.select(pdb_2, elety="CA")
xyz_2 <- fit.xyz(fixed=pdb_2$xyz, mobile=dcd_2,
               fixed.inds=ca_2.inds$xyz,
               mobile.inds=ca_2.inds$xyz)
rm(pdb_2,dcd_2,mydcdfile_2,mypdb_2)
rd_2 <- rmsd(xyz_2[1,ca_2.inds$xyz], xyz_2[,ca_2.inds$xyz])


###
mydcdfile_3 <- "/Users/jansiess/Desktop/main_repository/v/vik_lab/projects/ferredoxin_experiments/mutter_simulations/snn/single_chain/joint_analysis_in_R/snn_SC.thirdRun.wa.dyn.dcd"
dcd_3 <- read.dcd(mydcdfile_3)
mypdb_3 <- "/Users/jansiess/Desktop/main_repository/v/vik_lab/projects/ferredoxin_experiments/mutter_simulations/snn/single_chain/joint_analysis_in_R/snn_SC.thirdRun.wa.pdb"
pdb_3 <- read.pdb(mypdb_3)

ca_3.inds <- atom.select(pdb_3, elety="CA")
xyz_3 <- fit.xyz(fixed=pdb_3$xyz, mobile=dcd_3,
               fixed.inds=ca_3.inds$xyz,
               mobile.inds=ca_3.inds$xyz)
rm(pdb_3,dcd_3,mydcdfile_3,mypdb_3)
rd_3 <- rmsd(xyz_3[1,ca_3.inds$xyz], xyz_3[,ca_3.inds$xyz])




png('snn_SC.threeRuns_RMSD.png', pointsize=10, width=3000, height=2300, res=500)
plot(rd_1, ylim=range(c(rd_2,rd_3)), typ="l", ylab="RMSD", xlab="Frame No.", main="SNN Single Chain - Triplicate Runs",
     panel.first = grid(10,10))
points(rd_2, typ="l", col = "blue")
points(rd_3, typ="l", col = "red")
# points(lowess(c(rd_1,rd_2,rd_3)), typ="l", col="black", lty=2, lwd=3)
legend("bottomright", inset=.02, title="Simulations",
       legend=c("Run 1", "Run 2", "Run 3"), 
       col=c("black","blue", "red"), lty=1, cex=0.75
       )
dev.off()

rd_tot <- rd_1 + rd_2 + rd_3
rd_tot <- (rd_tot)/3

# RMSD Histogram
png('snn_SC.threeRuns_RMSD-histo.png', pointsize=10, width=3000, height=2300, res=500)
hist(rd_tot, breaks=40, freq=FALSE, main="RMSD Histogram - Average  ", xlab="RMSD")
lines(density(rd_tot), col="gray", lwd=3)
dev.off()


