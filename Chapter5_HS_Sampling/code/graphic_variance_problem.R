library(latex2exp)

## plot Qini curves
qini_table <- readRDS('results/sim_scen3/qini_table_unif.RDS')

## plot results
layout(matrix(c(1,2,3,4,4,4), ncol=3, byrow=TRUE), heights = c(0.9,0.1))

## plot ATE estimates
est_ate_ca <- readRDS('results/sim_scen3/est_ate_ca.RDS')
est_ate_het <- readRDS('results/sim_scen3/est_ate_het.RDS')
est_ate_hsca <- readRDS('results/sim_scen3/est_ate_hsca.RDS')
est_ate_unif <- readRDS('results/sim_scen3/est_ate_unif.RDS')

par(mar = c(3, 3, 3, 3), xpd = TRUE)
ate_points <- data.frame(x = rep(c(1,2,3,4), times = c(1000,1000,1000,1000)), y = c(est_ate_ca, est_ate_het, est_ate_hsca, est_ate_unif))
boxplot(y~x,ate_points, ylab = '', main = 'ATE estimators', xlab = '', xaxt = "n", cex.main=1.5, cex.axis = 1.5)
axis(1, at=1, labels=TeX('$\\hat{ATE}_{CA}$'),cex.axis=1.5, mgp=c(2,2,0))
axis(1, at=2, labels=TeX('$\\hat{ATE}_{HS}$'),cex.axis=1.5, mgp=c(2,2,0))
axis(1, at=3, labels=TeX('$\\hat{ATE}_{HSCA}$'),cex.axis=1.5, mgp=c(2,2,0))
axis(1, at=4, labels=TeX('$\\hat{ATE}$'),cex.axis=1.5, mgp=c(2,2,0))

plot((0:10)/10, c(0,apply(qini_table,2,mean)),type = 'l', ylim = c(-0.05,0.12), main = 'Qini curve', ylab = '', xlab = '', cex.main=1.5, cex.axis = 1.5)
for(i in 1:100){
  lines((0:10)/10, c(0,qini_table[i,]),col = 'red')
}
lines((0:10)/10, c(0,apply(qini_table,2,mean)), lwd = 3)

# legend("topright", inset = c(0.3, 1.2),                    
#        legend = c(TeX('$\\hat{ATE}_t\\cdot t$'), TeX('$\\ATE_t\\cdot t$')), 
#        lty = c(1,1), lwd = c(1,3), col = c('red', 'black'), ncol=2, bty = "n")

## plot AUQ values
auqcs_het <- readRDS('results/sim_scen3/auqcs_het_rf.RDS')
auqcs_unif <- readRDS('results/sim_scen3/auqcs_unif_rf.RDS')
auqct_het <- readRDS('results/sim_scen3/auqct_het_rf.RDS')
auqct_unif <- readRDS('results/sim_scen3/auqct_unif_rf.RDS')
auqcx_het <- readRDS('results/sim_scen3/auqcx_het_rf.RDS')
auqcx_unif <- readRDS('results/sim_scen3/auqcx_unif_rf.RDS')

auq_points <- data.frame(x = rep(c(1,2,3,4,5,6), times = c(1000,1000,1000,1000,1000,1000)), 
                         y = c(auqcs_het, auqcs_unif, auqct_het, auqct_unif, auqcx_het, auqcx_unif))
boxplot(y~x,auq_points, col = rep(c('orange', 'red'), 2), main = 'AUQ of uplift models', ylab = 'AUQ', xaxt = "n",
        xlab = '',mgp=c(2,2,0), cex.main=1.5, cex.axis = 1.5)
axis(1, at=1, labels=TeX('S_HS'),cex.axis=1.5, mgp=c(2,2,0))
axis(1, at=2, labels=TeX('S_Rand'),cex.axis=1.5, mgp=c(2,2,0))
axis(1, at=3, labels=TeX('T_HS'),cex.axis=1.5, mgp=c(2,2,0))
axis(1, at=4, labels=TeX('T_Rand'),cex.axis=1.5, mgp=c(2,2,0))
axis(1, at=5, labels=TeX('X_HS'),cex.axis=1.5, mgp=c(2,2,0))
axis(1, at=6, labels=TeX('X_Rand'),cex.axis=1.5, mgp=c(2,2,0))

##legend
par(mar=c(1, 1, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('red', 'black')
legend(x = "bottom",inset = 0,
       legend = c(TeX('$\\hat{ATE}_t\\cdot t$'),TeX('$\\ATE_t\\cdot t$')), 
       col=plot_colors, lty = c(1,1), lwd=c(1,3), cex=2, horiz = TRUE, bty = "n")


