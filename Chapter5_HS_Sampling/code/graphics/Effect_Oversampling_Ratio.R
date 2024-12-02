## plot results
layout(matrix(c(1,2), ncol=1, byrow=TRUE), heights = c(0.9,0.1))

par(mar = c(3, 3, 3, 3), xpd = TRUE)
## heteroskedasticity
pH <- 0.1
VH <- 6
VL <- 1
p <- (100:500)/1000
up_bound <- 1/(pH+(1-pH)*VL/VH)
in_green <- 1+(up_bound-1)/100*(0:100)
in_red <- up_bound + (max(p/pH)-up_bound)/100*(0:100)
default_variance <- pH*VH+(1-pH)*VL
opt_osratio <- sqrt(VH)/(pH*sqrt(VH)+(1-pH)*sqrt(VL))
est_osratio <- sqrt(VH*1.9)/(pH*sqrt(VH*1.9)+(1-pH)*sqrt(VL*0.6))

plot(p/pH,(pH^2)*VH/p+((1-pH)^2)*VL/(1-p), xaxt = "n", xlab = TeX(''), type = 'l', lwd = 3, ylab = TeX('$Var[\\hat{ATE}_{HS}]$'), 
     main = TeX('Variance of ATE estimate based on HS-sampling'), mgp=c(1,1,0), ylim = c(1,1.8), yaxt = 'n')
polygon(x = c(in_green, rev(in_green)), y = c((pH^2)*VH/in_green/pH+((1-pH)^2)*VL/(1-in_green*pH), rep(default_variance, length(in_green))), col = "#65BFFF")
polygon(x = c(in_red, rev(in_red)), y = c((pH^2)*VH/in_red/pH+((1-pH)^2)*VL/(1-in_red*pH), rep(default_variance, length(in_red))), col = "red")
lines(p/pH,(pH^2)*VH/p+((1-pH)^2)*VL/(1-p), lwd = 3)
lines(p/pH,(pH^2)*(VH*1.9)/p+((1-pH)^2)*(VL*0.6)/(1-p), lty = 1, lwd = 3, col = 'orange')
lines(c(1,1),c(0.95,default_variance), lty = 3, col = "#65BFFF")
axis(1, at=1, labels=TeX('1'), mgp=c(3,1.5,0))
lines(c(up_bound, up_bound),c(0.95,default_variance), lty = 3, col = "#65BFFF")
axis(1, at=up_bound, labels=TeX('$\\frac{1}{p_H+(1-p_H)\\cdot Q_{V}^{-1}}$'), mgp=c(3,2.5,0))
lines(c(opt_osratio,opt_osratio), c(0.95, (pH*sqrt(VH)+(1-pH)*sqrt(VL))^2), lty = 3)
axis(1, at=opt_osratio, labels=TeX('$R_H(Q_V)$'), mgp=c(3,1.5,0))
lines(c(est_osratio,est_osratio), c(0.95, (pH*sqrt(VH*1.9)+(1-pH)*sqrt(VL*0.6))^2), lty = 3)
axis(1, at=est_osratio, labels=TeX('$R_H(\\hat{Q}_V)$'), mgp=c(3,1.5,0))
arrows(opt_osratio, 1, est_osratio, 1)

##legend
par(mar=c(0, 0, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('black', 'orange')
legend(x = "bottom",inset = 0,
       legend = c(TeX('actual variance'),TeX('estimated variance')), 
       col=plot_colors, lty = c(1,1), lwd=c(3,3), cex=1, horiz = TRUE, bty = "n")