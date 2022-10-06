## Plots of summarised output dataframes produced by ecospat_output_summation.R File pathways made relative to ".Data/" need to be checked##

## Read in dataframes from ecospat_output_summation.R or have in environment to begin#

par(mar=c(2,2,2,2))

# line plots of buffer extent (x) against key metrics (y)- Sensitivity analyses for buffer size
tiff(filename = "/.Data//Europe_ebba_species/buffer_SA_Metrics_test_n=35species.tiff",width=1400,height=2000,res=300)

par(mfrow=c(3,2))

plot(bufferdf$mean.D ~ bufferdf$Buffer,type='l',xlab="Buffer/km",ylab="mean Schoener's D",ylim=c((min(bufferdf$mean.D)-2*max(bufferdf$D.sd)),(max(bufferdf$mean.D)+2*max(bufferdf$D.sd))))
lines(bufferdf$mean.D-bufferdf$D.sd ~ bufferdf$Buffer,lty=3)
lines(bufferdf$mean.D+bufferdf$D.sd ~ bufferdf$Buffer,lty=3)

plot(bufferdf$mean.I ~ bufferdf$Buffer,type='l',xlab="Buffer/km",ylab="mean I",ylim=c((min(bufferdf$mean.I)-2*max(bufferdf$I.sd)),(max(bufferdf$mean.I)+2*max(bufferdf$I.sd))))
lines(bufferdf$mean.I-bufferdf$I.sd ~ bufferdf$Buffer,lty=3)
lines(bufferdf$mean.I+bufferdf$I.sd ~ bufferdf$Buffer,lty=3)

plot(bufferdf$mean.stability ~ bufferdf$Buffer,type='l',xlab="Buffer/km",ylab="mean niche stability",ylim=c((min(bufferdf$mean.stability)-2*max(bufferdf$stability.sd)),(max(bufferdf$mean.stability)+2*max(bufferdf$stability.sd))))
lines(bufferdf$mean.stability-bufferdf$stability.sd ~ bufferdf$Buffer,lty=3)
lines(bufferdf$mean.stability+bufferdf$stability.sd ~ bufferdf$Buffer,lty=3)

plot(bufferdf$mean.unfilling ~ bufferdf$Buffer,type='l',xlab="Buffer/km",ylab="mean niche unfilling",ylim=c((min(bufferdf$mean.unfilling)-2*max(bufferdf$unfilling.sd)),(max(bufferdf$mean.unfilling)+2*max(bufferdf$unfilling.sd))))
lines(bufferdf$mean.unfilling-bufferdf$unfilling.sd ~ bufferdf$Buffer,lty=3)
lines(bufferdf$mean.unfilling+bufferdf$unfilling.sd ~ bufferdf$Buffer,lty=3)

plot(bufferdf$mean.expansion ~ bufferdf$Buffer,type='l',xlab="Buffer/km",ylab="mean niche expansion",ylim=c((min(bufferdf$mean.expansion)-2*max(bufferdf$expansion.sd)),(max(bufferdf$mean.expansion)+2*max(bufferdf$expansion.sd))))
lines(bufferdf$mean.expansion-bufferdf$expansion.sd ~ bufferdf$Buffer,lty=3)
lines(bufferdf$mean.expansion+bufferdf$expansion.sd ~ bufferdf$Buffer,lty=3)

plot(bufferdf$mean.niche.centroid.shift ~ bufferdf$Buffer,type='l',xlab="Buffer/km",ylab="mean niche centroid shift",ylim=c((min(bufferdf$mean.niche.centroid.shift)-2*max(bufferdf$niche.centroid.shift.sd)),(max(bufferdf$mean.niche.centroid.shift)+2*max(bufferdf$niche.centroid.shift.sd))))
lines(bufferdf$mean.niche.centroid.shift-bufferdf$niche.centroid.shift.sd ~ bufferdf$Buffer,lty=3)
lines(bufferdf$mean.niche.centroid.shift+bufferdf$niche.centroid.shift.sd ~ bufferdf$Buffer,lty=3)

dev.off()

tiff(filename = "/.Data//Europe_ebba_species/buffer_SA_p-value_test.tiff",width=1400,height=2000,res=300)

par(mfrow=c(3,2))

plot(bufferdf$mean.p.D ~ bufferdf$Buffer,type='l',xlab="Buffer/km",ylab="meanSchoener's D p-value",ylim=c((min(bufferdf$mean.p.D)-2*max(bufferdf$p.D.sd)),(max(bufferdf$mean.p.D)+2*max(bufferdf$p.D.sd))))
lines(bufferdf$mean.p.D-bufferdf$p.D.sd ~ bufferdf$Buffer,lty=3)
lines(bufferdf$mean.p.D+bufferdf$p.D.sd ~ bufferdf$Buffer,lty=3)

plot(bufferdf$mean.p.I ~ bufferdf$Buffer,type='l',xlab="Buffer/km",ylab="mean I p-value",ylim=c((min(bufferdf$mean.p.I)-2*max(bufferdf$p.I.sd)),(max(bufferdf$mean.p.I)+2*max(bufferdf$p.I.sd))))
lines(bufferdf$mean.p.I-bufferdf$p.I.sd ~ bufferdf$Buffer,lty=3)
lines(bufferdf$mean.p.I+bufferdf$p.I.sd ~ bufferdf$Buffer,lty=3)

plot(bufferdf$mean.p.stability ~ bufferdf$Buffer,type='l',xlab="Buffer/km",ylab="mean niche stability p-value",ylim=c((min(bufferdf$mean.p.stability)-2*max(bufferdf$p.stability.sd)),(max(bufferdf$mean.p.stability)+2*max(bufferdf$p.stability.sd))))
lines(bufferdf$mean.p.stability-bufferdf$p.stability.sd ~ bufferdf$Buffer,lty=3)
lines(bufferdf$mean.p.stability+bufferdf$p.stability.sd ~ bufferdf$Buffer,lty=3)

plot(bufferdf$mean.p.unfilling ~ bufferdf$Buffer,type='l',xlab="Buffer/km",ylab="mean niche unfilling p-value",ylim=c((min(bufferdf$mean.p.unfilling)-2*max(bufferdf$p.unfilling.sd)),(max(bufferdf$mean.p.unfilling)+2*max(bufferdf$p.unfilling.sd))))
lines(bufferdf$mean.p.unfilling-bufferdf$p.unfilling.sd ~ bufferdf$Buffer,lty=3)
lines(bufferdf$mean.p.unfilling+bufferdf$p.unfilling.sd ~ bufferdf$Buffer,lty=3)

plot(bufferdf$mean.p.expansion ~ bufferdf$Buffer,type='l',xlab="Buffer/km",ylab="mean niche expansion p-value",ylim=c((min(bufferdf$mean.p.expansion)-2*max(bufferdf$p.expansion.sd)),(max(bufferdf$mean.p.expansion)+2*max(bufferdf$p.expansion.sd))))
lines(bufferdf$mean.p.expansion-bufferdf$p.expansion.sd ~ bufferdf$Buffer,lty=3)
lines(bufferdf$mean.p.expansion+bufferdf$p.expansion.sd ~ bufferdf$Buffer,lty=3)



dev.off()



par(mfrow=c(1,1))
#plot niche vs range shifts#
tiff(filename = "/.Data//Europe_ebba_species/preliminary niche vs. range shift plot.tiff",width=1000,height=1000,res=150)

plot(df$RangeShiftMagnitude ~ df$niche.centroid.shift, main="All Species niche vs. range centroid shift",xlab="Niche centroid shift (PCA axis scale)",ylab="Range centroid shift (km)",pch=20)

dev.off()

#plot grouped niche vs range shifts#
tiff(filename = "/.Data//Europe_ebba_species/preliminary grouped niche vs. range shift plot (WRONG).tiff",width=1000,height=1000,res=150)

plot(GrpSmry$mean.niche.centroid.shift ~ GrpSmry$range.shift.mag,pch=20,xlim=c(0,250),ylim=c(0,2),xlab="Range centroid shift (km)",ylab="niche centroid shift (PCA units)",main="Niche vs. range shift of groups")

arrows(GrpSmry$range.shift.mag-GrpSmry$range.shift.mag.sd, GrpSmry$mean.niche.centroid.shift, GrpSmry$range.shift.mag+GrpSmry$range.shift.mag.sd, GrpSmry$mean.niche.centroid.shift, length=0.05, angle=90, code=3)#Add horizontal error bars

arrows(GrpSmry$range.shift.mag, GrpSmry$mean.niche.centroid.shift-GrpSmry$niche.centroid.shift.sd, GrpSmry$range.shift.mag, GrpSmry$mean.niche.centroid.shift+GrpSmry$niche.centroid.shift.sd, length=0.05, angle=90, code=3)#Add vertical error bars

dev.off()

#plot changes in expansion, contraction and stability in niche vs real space#
tiff(filename = "/.Data//Europe_ebba_species/Range vs. niche proportion changes.tiff",width=1000,height=1000,res=150)

par(mfrow=c(2,2))

plot(df$expansion ~ df$Range_Expansion,main="Expansion",ylab="Niche expansion proportion",xlab="Range expansion proportion",pch=20,col="green",xlim=c(0,1),ylim=c(0,1))

plot(df$stability ~ df$Range_Stability,main="Stability",ylab="Niche stability proportion",xlab="Range stability proportion",pch=20,col='blue',xlim=c(0,1),ylim=c(0,1))

plot(df$unfilling ~ df$Range_Contraction,main="Contraction/Unfilling",ylab="Niche unfilling proportion",xlab="Range contraction proportion",pch=20,col='red',xlim=c(0,1),ylim=c(0,1))

dev.off()