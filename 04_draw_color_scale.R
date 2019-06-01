par(mfrow=c(2,1))
outerLim <- 15
bkPrezero <- seq((-1 * outerLim)-0.01,to=-0.01,0.01)
bkPostzero <- seq(0.01,(outerLim+0.01),0.01)

bks <- c(bkPrezero,bkPostzero)
colPrezero <- colorRampPalette(c("darkblue", "aliceblue"))(length(bkPrezero)-1)
colPostzero <- colorRampPalette(c("salmon", "darkred"))(length(bkPostzero)-1)
colSeq <- c(colPrezero, "white",colPostzero)

plot(x = abs(bks)[-1], y = ifelse(bks[-1] > 0, 0.3, -0.3) , col = colSeq, ylim=c(-1,1), pch = 15, ylab = NA, 
     xlab = "A. Dudley et al. asymmetric color scale", 
     yaxt = 'n', xaxt = 'n'
     )

y <- c(0,5,10,15)
axis(3, at=y,labels=y)
axis(1, at=y,labels=-y)


colPrezero <- colorRampPalette(c("darkblue", "aliceblue"))(length(bkPrezero)-1)
colPostzero <- colorRampPalette(c("#fff8f0", "darkred"))(length(bkPostzero)-1)
colSeq <- c(colPrezero, "white",colPostzero)

plot(x = abs(bks)[-1], y = ifelse(bks[-1] > 0, 0.3, -0.3) , col = colSeq, ylim=c(-1,1), pch = 15, ylab = NA, yaxt = 'n',
     xlab = "B. Symmetric color scale", xaxt = 'n'
     )
y <- c(0,5,10,15)
axis(3, at=y,labels=y)
axis(1, at=y,labels=-y)

