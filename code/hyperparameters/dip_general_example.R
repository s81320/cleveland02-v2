#### handmade example : how the dip statistic works ####
########################################################

# creating coordinates 
# option 1) uniform distribution , bounded
#x1 <- runif(20,1,3)
#x2 <- runif(20,5,7)
#y1 <- runif(20,0,1)
#y2 <- runif(20,-1,0)

# option 2) normal distribution , overlapping
x1 <- rnorm(200,0,1)
x2 <- rnorm(200,5,1)
y1 <- rnorm(200,0,1)
y2 <- rnorm(200,3,1)

o1 <- order(sqrt((x1-mean(x1))^2+(y1-mean(y1))^2))
x1 <- x1[o1]
y1 <- y1[o1]
o2 <- order(sqrt((x2-mean(x2))^2+(y2-mean(y2))^2))
x2 <- x2[o2]
y2 <- y2[o2]

plot(x1,y1)
plot(x1,y1,type='n')
text(x1,y1,1:200)

# continue for both options
all.points <- rbind(cbind(x1,y1), cbind(x2,y2))
dm <- as.matrix(dist(rbind(cbind(x1,y1), cbind(x2,y2))))

#### example cluster , bimodal
plot(x1,y1, xlim=c(min(x1,x2)-0.2,max(x1,x2)+0.2)
     , ylim=c(min(y1,y2)-0.2, max(y1,y2)+0.2))
points(x2,y2)


col=rep('grey',nrow(all.points))
i <- 1
hist(dm[i,], breaks=20, main=i)
dip.test(dm[i,])

col=rep('grey',nrow(all.points))
col[dm[i,]<2] <- 'black'
col[i] <- 'red'
#col[c(27,29,30)] <- 'blue'
plot(all.points, col=col)

doc <- data.frame(matrix(0,nrow(dm),2))
for(i in 1:nrow(dm)){
  dt <-   dip.test(dm[i,])
  doc[i,] <-  c(dt$statistic, dt$p.value)
}
names(doc) <- c('dip','p')

dqu <- quantile(doc$dip,0.9)
dqu <- 0.1
dqu
plot(all.points 
     , col=ifelse(doc$dip>dqu,'blue','black')
     , main='large dip in blue, small p in green')
points(all.points+c(0.05,0.05) 
       , col=ifelse(doc$p<0.00001,'green','black') 
       , pch=ifelse(doc$p<0.00001,'x',''))
#summary(doc$p)

doc$p %>% hist(breaks=20)

which(doc$dip > dqu) %>% c(.,doc$p[.])
which(doc$dip > dqu)

j <- 1
dm[j,] %>% hist(breaks=20)
dm[j,] %>% ecdf %>% plot()
dip.test(dm[j,])
