# compare test result
# which model is best?
# plot success (=logloss) over size (=numer of trees)


res.test <- matrix(NA, ncol=5 , nrow=38)
# cols : 
# 1) method 1:4
# 2) dissim 1:4
# 3) parameter 0..1
# 4) size 1:500
# 5) logloss

# Chipman 1
{res.test[1,] <- c(1,1,0.2,NA,NA)
res.test[2,] <- c(1,1,0.4,NA,NA)
}
# Chipman 2

{res.test[3,] <- c(2,1,0,49.6,0.4284)
res.test[4,] <- c(2,2,0.3,40.2,0.4326)
res.test[5,] <- c(2,3,0.4,45.1,0.4446)
res.test[6,] <- c(2,4,0.4,49.3,0.4347)

res.test[7,] <- c(2,1,0.5,6.2,0.5871)
res.test[8,] <- c(2,2,0.9,5.7,0.6755)
res.test[9,] <- c(2,3,0.875,6.1,0.5841)
res.test[10,] <- c(2,4,0.875,6.9,0.7558)
}

# Meiner
{res.test[11,] <- c(3,1,0,49.9,0.4278)
res.test[12,] <- c(3,2,0.1,69.1,0.4266)
res.test[13,] <- c(3,4,0.1,69,0.4267)
# res.test[14,] <-c(3,2,0.2,NA,NA) # Meiner , d1 , 0.2
res.test[15,] <- c(3,3,0.2,44,0.4287)
res.test[16,] <- c(3,4,0.15,48.4,0.4281)
res.test[17,] <- c(3,1,0.5,5.7,0.5964)
res.test[18,] <- c(3,2,0.8,7.1,0.5450)
res.test[19,] <- c(3,3,0.7,7.1,0.5583)
res.test[20,] <- c(3,4,0.6,8.8,0.5201)
}

# high performers

{res.test[21,] <- c(4,0,NA,50,0.4289)
res.test[22,] <- c(4,0,NA,45,0.4299)
res.test[23,] <- c(4,0,NA,40,0.4315)

res.test[24,] <- c(4,0,NA,10,0.5297)
res.test[25,] <- c(4,0,NA,9,0.5466)
res.test[26,] <- c(4,0,NA,8,0.5672)
res.test[27,] <- c(4,0,NA,7,0.5979)
res.test[28,] <- c(4,0,NA,6,0.6469)
res.test[29,] <- c(4,0,NA,5,0.7135)
}
# regular small
{res.test[30,] <- c(5,0,0,5,0.7873)
res.test[31,] <- c(5,0,0,6,0.7007)
res.test[32,] <- c(5,0,0,7,0.6450)
res.test[33,] <- c(5,0,0,8,0.6041)
res.test[34,] <- c(5,0,0,9,0.5739)
res.test[35,] <- c(5,0,0,10,0.5503)

res.test[36,] <- c(5,0,0,40,0.43224932)
res.test[37,] <- c(5,0,0,45,0.43013630)
res.test[38,] <- c(5,0,0,50,0.4289)}

plot(res.test[,4:5]
     , xlab='size'
     , ylab='success'
     , col=res.test[,1]
     , pch=res.test[,2])
legend('topright', legend=c('d0','d1','d2','Shannon-Banks','Chipman 1','Chipman 2','Meiner','high performers','regular small'), col=c(1,1,1,1,1:5), pch=c(1:4,rep(15,5)))

plot(res.test[res.test[,4]<11,4:5]
     , xlab='size'
     , ylab='success'
     , col=res.test[res.test[,4]<11,1]
     , pch=res.test[res.test[,4]<11,2])

plot(res.test[res.test[,4]>10 & res.test[,4]< 60,4:5]
     , xlab='size'
     , ylab='success'
     , col=res.test[res.test[,4]>10 & res.test[,4]< 60,1]
     , pch=res.test[res.test[,4]>10 & res.test[,4]< 60,2])
