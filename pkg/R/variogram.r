svariog<-function(x, y, z, u) {
#********************************************************************
# Calculates the empirical semi-variogram based on smoothing with
# bins. An alternative function 'variogram' is provided in library
# geoR.
#********************************************************************

if (((length(x) == length(y)) & (length(x) == length(z)))== F)
     stop("Location and Z vectors must be of same length")

n<- length(x)
mid <- rep(0,length(u)-1)
dvec<-rep(0,n*(n-1)/2)
zvec<-rep(0,n*(n-1)/2)
zdist<-rep(0,n*(n-1)/2)
m <- 0
for (j in 2:n) {
for (k in 1:(j-1)) {
m<-m+1
dvec[m] <- sqrt((x[j] - x[k])**2 + (y[j] - y[k])**2)
zdist[m] = abs(z[j]-z[k])
zvec[m] <- ((z[j] - z[k])**2)/2 }}
cut.points <- cut(dvec,u)
vario <- tapply(zvec,cut.points,mean)
freq <- tapply(zvec,cut.points,length)
robvec=sqrt(zdist)
gr=tapply(robvec,cut.points,mean)
gr=gr**4/(0.457+0.494/freq); gr=gr/2
gm=tapply(robvec,cut.points,median)
gm=(gm**4)/(2*0.457)
for (l in 2:length(u)) {
mid[l-1] = (u[l]+u[l-1])/2 }
list(classical=vario, median=gm, robust=gr, freq=freq, mid=mid, zcloud=zvec, dcloud=dvec)
}

