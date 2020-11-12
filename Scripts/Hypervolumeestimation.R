install.packages("dynRB")
install.packages("sfsmisc")
library("dynRB")
library(ggplot2)
library(ggpubr)


data(finch)
head(finch)[,1:5]
finch2<-finch[1:32,]
table(finch2$Species)
View(finch)
# scatter plot
ggplot(finch2, aes(x=TarsusL, y=UBeakL, shape=Species, color=Species)) +
  geom_point()

# denisty plot
finch3<-finch2[1:19,]
finch4<-finch2[20:32,]



#subset by trait 
a<-finch3$BodyL
b<-finch4$BodyL
# simulate samples form 4 species which would be compared 



# define limits of a common grid, adding a buffer so that tails aren't cut off
lower <- min(c(a, b,d,e)) 
upper <- max(c(a, b,d,e))
plot(density(a,from = lower,to = upper,n=1000))
lines(density(a,from = lower,to = upper,n=1000))



lines(density(b,from = lower,to = upper,n=1000))
lines(density(d,from = lower,to = upper,n=1000))
lines(density(e,from = lower,to = upper,n=1000))

c<-runif(1000,lower,upper)
lines(density(c,from = lower,to = upper,n=1000))
# generate kernel densities from lowest observed to highest pbserved value
da <- density(a,from = lower,to = upper,n=1000)
db <- density(b,from = lower,to = upper,n=1000)
dc<- density(c,from = lower,to = upper,n=1000)
dd<- density(d,from = lower,to = upper,n=1000)
de<-density(e,from = lower,to = upper,n=1000)
d <- data.frame(x=da$x, a=da$y, b=db$y,c=dc$y,d=dd$y,e=de$y)



#calculate Volume integrating area under uniform curve for each species 
#D & E should be equal? potential for error propagation?
Volume_A<- sum(pmin(d$c, d$a) * diff(d$x[1:2]))
Volume_B<- sum(pmin(d$c, d$b) * diff(d$x[1:2]))
Volume_D<- sum(pmin(d$c, d$d) * diff(d$x[1:2]))
Volume_E<- sum(pmin(d$c, d$e) * diff(d$x[1:2]))

#calculate overlaps
#this is not port and not standardized to the volume metric 
overlap_AB<-sum(pmin(d$a, d$b) * diff(d$x[1:2]))
overlap_AD<-sum(pmin(d$a, d$d) * diff(d$x[1:2]))
overlap_AE<-sum(pmin(d$a, d$e) * diff(d$x[1:2]))
overlap_BD<-sum(pmin(d$b, d$d) * diff(d$x[1:2]))
overlap_BE<-sum(pmin(d$b, d$e) * diff(d$x[1:2]))
overlab_DE<-sum(pmin(d$d, d$e) * diff(d$x[1:2]))

Overlap_list<-c(overlap_AB,overlap_AD,overlap_AE,overlap_BD,overlap_BE,overlab_DE)
View(Overlap_list)

#redeveloping the portmetric

#integrate area under both curves and the uniform distribution
overlap_AB_C<-sum(pmin(d$a, d$b,d$c) * diff(d$x[1:2]))
overlap_AD_C<-sum(pmin(d$a, d$d,d$c) * diff(d$x[1:2]))
overlap_AE_C<-sum(pmin(d$a, d$e,d$c) * diff(d$x[1:2]))
overlap_BD_C<-sum(pmin(d$b, d$d,d$c) * diff(d$x[1:2]))
overlap_BE_C<-sum(pmin(d$b, d$e,d$c) * diff(d$x[1:2]))
overlab_DE_C<-sum(pmin(d$d, d$e,d$c) * diff(d$x[1:2]))

Overlap_list_C<-c(overlap_AB_C,overlap_AD_C,overlap_AE_C,overlap_BD_C,overlap_BE_C,overlab_DE_C)
View(cbind(Overlap_list,Overlap_list_C))
#port is overlap controlled by uniform prob, divided by volume 
overlap_AD_C/Volume_A
overlap_AD_C/Volume_D




















traildensity<-ggdensity(finch2, x = "BodyL",
          add = "mean", rug = TRUE,
          color = "Species", fill = "Species")

traildensity$coordinates$setup_data



ggdensity(finch, x = "WingL",
          add = "mean", rug = TRUE,
          color = "Species", fill = "Species")

#density function 
dens<-density(finch$WingL[136:146],n=1000)

dens$bw

dens<-density(finch$WingL[136:146],n=1000, bw="SJ-ste")


# dy/dx first derivative
first<-diff(dens$y)/diff(dens$x)
# Second derivative
second<-diff(first)/diff(dens$x[1:999])
# Condition for inflection point
flections<-c()
for(i in 2:length(second)){
  if(sign(second[i])!=sign(second[i-1])){
    flections<-c(flections,i)
  }
}
#condition for local maxima
maxima<-c()
for(i in 2:length(first)){
  if(sign(first[i])<sign(first[i-1])){
    maxima<-c(maxima,i)
  }
}



plot(density(finch$WingL[136:146],n=1000))

abline(v=dens$x[flections],lty=2)
abline(v=dens$x[maxima],lty=3,col="red")


flections
Nmaxima<- length(flections)/2


diff(flections)


# Show the contour only
ggplot(finch, aes(x=BodyL, y=WingL,group=Species,colour=Species) ) +
  geom_density_2d(aes(alpha=..level..),contour_var = "ndensity")+
  theme_bw()+
  geom_point()+
  labs(title="Plot 2")



