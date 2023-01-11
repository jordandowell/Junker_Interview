install.packages("dynRB")
library("dynRB")
data(finch)


View(finch)



#create a hypervolume object to store results 
Hypervolumeobject<-list()


#get species names 
SpeciesList<-unique(finch$Species)
#get trait names
TraitList<-names(finch[,-1])
Hypervolumeobject<-as.list(TraitList)
#for loop across traits

for (i in 1:length(TraitList)) {
  Trait<-TraitList[i]
  print(Trait)


#get upper and lower bound for the trait
lower<-min(finch[,Trait])
upper<-max(finch[,Trait])
#create uniform KDE bound 
density_data<-data.frame()
Uniform<-density(runif(n=1000,lower,upper),from = lower,to = upper,n=1000)
density_data<-data.frame(Uniform$x,Uniform$y)

#create lists to store volume
Volume<-list()
absoluteOverlap<-data.frame()
absoluteOverlap<-data.frame(matrix(0, ncol = length(SpeciesList), nrow = length(SpeciesList)),row.names = SpeciesList)
colnames(absoluteOverlap)<-SpeciesList

Uniformnormalized_Overlap<-absoluteOverlap

#create empty columns 

#for loop across species

for (j in 1:length(SpeciesList)) {
  Species_interest<-SpeciesList[j]
  print(Species_interest)
  #calculate denisty for each speciesxtrait
  density_data[,Species_interest] <- (density(finch[finch[,1]%in%Species_interest, Trait],from = lower,to = upper,n=1000))$y
  
  
  #calculate Volume integrating area under uniform curve for each species 
  #potential for error propogation
  Volume[Species_interest]<- sum(pmin(density_data[,Species_interest], density_data[,2]) * diff(density_data[,1][1:2]))
  



}
#calculate overlap and store
for (j in 1:length(SpeciesList)) {
  Species_interest<-SpeciesList[j]
for (l in 1:length(SpeciesList)) {
  
  Species_comparison<-SpeciesList[l]
  
  #calculate overlab
  absoluteOverlap[rownames(absoluteOverlap)%in%Species_interest,Species_comparison]<-sum(pmin(density_data[,Species_interest], density_data[,Species_comparison]) * diff(density_data[,1][1:2]))
  
}
}

#calculate unifor normalized overlap and store
for (j in 1:length(SpeciesList)) {
  Species_interest<-SpeciesList[j]
  for (l in 1:length(SpeciesList)) {
    
    Species_comparison<-SpeciesList[l]
    
    #calculate overlab The Diagonal is the volume of the species 
    Uniformnormalized_Overlap[rownames(Uniformnormalized_Overlap)%in%Species_interest,Species_comparison]<-sum(pmin(density_data[,2],density_data[,Species_interest], density_data[,Species_comparison]) * diff(density_data[,1][1:2]))
    
  }
}

##port is overlap controlled by uniform prob, divided by volume
Port<-Uniformnormalized_Overlap
Volumesfordivision<-diag(as.matrix(Uniformnormalized_Overlap))

for (m in 1:length(SpeciesList)) {
  Port[,m]<-Port[,m]/Volumesfordivision[m]
}

#save dataframes to hypervolume object 

Hypervolumeobject[[Trait]]$'Proportion'<-Port
Hypervolumeobject[[Trait]]$'Volume'<-Volumesfordivision
Hypervolumeobject[[Trait]]$'AbsoluteOverlap'<-absoluteOverlap
Hypervolumeobject[[Trait]]$'UniformOverlap'<-Uniformnormalized_Overlap
}

#gather Proportion to average 
Proportionforaveraging<-list()
for (m in 1:length(TraitList)) {
  Trait<-TraitList[m]
  Proportionforaveraging[[m]]<- Hypervolumeobject[[Trait]]$'Proportion'
}


#mean
View(apply(abind::abind(Proportionforaveraging, along = 3), 1:2, mean))

#geometric mean (na.rm==T)
View(apply(abind::abind(Proportionforaveraging, along = 3), 1:2, FUN = function(x){exp(mean(log(x)))}))

#product

View(apply(abind::abind(Proportionforaveraging, along = 3), 1:2, prod))

#best is sum divided by number of traits...each value is a probability so consider each independent and divide by total dimensions
View(apply(abind::abind(Proportionforaveraging, along = 3), 1:2, FUN= function(x){sum(x/length(TraitList))}  ))


#gather volumes to average 
Volumesforaveraging<-list()
for (m in 1:length(TraitList)) {
  Trait<-TraitList[m]
  Volumesforaveraging[[m]]<- Hypervolumeobject[[Trait]]$'Volume'
}


#mean best considering each value is area under the curve 
View(apply(abind::abind(Volumesforaveraging, along = 2), 1, mean))

#geometric mean (na.rm==T)
View(apply(abind::abind(Volumesforaveraging, along = 2), 1, psych::geometric.mean))

#product makes no sense considering surface has sparsity 
View(apply(abind::abind(Volumesforaveraging, along = 2), 1, prod))

#best is sum divided by number of traits...each value is a probability so consider each independent and divide by total dimensions
View(apply(abind::abind(Volumesforaveraging, along = 2), 1, FUN= function(x){sum(x/length(TraitList))}  ))


Hypervolumeobject$
