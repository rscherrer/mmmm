
#CARESSIM ######################################################

CAressim=function (size, tmax, bsus, dsus, bres, dres, mutrate, plotter, nrtomix){
  
  ######################################################################################## prework
  grid=data.frame()
  saver=data.frame()
  
  counter=1
  while(counter <=size){
    grid=rbind(grid, (rep(1,size)))
    counter=counter+1
  }
  
  xcoords=vector()
  counter=1
  while(counter<=size){
    line=rep (counter, size)
    xcoords=c(xcoords,line)
    counter=counter+1
  }
  ycoords=rep(1:size, size )
  
  translist=rbind(xcoords,ycoords)
  
  
  t=0 ######################## START THE SIMULATIONNNSSS   #############
  
  tobesaved=c(t, sum(grid == 1),sum(grid == 2))
  #print(tobesaved)
  saver=rbind(saver, tobesaved)
  
  
  if(plotter==1){m=data.matrix(grid) #plotter
  jpeg(paste(t, ".jpg"))
   par(mar=c(3, 3, 3,3))
  image(c(0:size), c(0:size), m, col=c("white","blue","red" ), useRaster=TRUE, axes=TRUE, zlim=c(0,2))
  dev.off()
  }
  
  
  
  while (t<tmax && (sum(grid == 1)+sum(grid==2))>0){ #add condition for extinction
    
    seqvec=sample (1:(size*size), size*size, replace=F) #sample sequence of cell updates (uniform probabilty no replacements, so all cells are updated)
    
    cell=1
    
    while(cell<=length (seqvec)){
      currentcell=seqvec[cell]
      coordcurcell=translist[,currentcell]
      #########################################
      #if (grid[coordcurcell[1],coordcurcell[2]]==0){print("empty")} #do nothing it is a dead cell what were you thinking
      
      #########################################
      if (grid[coordcurcell[1],coordcurcell[2]]==1)
      {
        event=sample ( c(-1,0,1,2),1, prob=c(dsus, (1-(dsus+bsus+mutrate)), bsus, mutrate))  #sample event
        
        if (event==-1){grid[coordcurcell[1],coordcurcell[2]]=0} #A death
        
        if (event==1) #a potential birth
        { 
          neighbor=sample(1:8, 1)
          if (neighbor==1){neighcoords=c(coordcurcell[1]-1,coordcurcell[2]+1)}
          if (neighbor==2){neighcoords=c(coordcurcell[1],coordcurcell[2]+1)}
          if (neighbor==3){neighcoords=c(coordcurcell[1]+1,coordcurcell[2]+1)}
          if (neighbor==4){neighcoords=c(coordcurcell[1]-1,coordcurcell[2])}  
          if (neighbor==5){neighcoords=c(coordcurcell[1]+1,coordcurcell[2])}
          if (neighbor==6){neighcoords=c(coordcurcell[1]-1,coordcurcell[2]-1)}
          if (neighbor==7){neighcoords=c(coordcurcell[1],coordcurcell[2]-1)}
          if (neighbor==8){neighcoords=c(coordcurcell[1]+1,coordcurcell[2]-1)}
          if(neighcoords[1]>=1 && neighcoords[1]<=size && neighcoords[2]>=1 && neighcoords[2]<=size){
            if (grid[neighcoords[1], neighcoords[2]]== 0){grid[neighcoords[1], neighcoords[2]]=1} # rejoice a suceptible is born}
          }
        }
        
        if (event==2) # a mutation
        {
          neighbor=sample(1:8, 1)
          if (neighbor==1){neighcoords=c(coordcurcell[1]-1,coordcurcell[2]+1)}
          if (neighbor==2){neighcoords=c(coordcurcell[1],coordcurcell[2]+1)}
          if (neighbor==3){neighcoords=c(coordcurcell[1]+1,coordcurcell[2]+1)}
          if (neighbor==4){neighcoords=c(coordcurcell[1]-1,coordcurcell[2])}  
          if (neighbor==5){neighcoords=c(coordcurcell[1]+1,coordcurcell[2])}
          if (neighbor==6){neighcoords=c(coordcurcell[1]-1,coordcurcell[2]-1)}
          if (neighbor==7){neighcoords=c(coordcurcell[1],coordcurcell[2]-1)}
          if (neighbor==8){neighcoords=c(coordcurcell[1]+1,coordcurcell[2]-1)}
          if(neighcoords[1]>=1 && neighcoords[1]<=size && neighcoords[2]>=1 && neighcoords[2]<=size){
            if (grid[neighcoords[1], neighcoords[2]]== 0){grid[neighcoords[1], neighcoords[2]]=2} # rejoice a resistant mutant is born}
          }
        }
      }
      
      #########################################
      
      if (grid[coordcurcell[1],coordcurcell[2]]==2)
      {
        event=sample ( c(-1,0,1),1, prob=c(dres, (1-(dres+bres)), bres))
        #########################################
        #if (grid[coordcurcell[1],coordcurcell[2]]==0){print("empty")} #do nothing it is a dead cell what were you thinking
        if (event==-1){grid[coordcurcell[1],coordcurcell[2]]=0} #A death
        
        if (event==1) #a potential birth
        { 
          neighbor=sample(1:8, 1)
          if (neighbor==1){neighcoords=c(coordcurcell[1]-1,coordcurcell[2]+1)}
          if (neighbor==2){neighcoords=c(coordcurcell[1],coordcurcell[2]+1)}
          if (neighbor==3){neighcoords=c(coordcurcell[1]+1,coordcurcell[2]+1)}
          if (neighbor==4){neighcoords=c(coordcurcell[1]-1,coordcurcell[2])}  
          if (neighbor==5){neighcoords=c(coordcurcell[1]+1,coordcurcell[2])}
          if (neighbor==6){neighcoords=c(coordcurcell[1]-1,coordcurcell[2]-1)}
          if (neighbor==7){neighcoords=c(coordcurcell[1],coordcurcell[2]-1)}
          if (neighbor==8){neighcoords=c(coordcurcell[1]+1,coordcurcell[2]-1)}
          if(neighcoords[1]>=1 && neighcoords[1]<=size && neighcoords[2]>=1 && neighcoords[2]<=size){
            if (grid[neighcoords[1], neighcoords[2]]== 0){grid[neighcoords[1], neighcoords[2]]=2} # rejoice a resistant mutant is born}
          }
        }
      }
      #########################################
      #print(cell)
      cell=cell+1
      
      
      
    }
    ######mixing######
    mixvec1=sample (1:(size*size), nrtomix, replace=F)
    mixvec2=sample (1:(size*size), nrtomix, replace=F)
    curmix=1
    
    while (curmix<=nrtomix){
      coordmix1=translist[,mixvec1[curmix]]
      coordmix2=translist[,mixvec2[curmix]]
      
      value1=grid[coordmix1[1],coordmix1[2]]
      value2=grid[coordmix2[1],coordmix2[2]]
      
      grid[coordmix1[1],coordmix1[2]]=value2
      grid[coordmix2[1],coordmix2[2]]=value1
      
      
      curmix=curmix+1
      
    }
    
    ###############
    
    t=t+1
    tobesaved=c(t, sum(grid == 1),sum(grid == 2))
    #print(tobesaved)
    saver=rbind(saver, tobesaved)
    
    
    if(plotter==1){m=data.matrix(grid) # plotter
    par(mar=c(3, 3, 3,3))
    jpeg(paste(t, ".jpg"))
    image(c(0:size), c(0:size), m, col=c("white","blue","red" ), useRaster=TRUE, axes=TRUE, zlim=c(0,2))
    dev.off()
    image(c(0:size), c(0:size), m, col=c("white","blue","red" ), useRaster=TRUE, axes=TRUE, zlim=c(0,2))
    }
    
  }
  print(tobesaved)
  return (saver)
} 



########################## END of FUNCTIONS ####################################################################



size= 50
tmax=500
bsus=0.08
dsus=0.1
bres=0.25
dres=0.1
plotter=1
maxtrails=500
mutrate=0.00007
saver=data.frame()

nrtomix=0 #no mixing

setwd("C:/Users/Timo v. Eldijk/Desktop/catest") #set directory appropriately, warning THIS CODE GENERATES A METRIC SHITTON OF JPEGS!
CAressim(size, tmax, bsus, dsus, bres, dres, mutrate, plotter, nrtomix)

nrtomix=2500 #fully mixed
setwd("C:/Users/Timo v. Eldijk/Desktop/catest") #set directory appropriately, warning THIS CODE GENERATES A METRIC SHITTON OF JPEGS!
CAressim(size, tmax, bsus, dsus, bres, dres, mutrate, plotter, nrtomix)



