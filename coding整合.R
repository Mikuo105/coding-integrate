links_to_Inventor_list_Long <- read.csv("C:/Users/kuo/Desktop/r/Inventor_list_long_June2020_v2.csv", header=T, as.is=T)
#---------------------------------------------------------------------------------------------------
Intersil<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Intersil")
SMSC<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="SMSC")
LSI<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="LSI")
Broadcom<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Broadcom")
Zoran<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Zoran")
Altera<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Altera")
Xilinx<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Xilinx")
Lattice<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Lattice")
PMC_Sierra<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="PMC_Sierra")
Qualcomm<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Qualcomm")
DSP_Group<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="DSP_Group")
Cirrus_Logic<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Cirrus_Logic")
SanDisk<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="SanDisk")
ARM<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="ARM")
ESS_Technology<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="ESS_Technology")
ATI_Technologies<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="ATI_Technologies")
Genesis_Microchip<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Genesis_Microchip")
Marvell<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Marvell")
Nvidia<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Nvidia")
Semtech<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Semtech")
Realtek<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Realtek")
Silicon_Labs<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Silicon_Labs")
Conexant<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Conexant")
MediaTek<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="MediaTek")
VIA<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="VIA")
Dialog<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Dialog")
Agere<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Agere")
Atheros<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Atheros")
CSR<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="CSR")
Ali<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Ali")
Avago<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Avago")
ICS<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="ICS")
SST<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="SST")
Qlogic<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Qlogic")
Globespan<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Globespan")
Virata<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Virata")
VIA_Cyrix<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="VIA_Cyrix")
SiS<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="SiS")
Sunplus<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Sunplus")
GlobespanVirata<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="GlobespanVirata")#為併購公司
Solomon<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Solomon")
Mstar<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Mstar")
Himax<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Himax")
Novatek<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Novatek")
Spreadtrum<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Spreadtrum")
Hisilicon<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="Hisilicon")
MegaChips<-subset(links_to_Inventor_list_Long,links_to_Inventor_list_Long[,11]=="MegaChips")
#-----------------------------------------------------------------------------------

x<-c("Intersil","SMSC","LSI","Broadcom","Zoran","Altera","Xilinx","Lattice","PMC_Sierra","Qualcomm","DSP_Group","Cirrus_Logic","SanDisk","ARM","ESS_Technology","ATI_Technologies","Genesis_Microchip","Marvell","Nvidia","Semtech","Realtek","Silicon_Labs","Conexant","MediaTek","VIA","Dialog","Agere","Atheros","CSR","Ali","Avago","ICS","SST","Qlogic","Globespan","Virata","VIA_Cyrix","SiS","Sunplus","GlobespanVirata","Solomon","Mstar","Himax","Novatek","Spreadtrum","Hisilicon","MegaChips")
#-----------------------------------------------------------------------------------
dir.create("C:\\Users\\kuo\\Desktop\\r\\data_firm_year_Oct2020")
for(i in 1:47) {
  dir.create(paste0("C:\\Users\\kuo\\Desktop\\r\\data_firm_year_Oct2020\\",as.character(x[i])))
}


for(i in 1:47) {
  y<-subset(Qualcomm,Qualcomm[,9]>=(1970+i-1) & Qualcomm[,9]<=(1972+i-1))
  z_a<-y[c(1,4)]
  z_b<-merge(z_a,z_a,by="Document.Number")
  z_c<-z_b[c(2,3)]
  write.csv(z_c, file = paste0("C:\\Users\\kuo\\Desktop\\r\\data_firm_year_Oct2020\\Qualcomm\\temp (",i,").csv"),row.names=F)
  cat("dl n=10000,","format=edgelist1\n","labels embedded\n","data:\n",file =paste0("C:/Users/kuo/Desktop/r/data_firm_year_Oct2020/Qualcomm/Qualcomm (",i,").txt"),append=T)
  write.table(z_c,file = paste0("C:\\Users\\kuo\\Desktop\\r\\data_firm_year_Oct2020\\Qualcomm\\Qualcomm (",i,").txt"),row.names=F,col.names=F,quote=F,append=T )
}





library(igraph)
#--------------------------Plotting Networks: Weighted Edges---------------------
#import the sample_weighted_adjmatrix file from bottom of the page:
dat=read.csv(file.choose(),header=TRUE,row.names=1,check.names=FALSE) # read .csv file
m=as.matrix(dat)
net=graph.adjacency(m,mode="undirected",weighted=TRUE,diag=FALSE) 
#--------------------------weight and edge name-----------------------------
summary(net)
E(net)$weight
V(net)$name

#--------------------------import the sample_attributes file---------------------
#a=read.csv("C:/Users/P65/Desktop/R教學/output_y_int_CNST.csv")
#V(net)$CNST=as.character(a$CNST[match(V(net)$name,a$Inventor_Parse3)])


#-----------------------------------
plot(net,
     vertex.size=3.5,
     vertex.label.cex=0.8,
     vertex.label="",
     vertex.label.dist=0.8,
     vertex.color=rainbow(10),
     edge.color="black",
     edge.width=E(net)$weight)
#------------------------------------
plot(net,
     vertex.size=3.5,
     vertex.label.cex=0.8,
     vertex.label=V(net)$CNST,
     vertex.label.dist=0.8,
     edge.color="black",
     edge.width=E(net)$weight,
     layout=layout.sphere)

#-------------------------------------
ceb <- cluster_edge_betweenness(net) 
dendPlot(ceb, mode="hclust")
plot(ceb, net,vertex.size=3.5,
     vertex.label="",
     vertex.label.cex=0.8,
     vertex.label=V(net)$CNST,
     vertex.label.dist=0.8,
     edge.color="black",
     edge.width=E(net)$weight)



#變數輸出
library(igraph)
path <- "C:/Users/P65/Desktop/R語言練習/degree_cluster_constraint/Silicon_Labs/"
files <- list.files(path = path, pattern = "Silicon_Labs*")
df1 <- data.frame()
for(file in files) {
  df1 <- read.csv(paste(path, file, sep=""),header=TRUE,row.names=1)
  m=as.matrix(df1)
  net=graph.adjacency(m,mode="undirect",weighted=TRUE,diag=FALSE) 
  summary(net)
  E(net)$weight
  V(net)$name
  #-------------------------------------------------------------------------------
  #import the sample_attributes file:
  a=read.csv("C:/Users/P65/Desktop/R語言練習/output_y_int_CNST.csv")
  V(net)$CNST=as.character(a$CNST[match(V(net)$name,a$Inventor_Parse3)])
  #-------------------------------------------------------------------------------
  y_degree<-degree(net,normalized = FALSE)
  y_degree_nor<-degree(net,normalized = TRUE)
  y_contraint<-constraint(net)
  y_transitivity_global<-transitivity(net, type="global",isolates = "nan")  # net is treated as an undirected network
  y_transitivity_local<-transitivity(net, type="local",isolates = "nan")
  CompareDegree <- cbind(y_degree,y_degree_nor,y_transitivity_global,y_transitivity_local,y_contraint,V(net)$CNST,1970+as.integer(substring(file,as.integer(nchar("Silicon_Labs"))+3,as.integer(nchar("Silicon_Labs"))+4))-1,1972+as.integer(substring(file,as.integer(nchar("Silicon_Labs"))+3,as.integer(nchar("Silicon_Labs"))+4))-1)
  head(CompareDegree)
  write.table(CompareDegree, file =("C:\\Users\\P65\\Desktop\\R語言練習\\degree_cluster_constraint\\Silicon_Labs_degree_cluster_constraint.CSV"),col.names=NA,append=TRUE,sep=",")
}




#---------------------------地圖加上網絡-------------------------------
library(network)
library(maptools)
data(wrld_simpl)
rawnodes<- read.csv("C:/Users/P65/Desktop/R教學/Country_terms_FREQ.csv", header=T, as.is=T)
names(rawnodes)
rawedges<- read.csv("C:/Users/P65/Desktop/R教學/Country_terms_COOC.csv", header=T, as.is=T)
names(rawedges)
nrow(rawedges)
Edgelist<-rawedges[rawedges$Tot_cooc>500,c('Source','Target',"Tot_cooc")]
head(Edgelist) 
nrow(Edgelist)

coocNet<-network(Edgelist,
                 matrix.type='edgelist',
                 directed=FALSE,  # this will be an undirected network
                 ignore.eval=FALSE,  # confusingly, this tells it to include edge weights
                 names.eval='Tot_cooc')  # names for the edge weight)


# attach the appropriate lat and long coordinates
# need to subset to the vertices actually in the network
coocNet%v%'lon'<-sapply(network.vertex.names(coocNet),function(name){
  rawnodes[rawnodes$ID==name,]$lon
})

coocNet%v%'lat'<-sapply(network.vertex.names(coocNet),function(name){
  rawnodes[rawnodes$ID==name,]$lat
})


# plot the map for the background
plot(wrld_simpl)

# plot the network using the geo coordinates
plot.network(coocNet,  # pass in the network
             # don't erase the map before drawing the network
             new=FALSE, 
             # get coordiantes from vertices and pass in as 2-col matrix
             coord=cbind(coocNet%v%'lon',coocNet%v%'lat'),  
             # ---- all the rest of these are optional to make it look nice ------
             # set a semi-transparent edge color
             edge.col='#AA555555',
             # specifiy an edge width scaled as fraction of total co-occurence
             edge.lwd=coocNet%e%'Tot_cooc'/500,
             # set the vertex size
             vertex.cex=0.4,
             # 節點與文字的距離
             label.pos=3,
             # 輸出節點名稱
             #label = network.vertex.names(coocNet),
             # 宣告可以彎曲節線，同時設定彎曲角度
             usecurve=TRUE,
             edge.curve=0.5,
             # set a semi transparent vertex color
             vertex.col='#AA555555',
             vertex.border='white',
             #edge attribute name
             #edge.label=TRUE
             # please don't jitter the points around
             jitter=FALSE,)

#------------------------一般網絡圖呈現
plot.network(coocNet,displaylabels=TRUE,boxed.labels=TRUE,
             vertex.cex=0,
             label.pos=5,
             label.cex=0.6,
             edge.lwd=coocNet%e%'Tot_cooc'/500,
             edge.col='#AA555555')




#=========================動態=========================
#this version of the script has been tested on igraph 1.0.1
#load libraries
library(igraph)
library(RColorBrewer)
display.brewer.all() 
#load the edges with time stamp
#there are three columns in edges: id1,id2,time
edges <- read.table("C:/Users/P65/Desktop/R教學/edges.csv",header=T)

#generate the full graph
g <- graph.data.frame(edges,directed=F)

#generate a cool palette for the graph (darker colors = older nodes)
YlOrBr.pal <- colorRampPalette(brewer.pal(8,"YlOrRd"))
#colors for the nodes are chosen from the very beginning
V(g)$color <- rev(YlOrBr.pal(vcount(g)))[as.numeric(V(g)$name)]

V(g)$color
vcount(g)
#time in the edges goes from 1 to 300. We kick off at time 3
ti <- 3
#remove edges which are not present
gt <- delete_edges(g,which(E(g)$time > ti))
#generate first layout using graphopt with normalized coordinates. This places the initially connected set of nodes in the middle. If you use fruchterman.reingold it will place that initial set in the outer ring.
layout.old <- norm_coords(layout.graphopt(gt), xmin = -1, xmax = 1, ymin = -1, ymax = 1)

V(g)$name
#total time of the dynamics
total_time <- max(E(g)$time)
#This is the time interval for the animation. In this case is taken to be 1/10
#of the time (i.e. 10 snapshots) between adding two consecutive nodes
dt <- 0.1
#Output for each frame will be a png with HD size 1600x900 :)
png(file="C:\\Users\\P65\\Desktop\\R教學\\Temporal_networks\\example%03d.png", width=1600,height=900)
#Time loop starts
for(time in seq(3,total_time,dt)){
  #remove edges which are not present
  gt <- delete_edges(g,which(E(g)$time > time))
  #with the new graph, we update the layout a little bit
  layout.new <- layout_with_fr(gt,coords=layout.old,niter=10,start.temp=0.05,grid="nogrid")
  #plot the new graph
  plot(gt,layout=layout.new,
       vertex.label=V(g)$name,vertex.size=1+2*log(degree(gt)),
       vertex.frame.color=V(g)$color,edge.width=1.5,
       asp=9/16,margin=-0.15)
  #use the new layout in the next round
  #use the new layout in the next round
  layout.old <- layout.new
}
dev.off()