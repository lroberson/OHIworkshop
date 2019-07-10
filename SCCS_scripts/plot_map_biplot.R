######################################
#### Code written by Felipe Suarez ####
#### Version 1:  01-04-2019       ####
######################################

####Clear workspace

rm(list = ls())

#install and call packages

library(SpaDES)
library(raster)
library(spatial.tools)
library(rgdal)
library(ggplot2)
library(rosm)
library(ggspatial)#requires ggplot2 version >= 2.2.1.9000
library("ggpubr")
library("gridExtra")


# Set working directory
setwd("//uq.edu.au/UQ-Research/LOGANROAD-A1321") #set your working directory


#Let's import some data
myExpl<-stack("./MyExpl_LCC2.grd") # you use "stack" to import a RasterStack, which is a collection of RasterLayer objects with the same spatial extent and resolution
#if you are working with only one raster and no a collection of rasters, you can use the function
#raster (e.g.raster("myraster.tif"))
LCC<-readOGR(dsn = ".", layer = "Logan_lga_GDA") # import a shapefile

#define your study area
extent_Logan<-extent(LCC) #use "extent" to define the extent of interest
myExplLCC<-crop(myExpl,extent_Logan) #crop your raster to your extent of interest
myExplLCC<-stack(myExplLCC) #transform to raster stack again (it was transformed to a rasterbrick when cropped)

#now let's extract some values from the raster dataset (this could take a few minutes)
myExplLCC$FPC[is.na(myExplLCC$FPC[])] <- 100 #no needed but I want to work with numeric values here
FPC<- raster::extract(myExplLCC$FPC, LCC, cellnumbers=TRUE)#extract FPC values
FPC<-as.data.frame(FPC)#convert to dataframe
elevation<-raster::extract(myExplLCC$elevation,LCC,cellnumbers = TRUE)#extract elevation values
elevation<-as.data.frame(elevation)

# extract the coordinates for each cell
coordinates_FPC<-as.data.frame(xyFromCell(myExplLCC$FPC,FPC[,1]))

#create a data frame with the coordinates and the cell values from your raster
toy_df<-data.frame(x=coordinates_FPC["x"],y=coordinates_FPC["y"],
                   FPC=as.numeric(FPC$value),elevation=elevation$value)

# Now let's add a column showing different associations between variables
# Here, we will categorize areas based on the association between temperature
# and elevation (high p- high e, high p low e, low p high e, low p low e)

#subset areas with low elevation (less than 100 masl)
lowelev<-toy_df[which(toy_df$elevation<100),]
#subset areas with high elevation (more than 100 masl)
highelev<-toy_df[which(toy_df$elevation>100),]
# now classify FPC values as high or low (based on the median) in each subset of elevation
lowelev$relationship<-with(lowelev,ifelse(FPC < median(toy_df$FPC), 1,2))
highelev$relationship<-with(highelev,ifelse(FPC < median(toy_df$FPC), 3,4))
toy_df<-rbind(lowelev,highelev) #your dataframe should have now a column specifiying the categories of association
    
head(toy_df)

#How about if we want to have more control with the colours
# to create a biplot showing a continous association?

# the trick is getting the colours exactly the way you want them
# set the rgb values of three corners and then interpolate the colour grid
#I have created a function for this.
# you need three vectors specifying the colors in RGB values. Just check the color code
# and divide it by 255. For example, the RGB for yellow is 255,255,0, so the vector
#you will use in the function is 1,1,0
#you also need to specify how big do you want your matrix of colours. This is controlled
#by the parameters cells.x and cells.y. Values of 20 for these parameters will create a
#matrix of 400 (20 X 20) cells with the gradient of colours you chose.

#call the function, no need to change anything in the next few lines
biplot_values<-function(data, # data frame, it must have at least two columns with the values of the variables you want to plot. Each row is a cell of your raster
                        minvaluex, #the minimum value of the variable you want to plot on the x axis
                        minvaluey, #the minimum value of the variable you want to plot on the y axis
                        maxvaluex, #the maximum value of the variable you want to plot on the x axis
                        maxvaluey, #the maximum value of the variable you want to plot on the y axis
                        col.ll, #vector of color that will be used in the bottom left corner of your plot
                        col.ul,#vector of color that will be used in the upper left corner of your plot
                        col.lr,#vector of color that will be used in the bottom right corner of your plot
                        cells.x, #number of cells you want in the x axis to show your gradient
                        cells.y) #number of cells you want in the y axis to show your gradient
                                  
{ 
  n <- cells.x * cells.y
  
  cols <- rep(NA, n)
  
  z <- 1
  for (i in 1:cells.x){
    for (j in 1:cells.y){
      w1 <- (i - 1) / (cells.x - 1)
      w2 <- (j - 1) / (cells.y - 1)
      if (w1 + w2 == 0) {
        v <- col.ll
      } else { 
        v <- (w1 * (((1 - w1) * col.ll + w1 * col.lr)) + w2 * ((1 - w2) * col.ll + w2 * col.ul)) / (w1 + w2)
      }
      cols[z] <- rgb(v[1], v[2], v[3])
      z <- z + 1
    }
  }
  
  #image(matrix(c(1:n), nrow=cells.x, ncol=cells.y, byrow=TRUE), col=cols, axes=F)
  
  # then plot the coordinates on top of this grid
  
  c.y <- seq(minvaluey, maxvaluey, length=cells.y)
  c.x <- seq(minvaluex, maxvaluex, length = cells.x)
  
  cols.df<-expand.grid(c.y,c.x)
  cols.df$color<-cols
  cols.df$rastvalue<-c(1:nrow(cols.df))
  
  pb <- txtProgressBar(min = 0, max = nrow(cols.df), style = 3) #To check the progress
  colorgrid<-list()
  
  for (i in 1:nrow(cols.df))
  {
    setTxtProgressBar(pb, i)
    colorgrid[[i]]<-subset(toy_df,toy_df$FPC >=cols.df[i,"Var1"] & 
                             toy_df$FPC<=cols.df[i + 1,"Var1"] &
                             toy_df$elevation >=cols.df[i,"Var2"] & 
                             toy_df$elevation<cols.df[i + cells.y,"Var2"])
    if (nrow(colorgrid[[i]])>0)
    {
      colorgrid[[i]]$color<-cols.df[i,"color"]
      colorgrid[[i]]$rastvalue<-cols.df[i,"rastvalue"]
      
    }else{
      colorgrid[[i]]$color<-as.character()
      colorgrid[[i]]$rastvalue<-as.numeric()
    }
  }
  
  return(do.call(rbind,colorgrid))
}

#run the function to create a column assigning colours to each cell (gradient based)
toy_df<-biplot_values(toy_df,
                      -0.5, #the minimum value of the variable you want to plot on the x axis
                      100,#the minimum value of the variable you want to plot on the y axis
                      450,#the maximum value of the variable you want to plot on the x axis
                      200,#the maximum value of the variable you want to plot on the y axis
                      c(253,184,99)/255,#vector of color that will be used in the bottom left corner of your plot
                      c(94, 60, 153)/255,#vector of color that will be used in the upper left corner of your plot
                      c(230,30,1)/255,#vector of color that will be used in the bottom right corner of your plot
                      20, #number of cells you want in the x axis to show your gradient
                      20 #number of cells you want in the y axis to show your gradient
                      )

  

#if you check your dataframe, it will have a column "color" with the color codes
#that have to be assigned to each point. In this case, The column "rastervalue" 
#specifies the values that will be used to generate the raster based on the association
#between the two variables
head(toy_df)


#Let's create a raster using the continuous gradient of colours from the data.frame using the function rasterFromXYZ
association_raster_cont<-rasterFromXYZ(toy_df[c(1,2,6)])
crs(association_raster)<-"+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80
+towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#Let's create a raster using the four categories from the data.frame using the function rasterFromXYZ
association_raster_cat<-rasterFromXYZ(toy_df[c(1,2,5)])
crs(association_raster)<-"+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80
+towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#save your raster
#writeRaster(association_raster,filename = "association_raster.tif")

##################################################################
######################TIME TO CREATE SOME PLOTS###################
#################################################################

#below there are different examples of how to plot the association among variables:

#sample some cells from your new dataframe to create a scatter plot
sample_df<-toy_df[sample(nrow(toy_df),10000),]


#simple plot elevation vs FPC (to get familiar with key components of ggplot)

ggplot(sample_df, aes(x=elevation,y=FPC)) +
  geom_point(size=1, color = "red")+#check ggplot cheat sheet
  #geom_smooth(aes(group = 1),method='lm',formula=y~x, se=TRUE)+#add a trendline
  
  theme_classic()+ #white background
  theme(axis.line.x = element_line(color="black", size = 0.5), #theme allows you to customize the non-data components of your plots: i.e. titles, labels, fonts, background, gridlines, and legends
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text=element_text(size=11),
        axis.title=element_text(size=11),
        plot.title = element_text(hjust = 0.5))+ 
  scale_y_continuous(breaks = seq(100, 200, 20))+#control y axis divisions
  scale_x_continuous(breaks = seq(0, 400, 50))+#control x axis divisions
  labs(x = " Elevation",
       y = "FPC", size=14)+
  ggtitle("Association between elevation and Foliage projective cover in the city of Logan")

#create plot showing the four categories of association.
#the difference here compared to the first plot is that we are basically using "scale_color_manual"
#to assign specific colurs

p1<-ggplot(sample_df, aes(x=elevation,y=FPC)) +
  geom_point(size=1, aes(colour= as.factor(sample_df$relationship)))+
  ggplot2::scale_color_manual(values = c('#fdb863','#e66101','#b2abd2','#5e3c99'), # check the links I sent to get these codes
                              name = "Association", labels = c("LL", "LH", "HL","HH"))+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text=element_text(size=11),
        axis.title=element_text(size=11),
        plot.title = element_text(hjust = 0.5))+ 
  scale_y_continuous(breaks = seq(100, 200, 20))+
  scale_x_continuous(breaks = seq(0, 450, 50))+
  labs(x = " Elevation",
       y = "FPC", size=14)+
  ggtitle("")+
 theme(legend.position="none",legend.box = "horizontal")#note that I am excluding the legend


#let's plot the continuos scheme we created using the function "biplot_values"
p2<-ggplot(sample_df, aes(x=elevation,y=FPC)) +
  geom_point(size=1, colour= sample_df$color)+
  theme_classic()+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text=element_text(size=11),
        axis.title=element_text(size=11))+ 
  scale_y_continuous(breaks = seq(100, 200, 100))+
  scale_x_continuous(breaks = seq(0, 400, 50))+
  labs(x = " Elevation",
       y = "FPC", size=14)+
  ggtitle("")+
  theme(legend.position="none",legend.box = "horizontal")


#another option is using the plot function
plot(sample_df$elevation,sample_df$FPC,col=sample_df$color,
     xlab = "Elevation",ylab = "FPC",cex.lab=1.5,cex.axis=1,
     main="Association between elevation and FPC in the Logan City Council")
#if you want to save this plot as an object you need to use the function recordPlot
p3<-recordPlot()

#############################################################
#let's create a figure with the map and the biplot together
#############################################################

cols<-unique(toy_df$color)#get vector of colors from dataframe

# here he's using rasterVis instead of ggplot2, and adding layers just like you would for geom_point or whatever
# then lattice extra is another package which adds the scale bar and north arrows with the layer function
map_logan_cont<- rasterVis::levelplot(association_raster_cont, 
                                 margin=FALSE,                       # suppress marginal graphics
                                 colorkey=FALSE,
                                 #colorkey=list(
                                 #space='none',                   # plot legend at bottom
                                 #labels=list(at=0:4, font=4)      # legend ticks and labels 
                                 #),    
                                 #par.settings=list(
                                 #axis.line=list(col='black') # suppress axes and legend outline
                                 #),
                                 #scales=list(draw=TRUE),            # suppress axis labels
                                 col.regions=cols,                  # colour ramp
                                 at=seq(min(toy_df$rastvalue), max(toy_df$rastvalue), len=length(cols)))+
          latticeExtra::layer(sp.polygons(LCC,col= "black"))+ #add border of LCC
          latticeExtra::layer({
            SpatialPolygonsRescale(layout.north.arrow(type = 2), #ADD NORTH ARROW
                                   offset = c(2020000,-3170000),
                                   scale = 5000)
          })+
          latticeExtra::layer({ #ADD SCALE BAR
            xs <- seq(2050000, 2060000, by=2500) #x coordinates for location within the map
            grid.rect(x=xs, y=-3196500,  #y coordinates for location within the map
                      width=2500, height=400, #size of the scale bar
                      gp=gpar(fill=rep(c('transparent', 'black'), 2)),
                      default.units='native')
            grid.text(x= xs - 1000, y=-3195500, seq(0, 10000, by=2500), #text coordinates
                      gp=gpar(cex=0.5), rot=0,
                      default.units='native')
          })

#see the map:   
map_logan_cont 

#create map showing only four categories:

map_logan_cat<- rasterVis::levelplot(association_raster_cat, 
                                      margin=FALSE,                       # suppress marginal graphics
                                      colorkey=FALSE,
                                      col.regions=c('#fdb863','#e66101','#b2abd2','#5e3c99'),                  # colour ramp
                                      at=seq(0, max(toy_df$relationship), len=5))+
  latticeExtra::layer(sp.polygons(LCC,col= "black"))+
  latticeExtra::layer({
    SpatialPolygonsRescale(layout.north.arrow(type = 2),
                           offset = c(2020000,-3170000),
                           scale = 5000)
  })+
  latticeExtra::layer({
    xs <- seq(2050000, 2060000, by=2500) #x coordinates for location within the map
    grid.rect(x=xs, y=-3196500,  #y coordinates for location within the map
              width=2500, height=400, #size of the scale bar
              gp=gpar(fill=rep(c('transparent', 'black'), 2)),
              default.units='native')
    grid.text(x= xs - 1000, y=-3195500, seq(0, 10000, by=2500), #text coordinates
              gp=gpar(cex=0.5), rot=0,
              default.units='native')
  })

#see the map:   
map_logan_cat 


#now let's combine the plots into one figure using the packages gridExtra and ggpubr

#combine plot and map for continous gradient
final_figure_cont<-grid.arrange(map_logan_cont, p2, #map + plot showing the association
             ncol = 2, nrow = 1, #the number of columns and rows in your figure
             widths = c(2.2, 1.5))
#add labels
figure1<-annotate_figure(final_figure_cont,
                top = text_grob("Association between Foliage Projective Cover and temperature in the Logan City Council LGA", color = "black", face = "bold", size = 14),
                #bottom = text_grob("I'm done, thanks :-)!", color = "blue",
                                   hjust = 1, x = 1, face = "italic", size = 10,
                #left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
                fig.lab = "Figure 1", fig.lab.face = "bold"
          )
figure1

#combine plot and map for categorical values of the association
final_figure_cat<-grid.arrange(map_logan_cat, p1, #map + plot showing the association
                                ncol = 2, nrow = 1, #the number of columns and rows in your figure
                                widths = c(2.2, 1.5))

figure2<-annotate_figure(final_figure_cat,
                         top = text_grob("Association between Foliage Projective Cover and temperature in the Logan City Council LGA", color = "black", face = "bold", size = 14),
                         bottom = text_grob("I'm done, thanks :-)!", color = "blue",
                                            hjust = 1, x = 1, face = "italic", size = 10),
                         #left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
                         fig.lab = "Figure 2", fig.lab.face = "bold"
)
figure2