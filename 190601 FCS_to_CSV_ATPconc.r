# FCS to CSV
    # Coverting .fcs file data into an .csv file with creating graphs
    # Hiroshi Kobayashi
    # 2018-08-14
    #
#  
#
##### USER INPUT #####
    

library('flowCore')
library(ggplot2)
library(viridis)
library(MASS)
library("mgcv")


PrimaryDirectory = getwd()

#ファイル中の+文字をすべて_に変更
files <- list.files(pattern="\\.fcs")

i<-1

newfiles <- gsub("\\+", "_", files)


#tmpfiles <- gsub("\\+", "_", files)
#for (file in tmpfiles){
#file <- gsub("Specimen_001_",paste("Specimen_002_",i,".",sep=""),file)
#if(i==1){ newfiles<- file}
#else { newfiles <-c(newfiles,file) }

#i <- i+1
#}

file.rename(files,newfiles)




#fractionの名前
#fraction_list <-  c("P12","P13","Q4")


# conditionの抽出
newfiles <- gsub("Specimen_001_","",newfiles)
newfiles <- gsub(".fcs","",newfiles)
condition_list <- levels(factor(newfiles))




# temporal order of each condition
file_order = list(c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1),c(1))
#file_order = list(c(2,1))

# for density estimation
grid_size= 100

get_density <- function(x, y, n) {
                  dens <- kde2d(x = x, y = y,n = n)
                  ix <- findInterval(x, dens$x)
                  iy <- findInterval(y, dens$y)
                  ii <- cbind(ix, iy)
                  return(dens$z[ii])
                }


get_diff <- function(time, FRET){
		time_0 <- c(time,0)
		time_1 <- c(0,time)
		FRET_0 <- c(FRET,0)
		FRET_1 <- c(0,FRET)

	diff_FRET <- FRET_0-FRET_1
	diff_time <- time_0-time_1
	diff <- diff_FRET/diff_time

	diff <- diff[-1]
	diff <- diff[-length(diff)]
	diff <- c(0,diff)
 return(diff)
}


l <- length(condition_list)
#For debug
#l <-1

for(condition_num in 1:l){
	condition = condition_list[condition_num]
    
# Find file names of .fcs files in the current working directory
    FileNames <- list.files(path=PrimaryDirectory, pattern = paste("Specimen_001_",condition,".fcs",sep=""))

# change the order of File names
  FileNames <- FileNames[unlist(file_order[condition_num])]
 

# data_tmp <- exprs(read.FCS("Specimen_001_#1_Glc+_Oligo+_Q4.fcs", transformation = FALSE))


number_of_Samples <- length(FileNames)
#試薬を加えた際の時間シフト（前の最後の検出時間+5000）
time_shift <- rep(0,number_of_Samples)
i=1
   for(sample in FileNames){
    
    data_tmp <- exprs(read.FCS(sample, transformation = FALSE))


#デバッグ用　生データCSV出力
#   write.csv(data_tmp,paste(condition,fraction,".csv",collapse=""))



    data_tmp_FRET <- data_tmp[,2]/data_tmp[,1]
    


    if (i == 1){
    data_tmp_TIME <- data_tmp[,3]
    time_shift[i] <- data_tmp[length(data_tmp_FRET),3] +5000
    }
    else{
    data_tmp_TIME <- data_tmp[,3] + time_shift[i-1]
    time_shift[i] <- data_tmp[length(data_tmp_FRET),3] +5000 +time_shift[i-1]
    
    }
    


    data_tmp_FRET <- cbind(data_tmp_TIME,data_tmp_FRET)
    #ATP conc high cellはのぞく

########################
####ATP conc conversion to avoid infinite, FRET >1.4 and FRET <7 is required
#########################


    keep = (data_tmp_FRET[,2] < 7) & (data_tmp_FRET[,2]>1.4)
    data_tmp_FRET <- data_tmp_FRET[keep,]
    
	a<- (data_tmp_FRET[,2]-1.4)/6
        data_tmp_FRET[,2]  <- exp((log(a/(1-a))-0.5187)/3.1234)
    
    l <- length(data_tmp_FRET)


    if(l == 1){next}
    
    if (i ==1) {
       data_FRET <- data_tmp_FRET
	i <- 2
	next
 }
    else{
    data_FRET <- rbind(data_FRET,data_tmp_FRET)
}

 i <- i+1

}

data_FRET <- data.frame(data_FRET)
#loessMod25 <- loess(data_tmp_FRET ~ data_tmp_TIME, data=data_FRET, span=0.25)
#predicted <- predict(loessMod25,newdata=data_FRET,se=TRUE)

# 25% smoothing span

#data_FRET <- data.frame(data_FRET,predicted$fit,predicted$se.fit)



#data_FRET$density <-kde(x=data_FRET[,1],y=data_FRET[,2])
density <-  get_density (x=data_FRET[,1],y=data_FRET[,2],n=grid_size)




ggplot(data_FRET)+
geom_point(size=0.5,aes(x=data_tmp_TIME, y=data_tmp_FRET,colour=density))+
stat_smooth(aes(x=data_tmp_TIME, y=data_tmp_FRET),colour="red",span = 0.3)+
scale_colour_viridis()+
theme_bw()+

theme(panel.border = element_rect(size=1,colour="black"),
      axis.title = element_text(size=18),
      axis.text.y = element_text(size=18),
      axis.text.x = element_text(size=14),
      axis.ticks.length = unit(0.5,"cm"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "none")+


xlab("Time")+
ylab("ATP")+
ylim(0,1.2)+

ggsave(paste(condition,"_ATP.png",collapse=""),dpi=600,width=max(data_tmp_TIME)/20000, height =3)


model.gam <- gam(data_tmp_FRET ~ s(data_tmp_TIME), data=data_FRET)
smoothed.gam <- predict(model.gam,se.fit = TRUE)
data_tmp_TIME <- unique(data_FRET[,1])
smoothed.fit <- smoothed.gam[[1]][!duplicated(smoothed.gam[[1]])]
smoothed.se <- smoothed.gam[[2]][!duplicated(smoothed.gam[[1]])]

up_CI  <- round(smoothed.fit + (2 * smoothed.se), digits=4)
lo_CI  <- round(smoothed.fit - (2 * smoothed.se), digits=4)
diff_1 <- get_diff(data_tmp_TIME,smoothed.fit)
diff_2 <- get_diff(data_tmp_TIME,diff_1)



csv_output= data.frame(round(data_tmp_TIME,digits=1),round(smoothed.fit,digits=4),up_CI,lo_CI,diff_1,diff_2)
names(csv_output) <- c("Time","ATP","Up_95%CI","Lo_95%CI","Diff1","Diff2")
write.csv(csv_output,paste(condition,"ATP_gam_diff.csv",collapse=""))







        }



