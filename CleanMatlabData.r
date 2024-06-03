library(dplyr)
library(tidyr)
library(janitor)
library(readxl)
library(writexl)
library(ggplot2)
library(ggbeeswarm)
library(viridis)

#data <- remove_empty(data,c("rows"))
tidy <- na.omit(data)
View(tidy)
#tidy$filename <- strsplit(tidy$filename, ' ')
tidy<-separate(data,filename,c('col1', 'col2','discs','radius','steps','runs','spring'),' ')
tidy <- na.omit(tidy)
tidy <- tidy[ -c(1,2) ]#remove first two columns 
tidy2<-tidy
View(tidy2)
tidy2$spring<-gsub(".csv"," ",as.character(tidy2$spring)) #remove file extension
View(tidy2)
tidy<-tidy2
View(tidy)
tidy2<-tidy2 %>% 
  mutate(
  spring = as.numeric(format(tidy2$spring, scientific = FALSE)),
  radius = as.numeric(format(tidy2$radius, scientific = FALSE))
)#remove sci not
View(tidy2)
write.csv(tidy2,file='SlopeOutputTidy3.csv')

springbig <- data[data$spring == 1, ]
springsmall<-data[data$spring == 0.1, ]

View(springsmall)
#plotting section 
p<- ggplot(springsmall,aes(x=discs, y=slope1,color=radius))+ geom_jitter()+ geom_line(aes(group = radius))+
  #(aes(group = radius),method = "lm", se = FALSE) +
  scale_color_viridis(option = "H")+
  ylim(-100,0)+
  #scale_color_manual(values = c("red", "blue", "green"))
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p 
p<-p + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p


