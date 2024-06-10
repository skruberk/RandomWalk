#Cleans output data from randomwalks where histograms of steps to hit are generated, for inpuut into fitting 
#for an exponential 

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
View(tidy)
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
write.csv(tidy2,file='SingleTermSlopeOutputTidy.csv')
data<-tidy2
tidy<-separate(data,beta,c('col1', 'col2','discs','radius','steps','runs','spring'),' ')
tidy <- na.omit(tidy)

calcfit <- function(data) {
  model <- lm(exponent ~ discs, data = data)
  summary(model)$r.squared
  #slope <- coef(model)["discs"]
}

# Calculate R-squared values for each group
linfit <- springbig %>%
  group_by(radius) %>%
  summarize(r_squared = calcfit(pick(exponent,discs))) %>%
  ungroup()

View(linfit)
# Convert radius to factor for consistent plotting
df$radius <- as.factor(df$radius)
linfit$radius <- as.factor(linfit$radius)
View(data)

springbig <- data[data$spring == 1, ]
springsmall<-data[data$spring == 0.1, ]
View(springsmall)
#plotting section 
p<- ggplot(springbig,aes(x=discs, y=exponent,color=radius))+ geom_jitter()+ geom_line(aes(group = radius))+
  geom_smooth(aes(group = radius),method = "lm", se = FALSE) +
  scale_x_continuous(breaks=seq(0,30,by=5))+
  scale_color_viridis(option = "H")+
  ylim(-0.75,0)+
  #geom_text(data = linfit, aes(label = sprintf("R-squared: %.2f", r_squared)))#, x = -Inf, y = -Inf, hjust = 0, vjust = 1)
  #scale_color_manual(values = c("red", "blue", "green"))
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p 

#clean data for filenames
data$filename<-gsub("%","",data$filename)#remove % from matlab
View(data)
data$filename <- paste0(data$filename, ".csv")
data <- data %>%
  mutate(filename = paste0('"', filename, '"')). #enclose in quotes

write.csv(data,file='filename.csv')
