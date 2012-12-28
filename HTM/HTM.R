
#R-script for HTM
setwd("/Users/stevenmartell1/Documents/IPHC/IPHC-project/HTM/")

require(Riscam)
require(ggplot2)


A <- read.admb("HTM")

A$area <- c("2A","2B","2C","3A","3B","4A","4B","4CDE")

plot(A$iyr,A$St)



# Biomass plot
bt 				<- data.frame(A$iyrs,A$bt)
colnames(bt)	<- c("year",A$area)
mdf 			<- melt(bt,id.var="year")
p.bt  			<- ggplot(mdf,aes(x=year,y=value)) + geom_line()
p.bt 			<- p.bt + facet_wrap(~variable,scales="fixed")


# Survey WPUE
yt 				<- data.frame(A$iyr,A$yt)
it  			<- data.frame(A$survey_wpue)
it[it<0] 		<- NA
colnames(yt) 	<- c("year",A$area)
colnames(it) 	<- c("year",A$area)
mdf 			<- melt(yt,id.var="year")
mdi 			<- melt(it,id.var="year",variable_name="WPUE")
p.yt 			<- ggplot(mdf,aes(x=year,y=value)) + geom_line() 
p.yt 			<- p.yt + geom_point(data=mdi,aes(x=year,y=value))
p.yt 			<- p.yt + facet_wrap(~variable,scales="fixed")


# Fishery WPUE
yt 				<- data.frame(A$iyr,A$it)
it  			<- data.frame(A$fishery_wpue)
it[it<0] 		<- NA
colnames(yt) 	<- c("year",A$area)
colnames(it) 	<- c("year",A$area)
mdf 			<- melt(yt,id.var="year")
mdi 			<- melt(it,id.var="year",variable_name="WPUE")
p.it 			<- ggplot(mdf,aes(x=year,y=value)) + geom_line() 
p.it 			<- p.it + geom_point(data=mdi,aes(x=year,y=value))
p.it 			<- p.it + facet_wrap(~variable,scales="fixed")



