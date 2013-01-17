
#R-script for HTM
setwd("/Users/stevenmartell1/Documents/IPHC/IPHC-project/HTM/")

require(Riscam)
require(ggplot2)


A <- read.admb("HTM")

A$area <- c("2A","2B","2C","3A","3B","4A","4B","4CDE")

plot(A$iyr,A$St)



# Biomass plot
bt 				<- data.frame(A$iyrs,A$bt)
ct 				<- data.frame(A$iyr,A$obs_ct)
colnames(bt)	<- c("year",A$area)
colnames(ct)	<- c("year",A$area)
mdf 			<- melt(bt,id.var="year")
mdi 			<- melt(ct,id.var="year")
p.bt  			<- ggplot(mdf,aes(x=year,y=value)) + geom_line()
p.bt 			<- p.bt + geom_point(data=mdi,aes(x=year,y=value),shape="-",col=4)
p.bt 			<- p.bt + facet_wrap(~variable,scales="fixed")
p.bt   			<- p.bt + labs(x="Year",y="6+ Biomass & removals (Mlb)")


# Fishing mortality rate
ft 				<- data.frame(A$iyr,A$ft)
colnames(ft)	<- c("year",A$area)
mdf 			<- melt(ft,id.var="year")
p.ft  			<- ggplot(mdf,aes(x=year,y=value)) + geom_line()
p.ft 			<- p.ft + facet_wrap(~variable,scales="fixed")
p.ft 			<- p.ft + labs(x="Year",y="Fishing mortality rate")


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
p.yt 			<- p.yt + labs(x="Year",y="Survey WPUE (lb per skate)")

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
p.it 			<- p.it + labs(x="Year",y="Fishery WPUE (lb per skate)")

# Average Weight in the catch
wbar 				<- data.frame(A$iyr,A$wt)
wt  			<- data.frame(A$fishery_wbar)
wt[wt<0] 		<- NA
colnames(wbar) 	<- c("year",A$area)
colnames(wt) 	<- c("year",A$area)
mdf 			<- melt(wbar,id.var="year")
mdi 			<- melt(wt,id.var="year",variable_name="WPUE")
p.wt 			<- ggplot(mdf,aes(x=year,y=value)) + geom_line() 
p.wt 			<- p.wt + geom_point(data=mdi,aes(x=year,y=value))
p.wt 			<- p.wt + facet_wrap(~variable,scales="fixed")
p.wt 			<- p.wt + labs(x="Year",y="Catch average weight (lbs)")

# Stock recruitment plot
ii   <- 1:(length(A$St)-A$agek)
S    <- A$St[ii]
R    <- A$Rt[ii+A$agek]
ss   <- seq(0,max(S),length=100)
rr   <- A$so*ss/(1+A$beta*ss)
p.sr <- qplot(S,R)+ylim(c(0,max(R)))+xlim(c(0,max(S)))
p.sr <- p.sr+geom_line(aes(x=ss,y=rr))
p.sr <- p.sr+labs(x="Spawning biomass")
p.sr <- p.sr+labs(y=paste("Recruitment (age ",A$agek,")",sep=""))


