# Work on Wisconsin within-host data
library(rethinking)
require(reshape2)
require(ggplot2)
library(ggridges)



# read in dataframe
df = read.csv("../data/WI-variants-vs-geo-2021-02-24.csv", header=TRUE)

# remove NaNs and drop Ct column
df = subset(df, select=-c(Ct_diff))
df = na.omit(df)

# let's first plot out a few relationships
ggplot(df, aes(x=prop_snvs_shared, y=divergence)) + geom_point()
ggplot(df, aes(x=prop_snvs_shared, y=great_circle_distance_km)) + geom_point()
ggplot(df, aes(x=prop_snvs_shared, y=Ct_diff)) + geom_point()
ggplot(df, aes(x=prop_snvs_shared, y=clades_same)) + geom_point()
ggplot(df, aes(x=divergence, y=great_circle_distance_km)) + geom_point()
ggplot(df, aes(y=divergence, x=clades_same)) + geom_point()


# using the built-in GLM function, we see estimated, small coefficients for everything
model.all = glm(prop_snvs_shared~divergence+clades_same+great_circle_distance_km+household_same,data=df,family = gaussian(link="identity"))
summary(model.all)


# I want to quantify the contributions of great circle distance, genetic divergence, Ct differnces, and clade on the proportion of within-host variation that is shared among samples.

# First, let's standardize some data 
df$prop_norm = (df$prop_snvs_shared - mean(df$prop_snvs_shared))/sd(df$prop_snvs_shared)
df$divergence_norm = (df$divergence - mean(df$divergence))/sd(df$divergence)
df$Ct_norm = (df$Ct_diff - mean(df$Ct_diff))/sd(df$Ct_diff)
df$circle_norm = (df$great_circle_distance_km - mean(df$great_circle_distance_km))/sd(df$great_circle_distance_km)

df$prop_norm = df$prop_snvs_shared
df$divergence_norm = df$divergence
df$Ct_norm = df$Ct_diff

# let's start with a simple regression model, where we want to use just an intercept to predict the proportion of variants shared among pairwise sets of tips
# we first fit the simplest model, one with intercept only; this tells us that given no other information, we expect a mean shared variatns of 28%, +- 6% ; that is sort of a lot
m1 <- map(
  alist(
    prop_snvs_shared ~ dnorm(mu, sigma),
    mu <- a,
    a ~ dnorm(0,1),
    sigma ~ dunif(0,10)
  ), data=df)
precis(m1)
plot(precis(m1))

# now fit a model with just 1 predictor, difference in Cts
m2 <- map(
  alist(
    prop_snvs_shared ~ dnorm(mu, sigma),
    mu <- a + betaCt*Ct_diff,
    a ~ dnorm(0,1),
    betaCt ~ dnorm(0,1),     # we set our prior here to have an ~ 1 log diff with an sd of 3 logs
    sigma ~ dunif(0,10)
  ), data=df)
precis(m2, digits=6, prob=0.95)
plot(precis(m2))


# now fit a model with just 1 predictor, clades
m3 <- map(
  alist(
    prop_snvs_shared ~ dnorm(mu, sigma),
    mu <- a + betaClade*clades_same,
    a ~ dnorm(0,1),
    betaClade ~ dnorm(0,1),    
    sigma ~ dunif(0,1)
  ), data=df)
precis(m3, digits=6, prob=0.95)
plot(precis(m3))

# now fit a model with just 1 predictor, divergence
m4 <- map(
  alist(
    prop_snvs_shared ~ dnorm(mu, sigma),
    mu <- a + betaDiv*divergence,
    a ~ dnorm(0,1),
    betaDiv ~ dnorm(0,1),    
    sigma ~ dunif(0,1)
  ), data=df)
precis(m4, digits=6, prob=0.95)
plot(precis(m4))


# now fit a model with just 1 predictor, great circle distance
m5 <- map(
  alist(
    prop_snvs_shared ~ dnorm(mu, sigma),
    mu <- a + betaGeo*great_circle_distance_km,
    a ~ dnorm(0,1),
    betaGeo ~ dnorm(0,1),    
    sigma ~ dunif(0,1)
  ), data=df)
precis(m5, digits=6, prob=0.95)
plot(precis(m5))

# now fit a model with just 1 predictor, great circle distance
mlocsame <- map(
  alist(
    prop_snvs_shared ~ dnorm(mu, sigma),
    mu <- a + betaLoc_same*location_same,
    a ~ dnorm(0,1),
    betaLoc_same ~ dnorm(0,1),    
    sigma ~ dunif(0,1)
  ), data=df)
precis(mlocsame, digits=6, prob=0.95)
plot(precis(mlocsame))

m6 <- map(
  alist(
    prop_snvs_shared ~ dnorm(mu, sigma),
    mu <- a + betaHousehold*household_same,
    a ~ dnorm(0,1),
    betaHousehold ~ dnorm(0,1),    
    sigma ~ dunif(0,1)
  ), data=df)
precis(m6, digits=6, prob=0.95)
plot(precis(m6))

m7 <- map(
  alist(
    prop_snvs_shared ~ dnorm(mu, sigma),
    mu <- a + betaTranspair*transmission_pair,
    a ~ dnorm(0,1),
    betaTranspair ~ dnorm(0,1),    
    sigma ~ dunif(0,1)
  ), data=df)
precis(m7, digits=6, prob=0.95)
plot(precis(m7))


# now fit a model with everything in it except Ct
m8 <- map(
  alist(
    prop_snvs_shared ~ dnorm(mu, sigma),
    mu <- a + betaDiv*divergence + betaClade*clades_same + betaHousehold*household_same + betaGeo*great_circle_distance_km,
    a ~ dnorm(0,1),
    betaDiv ~ dnorm(0,1),
    betaClade ~ dnorm(0,1),
    betaGeo ~ dnorm(0,1),
    betaHousehold ~ dnorm(0,1), 
    sigma ~ dunif(0,10)
  ), data=df)
precis(m8, digits=6, prob=0.95)
plot(precis(m8))

compare(m1,m3,m4,m5,m6,m7,m8)

# ok, so from these comparisons, it looks like household pair is a better predictor than transmission pair, and that they contribute the same information in the model. So we should probably just use household pair and not transmission pair 

# plot coefficients mean estimates with 95th percentiles
precis(m8, digits=6, prob=0.95)
plot(precis(m8, prob=0.95))

# plot!
# set output of precis to dataframe; this will let us plot the mean coefficient estimates
# with the 95% percentiles
df1 = data.frame(precis(m8, digits=6, prob=0.95))

# reindex dataframe
df2 <- cbind(measurement = rownames(df1), df1)
rownames(df2) <- 1:nrow(df2)
colnames(df2) = c("measurement","mean","stdev","lower95","upper95")

df2$measurement = gsub("betaHousehold","household",df2$measurement)
df2$measurement = gsub("betaGeo","distance",df2$measurement)
df2$measurement = gsub("betaDiv","divergence",df2$measurement)
df2$measurement = gsub("betaClade","clade",df2$measurement)
df2$measurement = gsub("sigma","variance",df2$measurement)

# rename a to intercept; format here is dataframe[rownumber, columnnumber] = new value
df2[1, 1] = "intercept"

df2$measurement_f = factor(df2$measurement, levels=c("intercept","variance","distance","divergence","clade","household"))


purple = "#5248AA"
blue = "#7C9BAC"
yellow = "#ECC58C"
red = "#551E32"
green = "#434A42"
brown = "#9E6240"
pink =  "#B4656F"
grey = "#939393"

p <- ggplot(data=df2, aes(y=measurement_f, x=mean, color=measurement_f)) + 
  geom_pointrange(data=df2, aes(xmin=lower95, xmax=upper95), size=1, alpha=0.7)+
  labs(x="\nmean",y="coefficient\n")+
  scale_color_manual(values=c(intercept=grey,variance=grey,distance=red,divergence=yellow,clade=blue,household=purple), guide=FALSE)+
  #scale_y_continuous(breaks=seq(0,0.8,0.2), limits=c(0,0.8))+
  scale_x_continuous(breaks=seq(0,0.3,0.1), limits=c(-0.05,0.3))+
  theme(panel.grid.major=element_line(colour=NA,size=NA))+    
  theme(panel.grid.minor=element_line(colour=NA,size=NA))+    
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(axis.line.x=element_line(colour="black"))+
  theme(axis.line.y=element_line(colour="black"))+
  theme(axis.title=element_text(size=24, vjust=5))+
  theme(axis.text.y=element_text(size=20, colour="black"))+
  theme(axis.text.x=element_text(size=20, colour="black", hjust=0.5))+
  theme(legend.text=element_text(size=20))+
  theme(legend.title=element_blank())+
  #theme(panel.margin=unit(1, "lines"))+
  #theme(plot.margin=unit(c(1,4,1,1),"cm"))+
  theme(legend.key.size=unit(0.7, "cm"))+
  theme(panel.background=element_rect(fill=NA))+
  theme(legend.key=element_rect(fill=NA))

ggsave("regression_coefficients-2020-02-25.pdf", width = 6, height = 4, device = "pdf", path = "../figures/individual-pdfs", dpi = 300)
p 


# now try plotting the same figure but plotting the posteriors instead
post_m8 <- extract.samples(m8)
colnames(post_m8) = c("intercept","divergence","clade","distance","household","variance")

# gather the 95% HPD estimate for each variable
hpdi_int = HPDI(post_m8$intercept, prob=0.95)
hpdi_var = HPDI(post_m8$variance, prob=0.95)
hpdi_div = HPDI(post_m8$divergence, prob=0.95)
hpdi_clade = HPDI(post_m8$clade, prob=0.95)
hpdi_dist = HPDI(post_m8$distance, prob=0.95)
hpdi_hh = HPDI(post_m8$household, prob=0.95)

# separate into high and low bounds of HPDI
lower_hpdi_int = hpdi_int[[1]]
upper_hpdi_int = hpdi_int[[2]]
lower_hpdi_var = hpdi_var[[1]]
upper_hpdi_var = hpdi_var[[2]]
lower_hpdi_div = hpdi_div[[1]]
upper_hpdi_div = hpdi_div[[2]]
lower_hpdi_clade = hpdi_clade[[1]]
upper_hpdi_clade = hpdi_clade[[2]]
lower_hpdi_dist = hpdi_dist[[1]]
upper_hpdi_dist = hpdi_dist[[2]]
lower_hpdi_hh = hpdi_hh[[1]]
upper_hpdi_hh = hpdi_hh[[2]]


# now, melt it 
post_m8_melt = melt(post_m8)
post_m8_melt$variable_f = factor(post_m8_melt$variable, levels=c("divergence","clade","household","distance","variance","intercept"))

# subest melted dataframe to include only vvalues in the HPDI
post_m8_melt_hpdi <- post_m8_melt[!(post_m8_melt$variable=="intercept" & post_m8_melt$value > upper_hpdi_int),]
post_m8_melt_hpdi <- post_m8_melt_hpdi[!(post_m8_melt_hpdi$variable=="intercept" & post_m8_melt_hpdi$value < lower_hpdi_int),]

post_m8_melt_hpdi <- post_m8_melt_hpdi[!(post_m8_melt_hpdi$variable=="variance" & post_m8_melt_hpdi$value > upper_hpdi_var),]
post_m8_melt_hpdi <- post_m8_melt_hpdi[!(post_m8_melt_hpdi$variable=="variance" & post_m8_melt_hpdi$value < lower_hpdi_var),]

post_m8_melt_hpdi <- post_m8_melt_hpdi[!(post_m8_melt_hpdi$variable=="clade" & post_m8_melt_hpdi$value > upper_hpdi_clade),]
post_m8_melt_hpdi <- post_m8_melt_hpdi[!(post_m8_melt_hpdi$variable=="clade" & post_m8_melt_hpdi$value < lower_hpdi_clade),]

post_m8_melt_hpdi <- post_m8_melt_hpdi[!(post_m8_melt_hpdi$variable=="distance" & post_m8_melt_hpdi$value > upper_hpdi_dist),]
post_m8_melt_hpdi <- post_m8_melt_hpdi[!(post_m8_melt_hpdi$variable=="distance" & post_m8_melt_hpdi$value < lower_hpdi_dist),]

post_m8_melt_hpdi <- post_m8_melt_hpdi[!(post_m8_melt_hpdi$variable=="divergence" & post_m8_melt_hpdi$value > upper_hpdi_div),]
post_m8_melt_hpdi <- post_m8_melt_hpdi[!(post_m8_melt_hpdi$variable=="divergence" & post_m8_melt_hpdi$value < lower_hpdi_div),]

post_m8_melt_hpdi <- post_m8_melt_hpdi[!(post_m8_melt_hpdi$variable=="household" & post_m8_melt_hpdi$value > upper_hpdi_hh),]
post_m8_melt_hpdi <- post_m8_melt_hpdi[!(post_m8_melt_hpdi$variable=="household" & post_m8_melt_hpdi$value < lower_hpdi_hh),]


p <- ggplot(data=post_m8_melt_hpdi, aes(x=value, y=variable_f, color=variable_f, fill=variable_f)) + 
  geom_density_ridges(alpha=0.6, scale=1, rel_min_height=0.01, size=0.25, position=position_points_jitter(width=0.05,height=0), pointsize=0.5, point_alpha = 0.5) +
  labs(x="\nposterior",y="coefficient\n")+
  scale_color_manual(values=c(intercept=grey,variance=grey,distance=red,divergence=yellow,clade=blue,household=purple), guide=FALSE)+
  scale_fill_manual(values=c(intercept=grey,variance=grey,distance=red,divergence=yellow,clade=blue,household=purple), guide=FALSE)+
  #scale_y_continuous(breaks=seq(0,10000,1000), limits=c(0,10000))+
  #scale_x_continuous(breaks=seq(0,0.3,0.1), limits=c(-0.05,0.3))+
  theme(panel.grid.major.y=element_line(colour="black",size=0.1))+    
  theme(panel.grid.major.x=element_line(colour=NA,size=NA))+
  #theme(panel.grid.minor=element_line(colour=NA,size=NA))+    
  #theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(axis.line.x=element_line(colour="black"))+
  theme(axis.line.y=element_line(colour="black"))+
  theme(axis.title=element_text(size=12, vjust=5))+
  theme(axis.text.y=element_text(size=10, colour="black"))+
  theme(axis.text.x=element_text(size=10, colour="black", hjust=0.5))+
  theme(legend.title=element_blank())+
  #theme(panel.margin=unit(1, "lines"))+
  #theme(plot.margin=unit(c(1,4,1,1),"cm"))+
  theme(legend.key.size=unit(0.7, "cm"))+
  theme(panel.background=element_rect(fill=NA))+
  theme(legend.key=element_rect(fill=NA))

ggsave("regression_coefficients-posteriors-2020-03-25.png", width =2.75, height = 2, path = "src/ncov-WI-within-host/figures/individual-pdfs/", dpi = 300)
p 



# Plot the full posterior distribution facetted, with alpha = HPDI

# add a column to store whether the value is within the HPDI; initialize with "yes"
post_m8_melt$hpdi <- "yes"

# add in "no"s for values outside HPDI
post_m8_melt$hpdi[which(post_m8_melt$variable == "intercept" & post_m8_melt$value < lower_hpdi_int)] = "no"
post_m8_melt$hpdi[which(post_m8_melt$variable == "intercept" & post_m8_melt$value > upper_hpdi_int)] = "no"

post_m8_melt$hpdi[which(post_m8_melt$variable == "variance" & post_m8_melt$value < lower_hpdi_var)] = "no"
post_m8_melt$hpdi[which(post_m8_melt$variable == "variance" & post_m8_melt$value > upper_hpdi_var)] = "no"

post_m8_melt$hpdi[which(post_m8_melt$variable == "clade" & post_m8_melt$value < lower_hpdi_clade)] = "no"
post_m8_melt$hpdi[which(post_m8_melt$variable == "clade" & post_m8_melt$value > upper_hpdi_clade)] = "no"

post_m8_melt$hpdi[which(post_m8_melt$variable == "divergence" & post_m8_melt$value < lower_hpdi_div)] = "no"
post_m8_melt$hpdi[which(post_m8_melt$variable == "divergence" & post_m8_melt$value > upper_hpdi_div)] = "no"

post_m8_melt$hpdi[which(post_m8_melt$variable == "distance" & post_m8_melt$value < lower_hpdi_dist)] = "no"
post_m8_melt$hpdi[which(post_m8_melt$variable == "distance" & post_m8_melt$value > upper_hpdi_dist)] = "no"

post_m8_melt$hpdi[which(post_m8_melt$variable == "household" & post_m8_melt$value < lower_hpdi_hh)] = "no"
post_m8_melt$hpdi[which(post_m8_melt$variable == "household" & post_m8_melt$value > upper_hpdi_hh)] = "no"


# plot it. I'm going to have to do some sort of fancy things to give each posterior its own x-axis.

# define the axis limits; this is calculated based on the min and max values in the posterior set, with a little bit of extra padding

require(grid)
library(grid)
require(gridExtra)
library(gridExtra)

post_m8_melt$variable_f = factor(post_m8_melt$variable, levels=c("divergence","clade","household","distance","variance","intercept"))

blank_data <- data.frame(variable_f = c("divergence","divergence","distance","distance","clade","clade","household","household","variance","variance","intercept","intercept"), x = c(-0.002,0.002,-0.0003,-0.0001,0.05,0.10,0.2,0.40,0.15,0.17,0.025,0.1), hpdi="yes")

p_facet <- ggplot(data=post_m8_melt, aes(x=value, color=variable_f, fill=variable_f, group=variable_f)) + 
  geom_density(alpha=0.8, size=0) +
  facet_wrap(~variable_f, scales="free")+
  geom_blank(data=blank_data, aes(x=x))+
  labs(x="\nposterior estimate",y="density\n")+
  scale_color_manual(values=c(intercept=grey,variance=grey,distance=red,divergence=yellow,clade=blue,household=purple), guide=FALSE)+
  scale_fill_manual(values=c(intercept=grey,variance=grey,distance=red,divergence=yellow,clade=blue,household=purple), guide=FALSE)+
  #scale_y_continuous(breaks=seq(0,10000,1000), limits=c(0,10000))+
  #scale_x_continuous(breaks=seq(seqs_list[[variable_f]]))+
  theme(panel.grid.major=element_line(colour=NA,size=NA))+    
  theme(panel.grid.minor=element_line(colour=NA,size=NA))+
  theme(panel.spacing.x=unit(4, "lines"))+
  theme(panel.spacing.y=unit(3, "lines"))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text=element_text(size=20))+
  theme(axis.line.x=element_line(colour="black"))+
  theme(axis.line.y=element_line(colour="black"))+
  theme(axis.title=element_text(size=20, vjust=5))+
  theme(axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank())+
  theme(axis.text.x=element_text(size=16, colour="black", hjust=0.5))+
  theme(legend.title=element_blank())+
  #theme(panel.margin=unit(1, "lines"))+
  #theme(plot.margin=unit(c(1,4,1,1),"cm"))+
  theme(legend.key.size=unit(0.7, "cm"))+
  theme(panel.background=element_rect(fill=NA))+
  theme(legend.key=element_rect(fill=NA))

ggsave("regression_coefficients-posteriors-facetted-2020-03-25.pdf", width = 12, height = 6, path = "src/ncov-WI-within-host/figures/individual-pdfs/", dpi = 300)
p_facet 

