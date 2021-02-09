# Work on Wisconsin within-host data
library(rethinking)
require(reshape2)
require(ggplot2)


# read in dataframe
df = read.csv("../data/WI-variants-vs-geo-2021-02-08.csv", header=TRUE)

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

# plot coefficients for model chosen
precis(m8, digits=6, prob=0.95)
plot(precis(m8, prob=0.95))

# plot!
# set output of precis to dataframe
df1 = data.frame(precis(m8, digits=6, prob=0.95))

# reindex dataframe
df2 <- cbind(measurement = rownames(df1), df1)
rownames(df2) <- 1:nrow(df2)
colnames(df2) = c("measurement","mean","stdev","lower95","upper95")

post_m8 <- extract.samples(m8)
post_m8 <- subset (post_m8, select = -betaGeo)

post_m8_melt = melt(post_m8)

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
  scale_x_continuous(breaks=seq(0,0.4,0.1), limits=c(-0.05,0.4))+
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

ggsave("regression_coefficients.pdf", width = 6, height = 4, device = "pdf", path = "../figures/individual-pdfs", dpi = 300)
p 

# UP TO HERE I AM UP TO DATE 
# extract samples from model m6, which has just intercept and Marshallese status
# use these samples to calculate the difference in having descendants between Marshallese and non-Marshallese tips 
# extract samples from the model fitting in m2 and m3; these samples form the posterior distributions for our variables we estimated

plot(df$prop_snvs_shared, df$household_same)


# to describe the distribution of predicted probability of having descendants for Marshallese and non-Marshallese tips, we need to apply the logistic function to the linear model that generated the data; Marshallese is coded as a 1, so it will be calculated using a sample from the posterior for the intercept + a sample from the posterior for the slope of Marshallese; for non-Marshallese, this is coded as a 0, so bm drops out
prop_shared.household_same.pred <- post_m8$a + post_m8$betaHousehold +  post_m8$betaDiv+ post_m8$betaClade
prop_shared.household_diff.pred <- post_m8$a + post_m8$betaDiv + post_m8$betaClade


df_households = melt(data.frame(prop_shared.household_diff.pred, prop_shared.household_same.pred))

ggplot(df_households)+ geom_density(aes(x=value, color=variable))+scale_color_manual(values=c("red","blue"))


# now do some implied predictions
# dummy data 
d.pred <- data.frame(
  household_same = c(0,1)  # no, yes, yes, yes, yes, no
)

mu <- link(m6, data=data.frame(household_same=d.pred))
mu.mean <- apply(mu,2,mean)
mu.HPDI <- apply(mu,2,HPDI, prob=0.89)

# summarize
pred.p <- apply(plants.ensemble$link, 2, mean)
pred.p.PI <- apply(plants.ensemble$link, 2, PI)