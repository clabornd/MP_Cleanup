source("requirements.R")

Rad.full <- read_csv("Rad_full_604.csv",col_types = cols(X1 = col_skip()))
Rad.full <- Rad.full[-c(9222, 9530, 9571, 9652, 9693, 9734, 9775, 9817, 9859, 9901),] 

###################################################
############## DATAFRAME CONSTRUCTION #############
###################################################

alpha <- .05

Rad.clean.small <- Rad.clean.small %>% group_by(Lat.dec, Long.dec) %>% summarise_each(funs(mean), Value)

#Method 1: Get "saturated" spatial points, calculate correlation coefficient, perform kriging on correlation coefficient.
locs <- unique(cbind(Rad.clean$Lat.dec, Rad.clean$Long.dec))

Rad.full <- data.frame()
#count <- 0

for(i in 1:dim(unique(cbind(Rad.clean$Lat.dec, Rad.clean$Long.dec)))[1]){
  if(length(which(Rad.clean$Lat.dec == locs[i,1] & Rad.clean$Long.dec == locs[i,2])) >= 4){
    #count <- count + 1
    Rad.full <- rbind(Rad.full, Rad.clean[which(Rad.clean$Lat.dec == locs[i,1] & Rad.clean$Long.dec == locs[i,2]),])
  }
}

###MAIN DATAFRAME, CONTAINS ALL COEFFICIENTS AND P-VALUES###
Rad.full <- Rad.full %>% slice_rows(c("Lat.dec", "Long.dec")) %>% 
  mutate(rank.val = rank(Value, ties.method = "random"),
         rank.time = rank(time.int),
         groupsize = n()) %>%
  unslice()%>%
  arrange(groupsize, Lat.dec)%>%
  group_by(groupsize, Lat.dec, Long.dec) %>%
  mutate(rho = cor(rank.time, rank.val), 
         pearson = cor(time.int, Value), 
         p_spearman = cor.test(rank.time, rank.val, method = "s",continuity = TRUE, alternative = "less")$p.value,
         p_pearson = cor.test(time.int, Value, method = "p",continuity = TRUE, alternative = "less")$p.value,
         CI_p_L = rep(cor.test(time.int, Value, method = "p",continuity = TRUE, alternative = "less")$conf.int[1],n()),
         CI_p_U = rep(cor.test(time.int, Value, method = "p",continuity = TRUE, alternative = "less")$conf.int[2],n())
  ) %>%
  mutate(CI_s_L = ((1+rho)/(1-rho)*exp(-2*qnorm(0.95,0,1)*sqrt(1/(n()-3)))-1)/((1+rho)/(1-rho)*exp(-2*qnorm(0.95,0,1)*sqrt(1/(n()-3)))+1),
         CI_s_U = ((1+rho)/(1-rho)*exp(2*qnorm(0.95,0,1)*sqrt(1/(n()-3)))-1)/((1+rho)/(1-rho)*exp(2*qnorm(0.95,0,1)*sqrt(1/(n()-3)))+1)
  ) %>%
  ungroup()

Rad.full["ID"] <- as.numeric(row.names(Rad.full))

#vector of group sizes for bayesian model
groupsizes <- Rad.full %>% group_by(Lat.dec, Long.dec) %>% summarise_each(funs(mean), groupsize) %>%
  ungroup() %>% dplyr::select(groupsize) %>% arrange(groupsize)

write_csv(Rad.full, "Rad_full_604.csv")

#collapse
Rad.ranks <- Rad.full %>% group_by(groupsize, Lat.dec, Long.dec) %>%
  summarise_each(funs(mean), rho, Lat.dec, Long.dec, 
                 dist, pearson, groupsize,p_spearman,p_pearson,
                 CI_p_L, CI_p_U, CI_s_U, CI_s_L) %>%
  ungroup()
initval <- Rad.full %>% group_by(Lat.dec, Long.dec) %>% slice(1) %>% ungroup() %>% dplyr::select(Value)

Rad.ranks["ID"] <- as.numeric(row.names(Rad.ranks))
Rad.ranks["init"] <- initval

#Random subset for bayesian model
randlocs <- sample_n(Rad.ranks[c("Lat.dec", "Long.dec")], 50)
testpoints <- data.frame()

for(i in 1:dim(Rad.full)[1]){
  for(j in 1:50){
    if(all(Rad.full[i,c("Lat.dec", "Long.dec")]==randlocs[j,])){
      testpoints <- rbind(testpoints, as.vector(Rad.full[i,]))
    }
  }
}

testpoints <- testpoints %>% group_by(Lat.dec, Long.dec) %>% 
  mutate(rank.val = rank(Value), rank.time = rank(time.int)) %>% ungroup()

write.csv(testpoints, "Rad_testpoints.csv")

################################################

##################
####REGRESSION####
##################

lmod <- lm(Value~dist*time.int, data = Rad.full)
summary(lmod)
qplot(fitted(lmod), residuals(lmod))

#nonlinear model
nlmod <- nls(data = Rad.full, control = nls.control(maxiter = 1000),
             log(Value) ~ z + 1/dist^a+b*time.int, 
             start = list(a = 2, z = 1, b=1))
summary(nlmod)
qplot((fitted(nlmod)), residuals(nlmod))

#mixed effect for time
mixed.log <- lme(fixed = log(Value) ~ dist*time.int, data = Rad.full, random = ~1|time.int)
vario.mixed <- Variogram(mixed.log, maxDist = 0.25)
mixed.geo <- as.geodata(cbind(Rad.full$Lat.dec, Rad.full$Long.dec, residuals(mixed.log)))
mixed.v <- variog(mixed.geo, max.dist = 0.25)
mixed.Exp <- variofit(mixed.v, fix.nugget = FALSE, cov.model='exponential')

plot(mixed.v)
lines(mixed.Exp)

rad.mmod <- lme(fixed = log(Value) ~ dist*time.int, data = Rad.full, random = ~1|time.int, method = "ML")
summary(rad.mmod)

rad.mmod.log <- lme(fixed = log(Value) ~ dist*time.int, data = Rad.full, 
                    random = ~1|time.int, 
                    correlation = corExp(0.1347, 
                    form = ~Lat.dec+Long.dec|time.int, fixed = TRUE), 
                    method = "ML"
                    )
vario.mixed <- Variogram(rad.mmod.log)


#nls for rho

gls.rho <-  gls(rho~dist*init, data = Rad.ranks, method = "ML")
summary(gls.rho)
qplot((fitted(gls.rho)), residuals(gls.rho))

gls.corr.rho <- update(gls.rho, correlation = corExp(1, form = ~Lat.dec+Long.dec), method = "ML")

Rad.ranks["resids"] <- residuals(gls.rho)

gls.corr <- update(gls.rho, correlation = corExp(1, form = ~Lat.dec + Long.dec), method = "ML")

gnls.rho = gnls(rho ~ z+1/dist^a+b*init, data = Rad.ranks, start = list(z=1, a=2, b=-1))
summary(gnls.rho)

gnls.corr <- update(gnls.rho, correlation = corExp(1, form = ~Lat.dec+Long.dec))

mod.pearson = gls(Value~dist*time.int, data = Rad.full, method = "ML")
summary(mod.pearson)

pearson.gau <- update(mod.pearson, correlation = corExp(1, form = ~Lat.dec+Long.dec+time.int), method = "ML")


#create bounding boxes#

#Large zone
coords = matrix(c(140.19, 38.174,
                  140.961, 38.174,
                  140.921, 38.003,
                  141.019, 37.737,
                  141.040, 37.358,
                  140.962, 36.958,
                  140.682, 36.612,
                  140.19, 36.612),
                ncol = 2, byrow = TRUE)
P1 <- Polygon(coords)
Ps1 = SpatialPolygons(list(Polygons(list(P1), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
plot(Ps1, axes = TRUE)




##### Get a p-value for an area ####

test.region <- expand.grid(gcoords1, gcoords2)

test.region <- test.region[sqrt(((grid[,2]-37.4213)*69)^2+((grid[,1]-141.0281)*cos(grid[,2])*69)^2)<5,]

############# KRIGING #############

#Small zone

coords.small <- read_csv("bbox_small.txt", col_names = FALSE)

P2 <- Polygon(cbind(coords.small$X2, coords.small$X1))
Ps2 = SpatialPolygons(list(Polygons(list(P2), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
plot(Ps2, axes = TRUE)

inner <- SpatialPoints(grid, proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

inner <- inner[Ps2,]

inner.df <- as.data.frame(inner)
inner.df <- inner.df %>% dplyr::rename(Long.dec = Var1, Lat.dec = Var2)

inner.df <- inner.df[sqrt(((inner.df[,2]-37.4213)*69)^2+((inner.df[,1]-141.0281)*cos(inner.df[,2])*69)^2)<18.6411 | inner.df[,2] > 37.5,]
#Get prediction locations
gcoords1 <- seq(min(Rad.ranks$Long.dec),max(Rad.ranks$Long.dec),0.01)
gcoords2 <- seq(min(Rad.ranks$Lat.dec),max(Rad.ranks$Lat.dec),0.01)

grid <- expand.grid(gcoords1, gcoords2)

grid <- grid[sqrt(((grid[,2]-37.4213)*69)^2+((grid[,1]-141.0281)*cos(grid[,2])*69)^2)<50,]

inner_grid <- locations.inside(grid, coords.small)


#Make geodata and get variogram estimate

#Raw Data
corr.full.geo <- as.geodata(Rad.ranks[,c("Long.dec", "Lat.dec", "rho")])

#Ranks
corr.geo <- as.geodata(Rad.ranks[,c("Long.dec", "Lat.dec", "resids")])
rad.v <- variog(corr.geo)
plot(rad.v)

fit1 <- variofit(rad.v, fix.nugget = F, weights = "npairs", cov.model="powered.exponential", fix.kappa = FALSE)
summary(fit1)

mlfit <- likfit(corr.geo, fix.kappa = FALSE, fix.nugget=FALSE, ini.cov.pars=fit1$cov.pars, nugget = fit1$nugget,
                kappa = fit1$kappa,
                lik.method = "REML", 
                cov.model="powered.exponential")

#krige correlation values on specified grid
pred<-krige.conv(corr.full.geo,locations=inner.df,
                 krige=krige.control(cov.model="powered.exponential",
                                     cov.pars=mlfit$cov.pars,
                                     nugget=mlfit$nugget,
                                     kappa = mlfit$kappa))

#plot
inner.df["est"] <- pred$predict

ggmap(map_image)+geom_tile(data = inner.df, aes(Long.dec, Lat.dec, fill = est))

   
###############################
####  BAYESIAN ANALYSIS  ######
###############################

library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

data.stan <- list(N = 9934, K =679, s = groupsizes$groupsize, y = Rad.full[,3:4])#, x = Rad.ranks[c("dist", "init")])

fit.rho.50 <- stan(file = 'Pearson_Bayes_hier_spearman.stan', data = data.stan, iter = 1000, chains = 2, refresh = 10, verbose = TRUE)

samples.full.pearson <- extract(fit.rho.full)
samples.full.hiervar <- extract(fit.rho.full) 
samples.rho <- load("rhosamples604_success2_small.RData")

save(samples.full.hiervar, file = "rhosamples_full_606141.RData")
save(fit.rho.hiervar, file = "rho_stanfit_full_606141.RData")
samples.full.hiervar <- load("rhosamples_full_606141.RData")


CI.bayes <- as.data.frame(t(apply(samples.full.pearson$rho, 2, FUN = quantile, probs = c(0.025,0.5, 0.95, 0.975))))
CI.bayes["ID"] <- as.numeric(row.names(CI.bayes))

for(i in 1:dim(unique(cbind(Rad.coarse$Lat.dec, Rad.coarse$Long.dec)))[1]){
  if(length(which(Rad.coarse$Lat.dec == locs.coarse[i,1] & Rad.coarse$Long.dec == locs.coarse[i,2])) == 6){
    Rad.full.coarse <- rbind(Rad.full.coarse, Rad.coarse[which(Rad.coarse$Lat.dec == locs.coarse[i,1] & Rad.coarse$Long.dec == locs.coarse[i,2]),])
  }
}

#highly aggregated test set.
Rad.test <- Rad.full %>% ungroup() %>% mutate(Lat.dec = round(Lat.dec,1), Long.dec=round(Long.dec,1)) %>% group_by(time.int, Lat.dec, Long.dec) %>%
                      summarise_each(funs(mean), Value) %>% 
                      group_by(Lat.dec, Long.dec) %>% 
                      mutate(groupsize=n(), rank.val = rank(Value), rank.time = rank(time.int)) %>%
                      ungroup() %>% arrange(groupsize, Lat.dec, Long.dec)

groupsizes <- Rad.test %>% group_by(Lat.dec, Long.dec) %>% summarise_each(funs(mean), groupsize) %>%
  ungroup() %>% dplyr::select(groupsize) %>% arrange(groupsize)


#misc
which(apply(samples.full.pearson$rho,2,median)>0)
which(Rad.ranks$pearson>0)
Rad.ranks[Rad.ranks$p_pearson > 0.05 & Rad.ranks$groupsize > 10,c("Lat.dec", "Long.dec")]$Lat.dec
agree <- intersect(intersect(which(Rad.ranks$CI_p_U>0), which(Rad.ranks$CI_s_U>0)), which(CI.bayes$`95%`>0))


##########################################
############## PLOTTING ##################
##########################################

diff.locs <- union(setdiff(which(Rad.ranks$CI_p_U>0 ),  which(CI.bayes$`95%`>0)), 
                   setdiff(which(CI.bayes$`95%`>0), which(Rad.ranks$CI_p_U>0)))

diff.locs <- union(union(which(((Rad.ranks$CI_p_U>0) != (CI.bayes$`95%`>0))),
                   which(((Rad.ranks$CI_p_U>0) != (Rad.ranks$CI_s_U>0)))),
                   which(((Rad.ranks$CI_s_U>0) != (CI.bayes$`95%`>0)))
                   )

df.pred <- cbind(inner_grid, pred$predict)

latlims <- c(max(Rad.ranks$Lat.dec), min(Rad.ranks$Lat.dec))
longlims <- c(max(Rad.ranks$Long.dec), min(Rad.ranks$Long.dec))

map_image <- get_map(location = c(lon = mean(longlims), lat = mean(latlims)), zoom = 10, source = "google")


Rad.clean <- Rad.clean %>% group_by(Lat.dec, Long.dec) %>% summarise_each(funs(mean), Value)
#initial plot of large zone
ggmap(map_image)+
  geom_tile(data = subset(Rad.clean, Value >= 629), aes(x = Long.dec, y = Lat.dec), fill = "red", color = "black")+
  geom_tile(data = subset(Rad.clean, Value < 629 & Value >= 50), 
            aes(x = Long.dec, y = Lat.dec, fill = Value), shape = 21, color = "black")+
  scale_fill_gradient("CPM", low = "blue", high = "orange")+
  geom_tile(data = subset(Rad.clean, Value < 50 & Value > 0), aes(x = Long.dec, y = Lat.dec), 
             fill = "green", shape = 21, color = "black")+
  geom_point(data = subset(Rad.clean, Value == 0), aes(x = Long.dec, y = Lat.dec), color = "black")+
  scale_size_continuous(name = "")+
  geom_point(aes(x = 141.03, y = 37.42), size = 6, color = "red")

#small zone
ggmap(map_image)+
  geom_tile(data = subset(Rad.clean.small, Value >= 629), aes(x = Long.dec, y = Lat.dec), fill = "red", color = "black")+
  geom_tile(data = subset(Rad.clean.small, Value < 629 & Value >= 50), 
            aes(x = Long.dec, y = Lat.dec, fill = Value), color = "black")+
  scale_fill_gradient("CPM", low = "blue", high = "orange")+
  geom_tile(data = subset(Rad.clean.small, Value < 50 & Value > 0), 
            aes(x = Long.dec, y = Lat.dec), 
                fill = "green", color = "black")+
  geom_point(data = subset(Rad.clean.small, Value == 0), aes(x = Long.dec, y = Lat.dec), color = "black")+
  scale_size_continuous(name = "")+
  geom_point(aes(x = 141.03, y = 37.42), size = 6, color = "black")

#Bad outlier
ggplot(data = Rad.full[1448:1455,], aes(time.int, Value))+geom_point(size = 4)+xlab("time")+ylab("CPM")
#Kriging maps
ggmap(map_image)+
  geom_tile(data = df.pred, aes(x = Var1, y = Var2, fill = pred$predict), shape = 21)+
  geom_point(aes(x = 141.03, y = 37.420), size = 6, color = "red")+
  scale_fill_gradient("correlation", low = "blue", high = "orange")

#Plot Correlation#

map_image <- get_map(location = c(lon = mean(longlims), lat = mean(latlims)), zoom = 10, source = "google")

ggmap(map_image)+
  geom_tile(data = Rad.ranks[Rad.ranks$pearson >= 0,], aes(x = Long.dec, y = Lat.dec),fill = "red", color = "black")+
  geom_tile(data = Rad.ranks[Rad.ranks$pearson < 0,], aes(x = Long.dec, y = Lat.dec, fill = pearson), color = "black")+
  geom_point(aes(x = 141.03, y = 37.42), size = 6, color = "red")+
  scale_fill_gradient2("Correlation", low = "green", mid = "blue", high = "orange", midpoint = -0.5)+
  geom_point(data = Rad.ranks[Rad.ranks$p_pearson > 0.05 & Rad.ranks$groupsize > 8,c("Lat.dec", "Long.dec")], 
             aes(x=Long.dec, y=Lat.dec), fill = "yellow", shape = 25, size = 4)+
  theme(legend.position = "none")

ggmap(map_image)+
  geom_tile(data = Rad.ranks[Rad.ranks$rho >= 0,], aes(x = Long.dec, y = Lat.dec),fill = "red", color = "black")+
  geom_tile(data = Rad.ranks[Rad.ranks$rho < 0,], aes(x = Long.dec, y = Lat.dec, fill = rho), color = "black")+
  geom_point(aes(x = 141.03, y = 37.42), size = 6, color = "red")+
  geom_point(data = Rad.ranks[Rad.ranks$p_spearman > 0.05 & Rad.ranks$groupsize > 8,c("Lat.dec", "Long.dec")], 
             aes(x=Long.dec, y=Lat.dec), fill = "yellow", shape = 25, size = 4)+
  scale_fill_gradient2("Spearman Correlation", low = "green", mid = "blue", high = "orange", midpoint = -0.5)+
  theme(legend.position = "none")

#points of disagreement
ggmap(map_image)+
  geom_tile(data = Rad.ranks[Rad.ranks$pearson >= 0,], aes(x = Long.dec, y = Lat.dec),fill = "red", color = "black")+
  geom_tile(data = Rad.ranks[Rad.ranks$pearson < 0,], aes(x = Long.dec, y = Lat.dec, fill = pearson), color = "black", shape = 21)+
  geom_point(aes(x = 141.03, y = 37.42), size = 6, color = "red")+
  geom_point(data = Rad.ranks[subset(Rad.ranks[diff.locs,], groupsize > 10)$ID,c("Lat.dec", "Long.dec")], 
             aes(x=Long.dec, y=Lat.dec), fill = "red", shape = 25, size = 4)+
 # geom_point(data = subset(Rad.ranks[agree,], groupsize >10), aes(x=Long.dec, y=Lat.dec), color = "green", size = 4)+
  scale_fill_gradient2("Spearman Correlation", low = "green", mid = "blue", high = "orange", midpoint = -0.5)+
  theme(legend.position = "none")

#points of agreement
ggmap(map_image)+
  geom_tile(data = Rad.ranks[Rad.ranks$pearson >= 0,], aes(x = Long.dec, y = Lat.dec),fill = "red", color = "black")+
  geom_tile(data = Rad.ranks[Rad.ranks$pearson < 0,], aes(x = Long.dec, y = Lat.dec, fill = pearson), color = "black", shape = 21)+
  geom_point(aes(x = 141.03, y = 37.42), size = 6, color = "red")+
  #geom_point(data = Rad.ranks[subset(Rad.ranks[diff.locs,], groupsize > 10)$ID,c("Lat.dec", "Long.dec")], 
  #           aes(x=Long.dec, y=Lat.dec), color = "red", shape = 17, size = 4)+
  geom_point(data = subset(Rad.ranks[agree,], groupsize >10), aes(x=Long.dec, y=Lat.dec), fill = "yellow", shape = 25, size = 4)+
  scale_fill_gradient2("Spearman Correlation", low = "green", mid = "blue", high = "orange", midpoint = -0.5)+
  theme(legend.position = "none")

#insignificant points for classical analysis
ggmap(map_image)+
  geom_tile(data = Rad.ranks, aes(x = Long.dec, y = Lat.dec, fill = rho), color = "black", shape = 21)+
  geom_point(aes(x = 141.03, y = 37.42), size = 6, color = "red")+
  geom_point(data = Rad.ranks[Rad.ranks$CI_s_U > 0 & Rad.ranks$groupsize > 8,c("Lat.dec", "Long.dec")], 
             aes(x=Long.dec, y=Lat.dec), color = "red", shape = 17, size = 4)+
  scale_fill_gradient2("correlation", low = "green", high = "yellow")

ggmap(map_image)+
  geom_tile(data = Rad.ranks, aes(x = Long.dec, y = Lat.dec, fill = rho), color = "black", shape = 21)+
  geom_point(aes(x = 141.03, y = 37.42), size = 6, color = "red")+
  geom_point(data = Rad.ranks[Rad.ranks$p_spearman > 0.05 & Rad.ranks$groupsize > 8,c("Lat.dec", "Long.dec")], 
             aes(x=Long.dec, y=Lat.dec), color = "red", shape = 17, size = 4)+
  scale_fill_gradient2("correlation", low = "green", high = "yellow")
#insignificant points for bayesian analysis
ggmap(map_image)+
  geom_tile(data = Rad.ranks, aes(x = Long.dec, y = Lat.dec, fill = pearson), shape = 21)+
  geom_point(aes(x = 141.03, y = 37.42), size = 6, color = "red")+
  geom_point(data = subset(Rad.ranks[which(CI.bayes$`95%`>0),c("Lat.dec", "Long.dec", "groupsize")], groupsize > 10), 
             aes(x=Long.dec, y=Lat.dec), color = "red", shape = 17, size = 4)+
  scale_fill_gradient("correlation", low = "blue", high = "yellow")

##
map_image <- get_map(location = c(lon = mean(longlims), lat = mean(latlims)), zoom = 9, source = "google")

ggmap(map_image)+
  geom_tile(data = Rad.full, aes(x = Long.dec, y = Lat.dec, fill = Value), shape = 21)+
  geom_point(aes(x = 141.03, y = 37.42), size = 6, color = "red")+
  scale_fill_gradient("p-values", low = "blue", high = "orange")

#CI's
ggplot(data = Rad.ranks, aes(x= ID, y = pearson))+geom_point(col="blue", size = 2)+
  geom_errorbar(aes(ymin = CI_p_L,ymax = CI_p_U),width = 0.2,colour = 'red')+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 224, col="green", size = 2)+ 
  ylab("Pearson Correlation")+
  theme(panel.background = element_rect(fill = 'grey', colour = 'grey'), panel.grid = element_line(color = "grey"))

ggplot(data = CI.bayes, aes(x= ID, y = `50%`))+geom_point(col="blue", size = 2)+
  geom_errorbar(aes(ymin = -1,ymax = `95%`),width = 0.2,colour = 'red')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = median(samples.full.hiervar$gamma))+
  geom_vline(xintercept = 224, col="green", size = 2)+
  ylab("Pearson Correlation")+
  theme(panel.background = element_rect(fill = 'grey', colour = 'grey'), panel.grid = element_line(color = "grey"))


boxplotdf <- melt(samples.full.pearson$rho[,seq(1,679,5)])
boxplotdf <- boxplotdf %>% group_by(Var.2) %>% mutate(upper = quantile(value, 0.95)) %>% mutate(sig = upper > 0)

ggplot(data = boxplotdf, aes(x = Var.2, y = value, group = Var.2, fill = sig))+geom_boxplot()+
  geom_hline(yintercept = 0)+
  #geom_hline(yintercept = median(samples.full.pearson$gamma), col = "red")+
  theme(legend.position = "none")


################
####DUMP########
################

mods <- Rad.full %>% slice_rows(c("Lat.dec", "Long.dec")) %>% do(model = lm(Value~time.int, data = .))
qplot(fitted(mods$model[[6]]), residuals(mods$model[[6]]))

####Nonlinear Mixed Effects####

#rad.gnlsmod <- gnls(Value ~ z+1/dist^a+b*time.int, data = Rad.full, start = list(z=1, a=2, b=-1))
#summary(rad.gnlsmod)

#gnlsmod.cor <- gnls(Value ~ z+1/dist^a+b*time.int, correlation = corExp(1,~Lat.dec+Long.dec|time.int), data = Rad.full, start = list(z=1, a=2, b=-1))
#

#rad.gls <- gls(log(Value)~dist*time.int, data = Rad.full)
#summary(rad.gls)

#vario <- Variogram(rad.gls, form = ~Lat.dec + Long.dec, resType = "pearson")

#####Method 2:  Calculate correlation at each prediction location with kriged surfaces#####

#Create geodata

#resids.rad <- as.geodata(Rad.clean[Rad.clean$year == 5,c("Long.dec", "Lat.dec", "resids")])

#rad.v <- variog(resids.rad, max.dist = 0.5)
#plot(rad.v)

#fit1 <- variofit(rad.v, fix.nugget = F, weights = "npairs", cov.model="powered.exponential", fix.kappa = FALSE)
#summary(fit1)

#mlfit <- likfit(resids.rad, fix.kappa = FALSE, fix.nugget=FALSE, ini.cov.pars=fit1$cov.pars, nugget = fit1$nugget,
#                kappa = fit1$kappa,
#                lik.method = "REML", 
#                cov.model="matern")

#pred<-krige.conv(resids.rad,locations=inner_grid,
#                 krige=krige.control(cov.model="powered.exponential",
#                                     cov.pars=fit1$cov.pars,
#                                     nugget=fit1$nugget,
#                                     kappa = fit1$kappa))

#inner_grid["estimates"] <- pred$predict

#ggplot(aes(Var1, Var2, fill = estimates), data = inner_grid) +geom_tile()+
#  scale_fill_distiller(palette = "Oranges", name = "Log Radiation", direction = 1)

#write.csv(inner_grid, "D:/Data/Kriging_y5.csv")

#Stack data from all kriging surfaces

#df.rho <- inner_grid %>% mutate(y_0 = read_csv("Data/Kriging_y0.csv")$estimates,
#                                y_1 = read_csv("Data/Kriging_y1.csv")$estimates,
#                                y_2 = read_csv("Data/Kriging_y2.csv")$estimates,
#                                y_3 = read_csv("Data/Kriging_y3.csv")$estimates,
#                                y_4 = read_csv("Data/Kriging_y4.csv")$estimates,
#                                y_5 = read_csv("Data/Kriging_y5.csv")$estimates)

