
###Set to wherever you're working, or use "here" library and skip
setwd("I:\\workspace\\fire-size-distribution")

library("RSQLite")
library(RODBC)
library(sf)
library(terra)
library(spatial)
library(sp)
library(spatstat)
library(VGAM)
library(fitdistrplus)
library(stringr)
library(readxl)

###Helper functions


# computing the two unique vectors with ties present: the function is tauuniq
tauuniq <- function(x,y) {
  n <- length(x)
  e <- 1:n
  xrr <- n+1 -rank(x)
  xtp <- x[order(y,x)]
  xtn <- x[order(y,xrr)]
  rkyp <- order(xtp,e)
  rkyn <- order(xtn,n:1)
  out <- cbind(rkyp,rkyn)
  out }
# calculation of Kendall's tau on unique max-min vectors
# the function is rtau
rtau <- function(x,y){
  ot <- tauuniq(x,y)
  rkyp <- ot[,1]
  rkyn <- ot[,2]
  dyp <- 0
  dyn <- 0
  n <- length(x)
  n2 <- ((n*(n-1))/2)
  n1 <- n-1
  for(i in 1:n1) {j <- i+1
  tempp <- rkyp[i]-rkyp[j:n]
  tempn <- rkyn[i]-rkyn[j:n]
  dyp <- dyp + sum(tempp<0)
  dyn <- dyn + sum(tempn<0)
  }
  out <- (dyp + dyn)/n2 -1
  out }
# output is slope and intercept, function name tauslp in positions 1 and 2
tauslp <- function(x,y) {
  rat <- c(outer(y,y,"-")/outer(x,x,"-"))
  ratv <- rat[!is.na(rat)]
  slp <- median(ratv)
  res <- y - slp * x
  aint <- median(res)
  res <- res - aint
  ck <- rtau(x,res)
  ck1 <- sum(res)
  ck2 <- median(res)
  out <- c(slp,aint,ck,ck1,ck2)
  out }

###Estimating pareto parameters using maximum likelihood

pareto.MLE <- function(X)
{
  n <- length(X)
  m <- min(X)
  a <- n/sum(log(X)-log(m))
  return( c(m,a) ) 
}

###Estimating truncated pareto parameters using maximum likelihood


bounded_powerlaw_mle <- function(data, xmin, xmax) {
  x <- data[data >= xmin & data <= xmax]
  n <- length(x)
  
  log_likelihood <- function(alpha) {
    if (alpha <= 0) return(-Inf)
    term1 <- n * log(alpha)
    term2 <- n * alpha * log(xmin)
    term3 <- -(alpha + 1) * sum(log(x))
    term4 <- -n * log(1 - (xmin / xmax)^alpha)
    return(term1 + term2 + term3 + term4)
  }
  
  opt <- optimize(function(a) -log_likelihood(a), interval = c(0.01, 10))
  return(opt$minimum)
}

# 
# 
# ## CA cleaned data
# 
# dta <- odbcConnectAccess2007("FOD_EXCERPT_CA_NM_NV_20250404.accdb")   #specifies the file path
# 
# sqlTables(dta, tableType = "TABLE")$TABLE_NAME
# 
# Fires <- sqlFetch(dta, "Fires")
# 
# 
# dta[[1]]
# 
# 
# #df1 <- sqlFetch(dta)  
# #odbcCloseAll()
# dta$TABLE_NAME
# con <- odbcConnect("FOD_EXCERPT_CA_NM_NV_20250404")
# 
# con <- odbcConnect("FOD_EXCERPT_CA_NM_NV_20250404")
# 
# db<-file.path("I:\\workspace\\fire-size-distribution\\FOD_EXCERPT_CA_NM_NV_20250404.accdb") #connect database.
# channel<-odbcConnectAccess2007(db)
# myConn <-odbcDriverConnect("Driver={Oracle in OraClient11g_home1};Dbq=NZSQL;Uid=cejacobson;Pwd=password;")
# 
# tbls <- sqlTables(dta)
# tbls$TABLE_NAME
# 
# qry <- "SELECT * FROM `Fires`"
# 
# fires <- sqlQuery(con, qry)
# con <- odbcConnect("FOD_EXCERPT_CA_NM_NV_20250404")
# 
# con <- dbConnect(drv=RSQLite::SQLite(), dbname="FPA_FOD_20221014.sqlite")
# 
# 
# ## connect to db
# con <- dbConnect(drv=RSQLite::SQLite(), dbname="FOD_EXCERPT_CA_NM_NV_20250404.accdb")


###Reading in FOD to get historical distribution

filename <- "FPA_FOD_20221014.sqlite"
sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver,
                dbname = filename)

## Some operations
dbListTables(db)
mytable <- dbReadTable(db,"Fires")

years_all <- mytable$FIRE_YEAR
acres_all <- mytable$FIRE_SIZE

yearvec <- sort(unique(mytable$FIRE_YEAR))
nyears <- length(yearvec)

slvec <- rep(NA, nyears)
parvec  <- rep(NA, nyears)

###Getting pyrome numbers for FOD

pyrome_shp  <- st_read("Pyromes_CONUS_20200206.shp")
pyrome_crs <- crs(pyrome_shp)
v <- vect("Pyromes_CONUS_20200206.shp")

v_crs <- crs(v)
#v_p4s <- proj4string(v_crs)
coords.latlong <- cbind(mytable$LONGITUDE,mytable$LATITUDE)
coords.sp <- SpatialPoints(coords.latlong)
P4S.latlon <- CRS("+proj=longlat +datum=WGS84")
proj4string(coords.sp) <- P4S.latlon
coords_albers <- spTransform(coords.sp,pyrome_crs )
p4s.pyrome <- st_crs(pyrome_shp)
nboot <- 1000

plot(coords_albers, pch = ".")
coord.vals <- coords_albers@coords

pyrome_vec <- extract(v, coord.vals)


### Get burnable areas

burnable_acres <- read.csv("fsim-pyrome-burnablearea.csv")



###Instantiate outputs

pyrome_vec_num <- pyrome_vec$PYROME
pyrome_vec_name<- pyrome_vec$NAME
pyromes_name_unique <- sort(unique(pyrome_vec_name))

pyromes_unique <- sort(unique(pyrome_vec_num))
npyromes <- length(pyromes_unique)

maxvec <- burnable_acres$burnable_ac


parmat <- matrix(NA, npyromes, 11)
years_unique <- sort(unique(mytable$FIRE_YEAR))
nfiremat <- matrix(NA, npyromes, length(years_unique))
length(years_unique)
parlist <- vector("list", npyromes)



nfiremat <- matrix(NA, npyromes, length(years_unique))
length(years_unique)
parlist <- vector("list", npyromes)
shapelist <- vector("list", npyromes)
nboot <- 1000

ci_list_hist <- vector("list", npyromes)
ci_list_sim <- vector("list", npyromes)

for(curpyrome in 1:npyromes)
{
  
  ###Get FSIM output
  setwd("I:\\workspace\\fsim-climate\\Calibration workbooks\\Calibration workbooks")
  cwd <- getwd()
  pynum <- str_pad(curpyrome, 3, pad = "0")
  
  fname <- paste("PY", pynum, "_FSimCalibration_v2022.12.xlsx", sep = "")
  if(file.exists(fname))
    {
    
    submat <- mytable[pyrome_vec_num == curpyrome,]
    navec <- is.na(submat$FIRE_SIZE)
    numvec <- as.logical(1-navec)
    submat <- submat[numvec,]
    yearvec <- submat$FIRE_YEAR
    #subsetting to only last half of dataset
    submat <- submat[yearvec > 2004,]
    yearvec <- submat$FIRE_YEAR
    
    years_unique <- sort(unique(yearvec))
    nyears <- length(years_unique)
    tempmax <- maxvec[curpyrome]
    
    parmat_year <- matrix(NA, length(years_unique), 11)
    
    burn_acres_row <- burnable_acres[burnable_acres$PYROME == curpyrome,]


    tempmax <- burn_acres_row$burnable_ac
    tempmin <- 0.1
    tempmin_fsim <- 18
    inc <- 1
    
    tempshape_vec <- rep(NA, length(years_unique))
    
    for(curyear in years_unique)
    {
      print(c(curpyrome, curyear))
      submat_year <- submat[yearvec == curyear,]
      submat_acres <- na.exclude(submat_year$FIRE_SIZE)
      ##Fixing all fires to 0.1 acres
      if(sum(submat_acres > 0))
      {
        submat_acres[submat_acres < tempmin] <- tempmin
        
        submat_gte18 <- submat_acres
        submat_gte18 <- submat_gte18[submat_gte18 >= 18]
        temp_shape <- bounded_powerlaw_mle(submat_gte18, tempmin_fsim, tempmax)
        tempshape_vec[inc] <- temp_shape
        # vgam_shape <- fitdist(submat_acres, "truncpareto",
        #                       start = list(lower=tempmin,
        #                                    upper=tempmax,
        #                                    shape=1.0))
        nfiremat[curpyrome, inc] <-dim(submat_year)[1]
      }
      inc <- inc + 1
    }
    
    shapelist[[curpyrome]] <- tempshape_vec
    
    ###Get shape parameter from FSIM output
    acres_hist <- submat$FIRE_SIZE
    acres_gte18 <- acres_hist[acres_hist >= 18]
    namevec <- "Hist"
    if(length(acres_gte18) > 30)
    {
      breaks <- 10^seq(log10(18), log10(tempmax), length.out = 30)
      h <- hist(acres_gte18, breaks = breaks, plot = FALSE)
      x <- h$mids
      y <- h$density
      valid <- y > 0
      x <- x[valid]
      y <- y[valid]
      ycum <-cumsum(y)
      lx <- log(x)
      ly <- log(y)
      lyc <- 1-cumsum(ycum)
      
      mx <- mean(lx)
      my <- mean(ly)
      
      #fname <- paste(cwd, fname, sep = "/")
      fsim_acre_df <- read_excel(fname, sheet = 2)
      fsim_acre_df <- data.frame(fsim_acre_df)
      first_row <- fsim_acre_df[1,]
      valid_vec <- as.logical(1-is.na(first_row))
      indvec <- 1:length(valid_vec)
      last_run <- max(indvec[valid_vec])
      
      alpha_vec_mean <- numeric(last_run)
      alpha_vec_se <- numeric(last_run)
      
      namevec <- c(namevec, names(first_row)[valid_vec])
      
      for(currun in 1:last_run)
      {
        fsim_acres <- fsim_acre_df[,currun]
        fsim_acres <- na.exclude(fsim_acres)
        fsim_shape <- bounded_powerlaw_mle(fsim_acres, tempmin_fsim, tempmax)
        alpha_boot_hist <- rep(NA, nboot)
        alpha_boot_sim <- rep(NA, nboot)
        
        for (b in 1:nboot) {
          histx <- sample(acres_gte18, replace = TRUE)
          simx <- sample(fsim_acres, length(acres_gte18), replace = TRUE)
          
          # fit MLE to simx
          alpha_boot_hist[b] <- bounded_powerlaw_mle(histx, xmin = tempmin_fsim, xmax = tempmax)
          alpha_boot_sim[b] <- bounded_powerlaw_mle(simx, xmin = tempmin_fsim, xmax = tempmax)
          
        }
        
        
        se_boot_hist <- sd(alpha_boot_hist, na.rm = TRUE)
        se_boot_sim <- sd(alpha_boot_sim, na.rm = TRUE)
        
        mean_boot_hist <- mean(alpha_boot_hist)
        mean_boot_sim <- mean(alpha_boot_sim)
        
        alpha_vec_mean[currun] <- mean_boot_sim
        alpha_vec_se[currun] <- se_boot_sim
        
        
      }
      
      alpha_mean_out <- c(mean_boot_hist, alpha_vec_mean)
      alpha_se_out <- c(se_boot_hist, alpha_vec_se)
      
      outmat_alpha <- cbind(alpha_mean_out, alpha_se_out)
      rownames(outmat_alpha) <- namevec
      fout <- paste( "outmat_alpha_", curpyrome, ".txt", sep = "")
      write.table(outmat_alpha, fout)
      
      
      rownames(outmat_alpha) <- c("Hist", as.character(1:last_run))
      z_test <- (mean_boot_hist - mean_boot_sim) / (se_boot_hist**2 + se_boot_sim**2)
      
      alpha_upper_hist <- mean_boot_hist + 1.96*se_boot_hist
      alpha_lower_hist <- mean_boot_hist - 1.96*se_boot_hist 
      
      alpha_upper_sim <- mean_boot_sim + 1.96*se_boot_sim
      alpha_lower_sim <- mean_boot_sim - 1.96*se_boot_sim 
      
      ci_list_hist[[curpyrome]] <- c(alpha_lower_hist, alpha_upper_hist)
      ci_list_sim[[curpyrome]] <- c(alpha_lower_sim, alpha_upper_sim)
      
      slope_fit_mean  <- -(mean_boot_hist+1)
      slope_fit_upper <- -(mean_boot_hist + 1.96*se_boot_hist + 1)
      slope_fit_lower <- -(mean_boot_hist - 1.96*se_boot_hist + 1)
      
      h_sim <- hist(fsim_acres, breaks = breaks, plot = FALSE)
      ###K-L divergence of densities
      dens_hist <- h$density
      dens_sim <- h_sim$density
      #x_dens <- rbind(dens_hist,dens_sim)
      #KL(x_dens, unit = "log10")
      #ks.test(acres_gte18, fsim_acres)
      
      
      
      fsim_slope <- -(fsim_shape + 1)
      
      int_mean_hist <- log(mean_boot_hist) + mean_boot_hist * log(tempmin_fsim) - log(1 - (tempmin_fsim / tempmax)^mean_boot_hist)
      int_upper_hist <- log(alpha_upper_hist) + alpha_upper_hist * log(tempmin_fsim) - log(1 - (tempmin_fsim / tempmax)^alpha_upper_hist)
      int_lower_hist <- log(alpha_lower_hist) + alpha_lower_hist * log(tempmin_fsim) - log(1 - (tempmin_fsim / tempmax)^alpha_lower_hist)
      int_fsim <- log(fsim_shape) + fsim_shape * log(tempmin_fsim) - log(1 - (tempmin_fsim / tempmax)^fsim_shape)
      
      
      fout <- paste("fsd-plot-", curpyrome, "-late.png", sep = "")
      png(fout)
      plot(lx,ly)
      
      
      abline(a = int_mean_hist, b = slope_fit_mean, col = "blue")
      abline(a = int_upper_hist, b = slope_fit_upper, col = "blue", lty = 2)
      abline(a = int_lower_hist, b = slope_fit_lower, col = "blue", lty = 2)
      abline(a = int_fsim, b = fsim_slope, col ="red", lty = 2)
      dev.off()
      ##Kendalls tau fit
      #plot(lx, ly)
      
      fout <- paste("shape_param_", curpyrome, "_late.png", sep = "")
      png(fout)
      plot(years_unique, shapelist[[curpyrome]], xlab  = "Year", ylab = "est-shape", ylim = c(0, 2))
      
      abline(h = fsim_shape, col = "red")
      abline(h = mean_boot_hist, col = "blue")
      abline(h = mean_boot_hist + 1.96*se_boot_hist, col = "blue", lty = 2)
      abline(h = mean_boot_hist - 1.96*se_boot_hist, col = "blue", lty = 2)
      
      dev.off()
    }
    
  }
  
}

plot(nfiremat[,1], xlab=  "Year")
plot(years_unique, shapelist[[128]], xlab  = "Year", ylab = "est-shape")


ci_mat_hist <- do.call("rbind", ci_list_hist)
ci_mat_sim<- do.call("rbind", ci_list_sim)

n_valid <- dim(ci_mat_hist)[1]
x <- 1:(2*n_valid)
x <- c(x, max(x) + 1)
x_hist <- seq(1, 2*n_valid, by = 2)
x_sim <- seq(2, 2*n_valid+1, by = 2)
col_seq <- rep(c("red", "blue"), n_valid)
#plot(ci_mat_hist[,1], ci_mat_sim[,1], col = "blue")
ci_mid_hist <- (ci_mat_hist[,2] - ci_mat_hist[,1])/2 + ci_mat_hist[,1]
ci_mid_sim <- (ci_mat_sim[,2] - ci_mat_sim[,1])/2 + ci_mat_sim[,1]
interleaved_vec <- as.vector(rbind(ci_mid_hist, ci_mid_sim))
interleaved_L <- as.vector(rbind(ci_mat_hist[,1], ci_mat_sim[,1]))
interleaved_U <- as.vector(rbind(ci_mat_hist[,2], ci_mat_sim[,2]))
inc <- 10
inc_upper <- inc + 10
x <- x/2

points(ci_mat_hist[,2], ci_mat_sim[,2], col = "red")

outmat <- cbind(ci_mat_hist, ci_mat_sim)
names(outmat) <- c("Hist-L", "Hist-U", "Sim-L", "Sim_U")
write.table(outmat, "CI-Hist-Sim.csv")

