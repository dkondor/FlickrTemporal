# flickr_ts2.r -- run time series analysis steps for all Flickr data

library(ggplot2)
library(reshape2)

setwd("~/senseable/Data/flickr_timeseries/")
cities = c(0,1,2,3,4,5,7,8,10,46)

# steps:
# 0. read file, aggregation (daily / 4-hour), check for missing dates, fill these with
#   either zero or NAs
# 1. Fourier-spectrum (daily / 4-hour data)
# 2. decomposition (additive / multiplicative, 7-day / 28-day periods)
# 3. temporal autocorrelation in the residuals
# 4. distribution of noise, test for normality
# 5. distribution of normalized residuals, average values for these

dir.create("run1")
setwd("run1")


# perform a multiplicative model (i.e. x = seasonal*trend + error) decomposition of a 
#   given time series (ts1 should be given as a vector)
# !!notes: this might not work if length(ts1) is not an exact multiple of the period,
#   (I have not tested that configuration)
# also, this might be a bit slow, as the code is not optimal (uses loops to go thorugh
#   the data as my R skills are still a bit limited in some aspects)
mdecompose <- function(ts1, period) {
  len = length(ts1)
  ws = period
  ts5 = data.frame(hour4 = 1:len, cnt = ts1)
  n1 = aggregate(ts5["cnt"], by = data.frame(x = floor((ts5$hour4-1)/ws)), FUN = sum)
  n12 = ts5
  for(i in 1:len) {
    n12$cnt[i] = n12$cnt[i] / n1$cnt[floor((i-1)/ws)+1]
  }
  a1 = aggregate(n12$cnt[1:len], by = data.frame(d = 0:(len-1)%%ws), FUN = mean)
  pr1 = data.frame(hour4 = 1:len, cnt = 1:len)
  for(i in 1:len) {
    pr1$cnt[i] = n1$cnt[floor((i-1)/ws+1)] * a1$x[ (i-1)%%ws + 1 ]
  }
  
  t1 = data.frame(hour4 = 1:len, cnt = 1:len)
  for(i in 1:len) {
    t1$cnt[i] = n1$cnt[floor((i-1)/ws)+1] / ws
  }
  s1 = data.frame(hour4 = 1:len, cnt = 1:len)
  for(i in 1:len) {
    s1$cnt[i] = a1$x[ (i-1)%%ws + 1 ]
  }
  
  tsall1 = ts5[1:len,]
  tsall1$type = "observed"
  t1$type = "trend"
  tsall1 = rbind(tsall1,t1)
  s1$type = "seasonal"
  s1$cnt = s1$cnt * mean(ts5$cnt[1:len]) * ws
  tsall1 = rbind(tsall1,s1)
  
  e1 = merge(ts5[1:len,],pr1,by="hour4")
  e1$err = e1$cnt.x - e1$cnt.y
  tsall1 = rbind(tsall1, data.frame(hour4 = e1$hour4, cnt = e1$err, type = "random"))
  pr1$type = "predicted"
  tsall1 = rbind(tsall1,pr1)
  return(tsall1)
}


for(c in cities) {
  for(b in c("R","T")) {
    # 0. read file, aggregate, check for missing dates / times
    fn = paste0("../c_",c,"_",b,".txt")
    t1 = read.table(fn)
    # !! note: I replaced all the commas separating the coordinates with spaces in the files !!
    if(b == "R")
      names(t1) = c("null","photoid","day","time","lat","lon")
    else
      names(t1) = c("null","photoid","day","time","lat","lon","country")
    
    dir1 = paste0("c_",c,"_",b,"_1day")
    dir2 = paste0("c_",c,"_",b,"_4hour")
    dir.create(dir1)
    dir.create(dir2)
    # daily aggregation
    ts1 = aggregate(t1["null"],by=t1["day"],FUN=length)
    ts1 = data.frame(day = as.POSIXct(ts1$day,tz="UTC"), cnt = ts1$null)
    # 4-hour aggregation
    t1$timestamp = paste0(t1$day,"T",t1$time)
    t1$timestamp2 = as.numeric(strptime(t1$timestamp,"%Y-%m-%dT%H:%M:%S",tz="UTC"))
    t1$hour4 = 14400*floor(t1$timestamp2 / 14400)
    ts2 = aggregate(t1["null"],by=t1["hour4"],FUN=length)
    names(ts2)[2] = "cnt"
    
    # check for missing dates / times
    # all ts start 2007-01-01 and finish 2009-12-31 (incl.)
    ts1start = as.numeric(strptime("2007-01-01","%Y-%m-%d",tz="UTC"))
    ts1end = as.numeric(strptime("2009-12-31","%Y-%m-%d",tz="UTC"))
    ts2end = as.numeric(strptime("2009-12-31T20:00","%Y-%m-%dT%H:%M",tz="UTC"))
    n1 = (ts1end - ts1start) / 86400
    n2 = (ts2end - ts1start) / 14400
    
    ts1t = data.frame(day = as.POSIXct(ts1start + 0:n1*86400, origin = '1970-01-01', tz = 'UTC'))
    ts2t = data.frame(hour4 = ts1start + 0:n2*14400)
    
    ts1 = merge(ts1t,ts1,by="day",all.x=TRUE)
    ts2 = merge(ts2t,ts2,by="hour4",all.x=TRUE)
    
    ts1[is.na(ts1$cnt),"cnt"] = 0
    ts1 = ts1[with(ts1,order(day)),]
    ts2[is.na(ts2$cnt),"cnt"] = 0
    ts2 = ts2[with(ts2,order(hour4)),]
    
    
    #####################################################
    # 1. calculate Fourier-transform, create figs of it
    fft1 = fft(ts1$cnt[1:1092])
    fft2 = fft(ts2$cnt[1:6552])
    # 6552/2 # 3276
    # 1092/2 # 546
    tmp12 = Re(fft1[1:546]*Conj(fft1[1:546]))
    tmp22 = Re(fft2[1:3276]*Conj(fft2[1:3276]))
    tmp13 = data.frame(f = 0:545, w = tmp12)
    tmp23 = data.frame(f = 0:3275, w = tmp22)
    
    p1 = ggplot(tmp13) + geom_line(aes(x=f,y=w),size=0.35) + scale_y_continuous(trans="log") +
      scale_x_continuous(trans="log", breaks = c(3,12,36.4,156,312),
      labels = c("1 year","1/4 year","1 month","1 week","1/2 week"))
    p1 = p1 + theme_bw(base_size=7) + theme(axis.text.y=element_blank()) + xlab("period") + 
      ylab("weight (arbitrary units, log scale)")
    
    p2 = ggplot(tmp23) + geom_line(aes(x=f,y=w),size=0.35) + scale_y_continuous(trans="log") +
      scale_x_continuous(trans="log", breaks = c(3,12,36.4,156,312,1092),
      labels = c("1 year","1/4 year","1 month","1 week","1/2 week","1 day"))
    p2 = p2 + theme_bw(base_size=7) + theme(axis.text.y=element_blank()) + xlab("period") + 
      ylab("weight (arbitrary units, log scale)")
    
    ggsave(paste0(dir1,"/fourier.pdf"),p1,width=3.2,height=2)
    ggsave(paste0(dir2,"/fourier.pdf"),p2,width=3.2,height=2)
    ggsave(paste0(dir1,"/fourier.png"),p1,width=3.2,height=2,dpi=220)
    ggsave(paste0(dir2,"/fourier.png"),p2,width=3.2,height=2,dpi=220)
    
    # write out top 10 components to files
    tmp14 = head(tmp13[order(tmp13$w,decreasing = TRUE),],n = 10)
    tmp14$f = 1092 / tmp14$f
    names(tmp14) = c("period (days)","weight")
    write.table(tmp14,paste0(dir1,"/fouriercomponents.txt"),col.names=TRUE,row.names=FALSE)
    tmp24 = head(tmp23[order(tmp23$w,decreasing = TRUE),],n = 10)
    tmp24$f = 1092 / tmp24$f
    names(tmp24) = c("period (days)","weight")
    write.table(tmp24,paste0(dir2,"/fouriercomponents.txt"),col.names=TRUE,row.names=FALSE)
    
    #########################################################
    ## 2.1 additive decomposition, 7-day / 28-day period
    ts21 = ts(ts2$cnt[1:6552], frequency = 42, start = c(1,1))
    ts22 = ts(ts2$cnt[1:6552], frequency = 168, start = c(1,1)) # 28-day period
    ts11 = ts(ts1$cnt[1:1092], frequency = 7, start = c(1,1))
    ts12 = ts(ts1$cnt[1:1092], frequency = 28, start = c(1,1)) # 28-day period
    
    c11 = decompose(ts11)
    c12 = decompose(ts12)
    c21 = decompose(ts21)
    c22 = decompose(ts22)
    
    # plots for the daily aggregated data
    c11 = data.frame(day = ts1$day[1:1092], observed = c11$x, seasonal = c11$seasonal,
      trend = c11$trend, random = c11$random, predicted = c11$seasonal + c11$trend)
    c111 = melt(c11,id.vars = "day")
    p1 = ggplot(c111) + geom_line(aes(x=day,y=value)) + facet_grid(variable~.) + theme_bw()
    ggsave(paste0(dir1,"/additive_model_7day.pdf"),width=6,height=6)
    ggsave(paste0(dir1,"/additive_model_7day.png"),width=6,height=6,dpi = 150)
    c12 = data.frame(day = ts1$day[1:1092], observed = c12$x, seasonal = c12$seasonal,
      trend = c12$trend, random = c12$random, predicted = c12$seasonal + c12$trend)
    c121 = melt(c12,id.vars = "day")
    p1 = ggplot(c121) + geom_line(aes(x=day,y=value)) + facet_grid(variable~.) + theme_bw()
    ggsave(paste0(dir1,"/additive_model_28day.pdf"),width=6,height=6)
    ggsave(paste0(dir1,"/additive_model_28day.png"),width=6,height=6,dpi = 150)
    
    # plot the seasonal TS
    p1 = ggplot(data.frame(day = 1:7, activity = c11$seasonal[1:7])) + geom_line(aes(
      x=day,y=activity)) + theme_bw(base_size = 8) +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    # note: ts starts at 2007-01-01, which was a Monday
    p1 = p1 + scale_x_continuous(breaks=1:7,labels=c("Mon","Tue","Wed","Thu","Fri","Sat","Sun"))
    ggsave(paste0(dir1,"/additive_seasonal7.pdf"),width=2.5,height=2)
    ggsave(paste0(dir1,"/additive_seasonal7.png"),width=2.5,height=2,dpi=220)
    
    p1 = ggplot(data.frame(day = 1:28, activity = c12$seasonal[1:28])) + geom_line(aes(
      x=day,y=activity)) + theme_bw(base_size = 8) +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    ggsave(paste0(dir1,"/additive_seasonal28.pdf"),width=2.5,height=2)
    ggsave(paste0(dir1,"/additive_seasonal28.png"),width=2.5,height=2,dpi=220)
    
    # plots for the 4-hour aggregated data
    c21 = data.frame(hour4 = as.POSIXct(ts2$hour4[1:6552],origin='1970-01-01',tz='UTC'),
              observed = c21$x, seasonal = c21$seasonal, trend = c21$trend,
              random = c21$random, predicted = c21$seasonal + c21$trend)
    c211 = melt(c21,id.vars = "hour4")
    p1 = ggplot(c211) + geom_line(aes(x=hour4,y=value)) + facet_grid(variable~.) + theme_bw()
    ggsave(paste0(dir2,"/additive_model_7day.pdf"),width=6,height=6)
    ggsave(paste0(dir2,"/additive_model_7day.png"),width=6,height=6,dpi = 150)
    c22 = data.frame(hour4 = as.POSIXct(ts2$hour4[1:6552],origin='1970-01-01',tz='UTC'),
              observed = c22$x, seasonal = c22$seasonal, trend = c22$trend,
              random = c22$random, predicted = c22$seasonal + c22$trend)
    c221 = melt(c22,id.vars = "hour4")
    p1 = ggplot(c221) + geom_line(aes(x=hour4,y=value)) + facet_grid(variable~.) + theme_bw()
    ggsave(paste0(dir2,"/additive_model_28day.pdf"),width=6,height=6)
    ggsave(paste0(dir2,"/additive_model_28day.png"),width=6,height=6,dpi = 150)
    
    # plot the seasonal TS
    p1 = ggplot(data.frame(hour4 = 1:42, activity = c21$seasonal[1:42])) + geom_line(aes(
      x=hour4,y=activity)) + theme_bw(base_size = 8) +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    # note: ts starts at 2007-01-01, which was a Monday
    p1 = p1 + scale_x_continuous(breaks=6*0:6+1,labels=c("Mon","Tue","Wed","Thu","Fri","Sat","Sun"))
    ggsave(paste0(dir2,"/additive_seasonal7.pdf"),width=2.5,height=2)
    ggsave(paste0(dir2,"/additive_seasonal7.png"),width=2.5,height=2,dpi=220)
    
    p1 = ggplot(data.frame(hour4 = 0:167/6, activity = c22$seasonal[1:168])) +
      geom_line(aes(x=hour4,y=activity)) + theme_bw(base_size = 8) + xlab("day") +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    ggsave(paste0(dir2,"/additive_seasonal28.pdf"),width=2.5,height=2)
    ggsave(paste0(dir2,"/additive_seasonal28.png"),width=2.5,height=2,dpi=220)
    
    
    ###################################################################
    ## 2.2 multiplicative decomposition with 7/28 day period
    m11 = mdecompose(ts1$cnt[1:1092],period = 7)
    names(m11)[1] = "day"
    m11$day2 = ts1$day[1:1092]
    m12 = mdecompose(ts1$cnt[1:1092],period = 28)
    names(m12)[1] = "day"
    m12$day2 = ts1$day[1:1092]
    m21 = mdecompose(ts2$cnt[1:6552],period = 42)
    m22 = mdecompose(ts2$cnt[1:6552],period = 168)
    m21$day2 = as.POSIXct(ts2$hour4[1:6552],origin='1970-01-01',tz='UTC')
    m22$day2 = as.POSIXct(ts2$hour4[1:6552],origin='1970-01-01',tz='UTC')
    
    # plots of the daily aggregated data
    p1 = ggplot(m11) + geom_line(aes(x=day2,y=cnt)) + facet_grid(type~.) + theme_bw()
    ggsave(paste0(dir1,"/multiplicative_model_7day.pdf"),width=6,height=6)
    ggsave(paste0(dir1,"/multiplicative_model_7day.png"),width=6,height=6,dpi = 150)
    p1 = ggplot(m12) + geom_line(aes(x=day2,y=cnt)) + facet_grid(type~.) + theme_bw()
    ggsave(paste0(dir1,"/multiplicative_model_28day.pdf"),width=6,height=6)
    ggsave(paste0(dir1,"/multiplicative_model_28day.png"),width=6,height=6,dpi = 150)
    
    # plot the seasonal TS separately
    p1 = ggplot(data.frame(day = 1:7, activity = m11[m11$type == "seasonal","cnt"][1:7])) +
      geom_line(aes(x=day,y=activity)) + theme_bw(base_size = 8) +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    # note: ts starts at 2007-01-01, which was a Monday
    p1 = p1 + scale_x_continuous(breaks=1:7,labels=c("Mon","Tue","Wed","Thu","Fri","Sat","Sun"))
    ggsave(paste0(dir1,"/multiplicative_seasonal7.pdf"),width=2.5,height=2)
    ggsave(paste0(dir1,"/multiplicative_seasonal7.png"),width=2.5,height=2,dpi=220)
    
    p1 = ggplot(data.frame(day = 1:28, activity = m12[m12$type == "seasonal","cnt"][1:28])) +
      geom_line(aes(x=day,y=activity)) + theme_bw(base_size = 8) +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    ggsave(paste0(dir1,"/multiplicative_seasonal28.pdf"),width=2.5,height=2)
    ggsave(paste0(dir1,"/multiplicative_seasonal28.png"),width=2.5,height=2,dpi=220)
    
    # plots for the 4-hour aggregated data
    p1 = ggplot(m21) + geom_line(aes(x=day2,y=cnt)) + facet_grid(type~.) + theme_bw()
    ggsave(paste0(dir2,"/multiplicative_model_7day.pdf"),width=6,height=6)
    ggsave(paste0(dir2,"/multiplicative_model_7day.png"),width=6,height=6,dpi = 150)
    p1 = ggplot(m22) + geom_line(aes(x=day2,y=cnt)) + facet_grid(type~.) + theme_bw()
    ggsave(paste0(dir2,"/multiplicative_model_28day.pdf"),width=6,height=6)
    ggsave(paste0(dir2,"/multiplicative_model_28day.png"),width=6,height=6,dpi = 150)
    
    # plot the seasonal TS separately
    p1 = ggplot(data.frame(day = 0:41/6+1, activity = m21[m21$type == "seasonal","cnt"][1:42])) +
      geom_line(aes(x=day,y=activity)) + theme_bw(base_size = 8) +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    # note: ts starts at 2007-01-01, which was a Monday
    p1 = p1 + scale_x_continuous(breaks=1:7,labels=c("Mon","Tue","Wed","Thu","Fri","Sat","Sun"))
    ggsave(paste0(dir2,"/multiplicative_seasonal7.pdf"),width=2.5,height=2)
    ggsave(paste0(dir2,"/multiplicative_seasonal7.png"),width=2.5,height=2,dpi=220)
    
    p1 = ggplot(data.frame(day = 0:167/6+1, activity = m22[m22$type == "seasonal","cnt"][1:168])) +
      geom_line(aes(x=day,y=activity)) + theme_bw(base_size = 8) +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
    ggsave(paste0(dir2,"/multiplicative_seasonal28.pdf"),width=2.5,height=2)
    ggsave(paste0(dir2,"/multiplicative_seasonal28.png"),width=2.5,height=2,dpi=220)
    
    
    # note: the decompositions are stored in c11,c12,c21,c22 and m11,m12,m21,m22
    #   for the additive and multiplicative models respectively
    #   first digit: 1: daily aggregation, 2: 4-hour aggregation (output should match this
    #   with dir1/dir2);
    #   second digit: 1: 7-day period, 2: 28-day period
    names(c111) = c("timestamp","type","cnt")
    names(c121) = c("timestamp","type","cnt")
    names(c211) = c("timestamp","type","cnt")
    names(c221) = c("timestamp","type","cnt")
    m111 = data.frame(timestamp = m11$day2, type = m11$type, cnt = m11$cnt)
    m121 = data.frame(timestamp = m12$day2, type = m12$type, cnt = m12$cnt)
    m211 = data.frame(timestamp = m21$day2, type = m21$type, cnt = m21$cnt)
    m221 = data.frame(timestamp = m22$day2, type = m22$type, cnt = m22$cnt)
    models = list(c111,c121,c211,c221,m111,m121,m211,m221)
    dirs = c(dir1,dir1,dir2,dir2,dir1,dir1,dir2,dir2)
    periods = c("7day","28day","7day","28day","7day","28day","7day","28day")
    periods2 = c(7,28,42,168,7,28,42,168)
    aggr1 = c(1,1,6,6,1,1,6,6)
    modeltypes = c("additive","additive","additive","additive","multiplicative",
                   "multiplicative","multiplicative","multiplicative")
    
    ###################################################################################
    # 3. temporal autocorrelation in the residuals, all cases
    for(i in 1:8) {
      fnbase = paste0(dirs[i],"/",modeltypes[i],"_acf_",periods[i])
      a1 = acf(na.omit(models[[i]][models[[i]]$type == "random","cnt"]),lag.max = aggr1[i]*370)
      a1 = data.frame(lag = a1$lag / aggr1[i], acf = a1$acf)
      p1 = ggplot(a1) + geom_line(aes(x=lag,y=acf)) + theme_bw(base_size=8) +
        xlab("lag (days)")
      ggsave(paste0(fnbase,".pdf"),p1,width=3.2,height=2)
      ggsave(paste0(fnbase,".png"),p1,width=3.2,height=2,dpi=220)
      p1 = p1 + scale_x_continuous(limits=c(0,7))
      ggsave(paste0(fnbase,"_1week.pdf"),p1,width=3.2,height=2)
      ggsave(paste0(fnbase,"_1week.png"),p1,width=3.2,height=2,dpi=220) 
    }
    
    
    #############################################################################
    # 4. distribution of noise / random component, test for normal distribution
    for(i in 1:8) {
      fnbase = paste0(dirs[i],"/",modeltypes[i],"_errdist_",periods[i])
      e1 = subset(na.omit(models[[i]]),type=="random","cnt")
      p1 = ggplot(e1) + geom_histogram(aes(x=cnt),binwidth=10,color="red",fill="red",alpha=0.2,size=0.3) +
        theme_bw(base_size = 8) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
        scale_x_continuous(limits=quantile(e1$cnt,c(0.01,0.99)))
      ggsave(paste0(fnbase,".pdf"),p1,width=3.55,height=2.35)
      ggsave(paste0(fnbase,".png"),p1,width=3.55,height=2.35,dpi=220)
      p1 = p1 + scale_y_continuous(trans="log",breaks=c(1,10,100,1000),limits=c(1,1000)) +
        scale_x_continuous(limits=quantile(e1$cnt,c(0,1)))
      ggsave(paste0(fnbase,"_log.pdf"),p1,width=3.55,height=2.35)
      ggsave(paste0(fnbase,"_log.png"),p1,width=3.55,height=2.35,dpi=220)
      
      fnbase = paste0(dirs[i],"/",modeltypes[i],"_normtest_",periods[i])
      nr1 = nrow(e1)
      v = data.frame(x = qnorm((1:nr1)/(nr1+1)), err = sort(e1$cnt))
      f1 = lm(err~x,v)
      a = f1$coefficients[2]
      b = f1$coefficients[1]
      
      p1 = ggplot(v) + geom_line(aes(x=x,y=err),color="red") + theme_bw(base_size = 8)
      p1 = p1 + geom_line(data = data.frame(x = c(-4,4),y=b + a*c(-4,4)),aes(x=x,y=y),color="black")
      p1 = p1 + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
      p1 = p1 + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
      ggsave(paste0(fnbase,".pdf"),p1,width=1.5,height=1)
      ggsave(paste0(fnbase,".png"),p1,width=1.5,height=1,dpi=220)
      
      fnbase = paste0(dirs[i],"/",modeltypes[i],"_kstest_",periods[i])
      v2 = data.frame(cdf = (0:(nr1-1))/(nr1-1),x = v$x,y = (v$err - b)/a)
      p2 = ggplot(v2) + geom_line(aes(x=x,y=cdf),color="blue") + geom_line(aes(x=y,y=cdf),color="red")
      p2 = p2 + theme_bw(base_size = 8) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
      p2 = p2 + scale_x_continuous(limits=c(-4,4))
      p2 = p2 + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
      ggsave(paste0(fnbase,".pdf"),p2,width=1.5,height=1)
      ggsave(paste0(fnbase,".png"),p2,width=1.5,height=1,dpi=220)
    }
    
    
    #####################################################################
    # 5. distribution of normalized residuals, average values for these
    for(i in 1:8) {
      fnbase = paste0(dirs[i],"/",modeltypes[i],"_epshist_",periods[i])
      tmp1 = na.omit(dcast(models[[i]],timestamp~type,value.var = "cnt"))
      e1 = as.numeric(na.omit(2*abs(tmp1$random) / (tmp1$observed + abs(tmp1$predicted))))
      p1 = ggplot(data.frame(err = e1)) + geom_histogram(aes(x=err),binwidth=0.05,
               color="red",fill="red",alpha=0.2,size=0.3)
      p1 = p1 + theme_bw(base_size = 8) + theme(panel.grid.major=element_blank(),
                                                panel.grid.minor=element_blank())
      avg1 = mean(e1)
      m1 = tmp1$observed
      avg2 = sum(m1*e1) / sum(m1)
      p1 = p1 + ggtitle(paste0("avg.: ",avg1,", weighted avg.: ",avg2))
      
      ggsave(paste0(fnbase,".pdf"),width=3.55,height=2.35)
      ggsave(paste0(fnbase,".png"),width=3.55,height=2.35,dpi=220)
    }
    
  }
}




#############################################################################
# more "basic" analysis:
# fit a lognormal distribution to the daily number of photos
# for New York (c = 0) it seems a good fit

dir.create("lndist")
setwd("lndist")

for(c in cities) {
  for(b in c("R","T")) {
    fn = paste0("../c_",c,"_",b,".txt")
    t1 = read.table(fn)
    if(b == "R")
      names(t1) = c("null","photoid","day","time","lat","lon")
    else
      names(t1) = c("null","photoid","day","time","lat","lon","country")
    
    ts1 = aggregate(t1["null"],by=t1["day"],FUN=length)
    
    t1 = data.frame(x=1:nrow(ts1),y=sort(ts1$null))
    mn = mean(log(ts1$null))
    s1 = sd(log(ts1$null))
    t1$x = t1$x/(nrow(ts1)+1)
    t1$x2 = plnorm(t1$y,mn,s1)
    p1 = ggplot(t1) + geom_line(aes(x=y,y=x),color="blue") + geom_line(aes(x=y,y=x2),color="red")
    p1 = p1 + theme_bw(base_size=8) + xlab("count") + ylab("CDF")
    fn1 = paste0("c_",c,"_",b,"_cdf")
    ggsave(paste0(fn1,".pdf"),width=3.55,height=2.35)
    ggsave(paste0(fn1,".png"),width=3.55,height=2.35,dpi=220)
    p1 = p1 + scale_x_continuous(trans="log",
        breaks=c(min(ts1$null),1,10,100,1000,max(ts1$null)),
        limits=c(min(ts1$null),max(ts1$null)))
    ggsave(paste0(fn1,"_log.pdf"),width=3.55,height=2.35)
    ggsave(paste0(fn1,"_log.png"),width=3.55,height=2.35,dpi=220)
    
    # quantile plot
    t1$y2 = qlnorm(t1$x,mn,s1)
    p1 = ggplot(t1) + geom_line(aes(x=y,y=y2)) + theme_bw(base_size=8) + xlab("counts") +
      ylab("distribution quantiles")
    
    # K-S distance
    t2 = ks.test(t1$y,plnorm,mn,s1)
    p1 = p1 + ggtitle(paste0("p = ",t2$p.value))
    
    fn1 = paste0("c_",c,"_",b,"_qq")
    ggsave(paste0(fn1,".pdf"),width=3.55,height=2.35)
    ggsave(paste0(fn1,".png"),width=3.55,height=2.35,dpi=220)
    p1 = p1 + scale_x_continuous(trans="log",
          breaks=c(min(ts1$null),1,10,100,1000,max(ts1$null)),
          limits=c(min(ts1$null),max(ts1$null))) + 
      scale_y_continuous(trans="log",
          breaks=c(min(t1$y2),1,10,100,1000,max(t1$y2)),
          limits=c(min(t1$y2),max(t1$y2)))
    ggsave(paste0(fn1,"_log.pdf"),width=3.55,height=2.35)
    ggsave(paste0(fn1,"_log.png"),width=3.55,height=2.35,dpi=220)
    
    print(paste0(c,"  --  ",b))
  }
}



##########################################################################
# similar analysis, print out a list of "outliers"
# simply take the logarithms, and select days with more than 2 sigma deviation
setwd("/home/dkondor/senseable/Data/flickr_timeseries/lndist/")

# problem: trend and weekly pattern
# try to "normalize" the daily values?
# i.e. totally multiplicative decomposition: A_t = s_t * d_t * u_t
setwd("/home/dkondor/senseable/Data/flickr_timeseries/")
dir.create("lndist2")
setwd("lndist2")


for(c in cities) {
  for(b in c("R","T")) {
    # 0. read file, aggregate, check for missing dates / times
    fn = paste0("../c_",c,"_",b,".txt")
    t1 = read.table(fn)
    # !! note: I replaced all the commas separating the coordinates with spaces in the files !!
    if(b == "R")
      names(t1) = c("null","photoid","day","time","lat","lon")
    else
      names(t1) = c("null","photoid","day","time","lat","lon","country")
    
    # daily aggregation
    ts1 = aggregate(t1["null"],by=t1["day"],FUN=length)
    ts1 = data.frame(day = as.POSIXct(ts1$day,tz="UTC"), cnt = ts1$null)
    
    # check for missing dates / times
    # all ts start 2007-01-01 and finish 2009-12-31 (incl.)
    ts1start = as.numeric(strptime("2007-01-01","%Y-%m-%d",tz="UTC"))
    ts1end = as.numeric(strptime("2009-12-31","%Y-%m-%d",tz="UTC"))
    n1 = (ts1end - ts1start) / 86400
    
    ts1t = data.frame(day = as.POSIXct(ts1start + 0:n1*86400, origin = '1970-01-01', tz = 'UTC'))
    
    ts1 = merge(ts1t,ts1,by="day",all.x=TRUE)
    
    ts1[is.na(ts1$cnt),"cnt"] = 0
    ts1 = ts1[with(ts1,order(day)),]
    
    ts11 = ts(ts1$cnt[1:1092], frequency = 7, start = c(1,1))
    c11 = decompose(ts11,"multiplicative")
    fn1 = paste0("c_",c,"_",b,"_mdc")
    pdf(paste0(fn1,".pdf"),width=6,height=6)
    plot(c11)
    dev.off()
    system(paste0("convert -density 150 ",fn1,".pdf -bordercolor white -border 0x0 ",
                  fn1,".png"))
    
    r1 = na.omit(c11$random)
    r1 = r1[r1 > 0]
    t1 = data.frame(x=1:length(r1),y=sort(r1))
    mn = mean(log(t1$y))
    s1 = sd(log(t1$y))
    
    # K-S distance
    t2 = ks.test(t1$y,plnorm,mn,s1)
    
    t1$x = t1$x/(nrow(t1)+1)
    t1$x2 = plnorm(t1$y,mn,s1)
    p1 = ggplot(t1) + geom_line(aes(x=y,y=x),color="blue") + geom_line(aes(x=y,y=x2),color="red")
    p1 = p1 + theme_bw(base_size=8) + xlab("random") + ylab("CDF")
    p1 = p1 + ggtitle(paste0("p = ",t2$p.value))
    fn1 = paste0("c_",c,"_",b,"_cdf")
    ggsave(paste0(fn1,".pdf"),width=3.55,height=2.35)
    ggsave(paste0(fn1,".png"),width=3.55,height=2.35,dpi=220)
    p1 = p1 + scale_x_continuous(trans="log")
    ggsave(paste0(fn1,"_log.pdf"),width=3.55,height=2.35)
    ggsave(paste0(fn1,"_log.png"),width=3.55,height=2.35,dpi=220)
    
    # quantile plot
    t1$y2 = qlnorm(t1$x,mn,s1)
    p1 = ggplot(t1) + geom_line(aes(x=y,y=y2)) + theme_bw(base_size=8) + xlab("random") +
      ylab("distribution quantiles")
    
    
    fn1 = paste0("c_",c,"_",b,"_qq")
    ggsave(paste0(fn1,".pdf"),width=3.55,height=2.35)
    ggsave(paste0(fn1,".png"),width=3.55,height=2.35,dpi=220)
    p1 = p1 + scale_x_continuous(trans="log") + 
      scale_y_continuous(trans="log")
    ggsave(paste0(fn1,"_log.pdf"),width=3.55,height=2.35)
    ggsave(paste0(fn1,"_log.png"),width=3.55,height=2.35,dpi=220)
    
  }
}

setwd('..')
dir.create('lndist2trend')
setwd('lndist2trend')

for(c in cities) {
  for(b in c("R","T")) {
    # 0. read file, aggregate, check for missing dates / times
    fn = paste0("../c_",c,"_",b,".txt")
    t1 = read.table(fn)
    # !! note: I replaced all the commas separating the coordinates with spaces in the files !!
    if(b == "R")
      names(t1) = c("null","photoid","day","time","lat","lon")
    else
      names(t1) = c("null","photoid","day","time","lat","lon","country")
    
    # daily aggregation
    ts1 = aggregate(t1["null"],by=t1["day"],FUN=length)
    ts1 = data.frame(day = as.POSIXct(ts1$day,tz="UTC"), cnt = ts1$null)
    
    # check for missing dates / times
    # all ts start 2007-01-01 and finish 2009-12-31 (incl.)
    ts1start = as.numeric(strptime("2007-01-01","%Y-%m-%d",tz="UTC"))
    ts1end = as.numeric(strptime("2009-12-31","%Y-%m-%d",tz="UTC"))
    n1 = (ts1end - ts1start) / 86400
    
    ts1t = data.frame(day = as.POSIXct(ts1start + 0:n1*86400, origin = '1970-01-01', tz = 'UTC'))
    
    ts1 = merge(ts1t,ts1,by="day",all.x=TRUE)
    
    ts1[is.na(ts1$cnt),"cnt"] = 0
    ts1 = ts1[with(ts1,order(day)),]
    
    ts11 = ts(ts1$cnt[1:1092], frequency = 7, start = c(1,1))
    c11 = decompose(ts11,"multiplicative")
    
    r1 = na.omit(c11$trend)
    r1 = r1[r1 > 0]
    t1 = data.frame(x=1:length(r1),y=sort(r1))
    mn = mean(log(t1$y))
    s1 = sd(log(t1$y))
    
    # K-S distance
    t2 = ks.test(t1$y,plnorm,mn,s1)
    
    t1$x = t1$x/(nrow(t1)+1)
    t1$x2 = plnorm(t1$y,mn,s1)
    p1 = ggplot(t1) + geom_line(aes(x=y,y=x),color="blue") + geom_line(aes(x=y,y=x2),color="red")
    p1 = p1 + theme_bw(base_size=8) + xlab("trend") + ylab("CDF")
    p1 = p1 + ggtitle(paste0("p = ",t2$p.value))
    fn1 = paste0("c_",c,"_",b,"_cdf")
    ggsave(paste0(fn1,".pdf"),width=3.55,height=2.35)
    ggsave(paste0(fn1,".png"),width=3.55,height=2.35,dpi=220)
    p1 = p1 + scale_x_continuous(trans="log")
    ggsave(paste0(fn1,"_log.pdf"),width=3.55,height=2.35)
    ggsave(paste0(fn1,"_log.png"),width=3.55,height=2.35,dpi=220)
    
    # quantile plot
    t1$y2 = qlnorm(t1$x,mn,s1)
    p1 = ggplot(t1) + geom_line(aes(x=y,y=y2)) + theme_bw(base_size=8) + xlab("trend") +
      ylab("distribution quantiles")
    
    
    fn1 = paste0("c_",c,"_",b,"_qq")
    ggsave(paste0(fn1,".pdf"),width=3.55,height=2.35)
    ggsave(paste0(fn1,".png"),width=3.55,height=2.35,dpi=220)
    p1 = p1 + scale_x_continuous(trans="log") + 
      scale_y_continuous(trans="log")
    ggsave(paste0(fn1,"_log.pdf"),width=3.55,height=2.35)
    ggsave(paste0(fn1,"_log.png"),width=3.55,height=2.35,dpi=220)
    
  }
}



#################################################
# write out a list of outlier days (top 1% for all cases)
setwd("~/senseable/Data/flickr_timeseries/")
dir.create("outliers")
setwd("outliers")
cities = c(0,1,2,3,4,5,7,8,10,46)
for(c in cities) {
  for(b in c("R","T")) {
    # 0. read file, aggregate, check for missing dates / times
    fn = paste0("../c_",c,"_",b,".txt")
    t1 = read.table(fn)
    # !! note: I replaced all the commas separating the coordinates with spaces in the files !!
    if(b == "R")
      names(t1) = c("null","photoid","day","time","lat","lon")
    else
      names(t1) = c("null","photoid","day","time","lat","lon","country")
    
    # daily aggregation
    ts1 = aggregate(t1["null"],by=t1["day"],FUN=length)
    ts1 = data.frame(day = as.POSIXct(ts1$day,tz="UTC"), cnt = ts1$null)
    
    # check for missing dates / times
    # all ts start 2007-01-01 and finish 2009-12-31 (incl.)
    ts1start = as.numeric(strptime("2007-01-01","%Y-%m-%d",tz="UTC"))
    ts1end = as.numeric(strptime("2009-12-31","%Y-%m-%d",tz="UTC"))
    n1 = (ts1end - ts1start) / 86400
    
    ts1t = data.frame(day = as.POSIXct(ts1start + 0:n1*86400, origin = '1970-01-01', tz = 'UTC'))
    
    ts1 = merge(ts1t,ts1,by="day",all.x=TRUE)
    
    ts1[is.na(ts1$cnt),"cnt"] = 0
    ts1 = ts1[with(ts1,order(day)),]
    
    ts11 = ts(ts1$cnt[1:1092], frequency = 7, start = c(1,1))
    c11 = decompose(ts11,"multiplicative")
    
    r1 = sort(ts1$cnt)
    n1 = length(r1)
    n2 = ceiling(n1*0.01)
    t1 = r1[n1-n2]
    days2 = ts1[ts1$cnt >= t1,"day"]
    fn = paste0("c_",c,"_",b,"_cnt.txt")
    write.table(days2,fn,row.names = FALSE,col.names = FALSE)
    
    r1 = sort(na.omit(c11$random))
    n1 = length(r1)
    n2 = ceiling(n1*0.01)
    t1 = r1[n1-n2]
    days2 = na.omit(ts1[c11$random >= t1,"day"])
    fn = paste0("c_",c,"_",b,"_random.txt")
    write.table(days2,fn,row.names = FALSE,col.names = FALSE)
    
    r1 = sort(na.omit(c11$trend))
    n1 = length(r1)
    n2 = ceiling(n1*0.01)
    t1 = r1[n1-n2]
    days2 = na.omit(ts1[c11$trend >= t1,"day"])
    fn = paste0("c_",c,"_",b,"_trend.txt")
    write.table(days2,fn,row.names = FALSE,col.names = FALSE)
  }
}


