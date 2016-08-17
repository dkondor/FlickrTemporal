# flickr_ts3.r -- new processing steps for Flickr time series data
# "missing": relative weights of components in the power spectrum
#   normalized residuals

library(ggplot2)
library(reshape2)

setwd("~/senseable/Data/flickr_timeseries/")
cities = c(0,1,2,3,4,5,7,8,10,46)
cities2 = data.frame(ids = cities, names = c("New York","London","Paris","San Francisco","Berlin","Washington","Barcelona","Rome","Chicago","Los Angeles"))


# steps:
# 0. read file, aggregation (daily), check for missing dates, fill these with
#   either zero or NAs
# 1. Fourier-spectrum, plot of "selected" frequencies
# 2. decomposition (additive / multiplicative noise, 7-day period)
# 3. average "predictability score" / relative noise

dir.create("run3")
setwd("run3")

colors1 = c("#5aa732", "#e75216", "#009ace", "#ffd500", "#e78502", "#db6aa2", "#007939", "#8671a7", "#005c96", "#815234", "#9a2846")

# perform a multiplicative model (i.e. x = seasonal*trend + error) decomposition of a 
#   given time series (ts1 should be given as a vector)
# !!notes: this might not work if length(ts1) is not an exact multiple of the period,
#   (I have not tested that configuration)
# also, this might be a bit slow, as the code is not optimal (uses loops to go thorugh
#   the data as my R skills are still a bit limited in some aspects)
mdecompose <- function(ts1, period) {
  len = length(ts1)
  ws = period
  ts5 = data.frame(day = 1:len, cnt = ts1)
  n1 = aggregate(ts5["cnt"], by = data.frame(x = floor((ts5$day-1)/ws)), FUN = sum)
  n12 = ts5
  for(i in 1:len) {
    n12$cnt[i] = n12$cnt[i] / n1$cnt[floor((i-1)/ws)+1]
  }
  a1 = aggregate(n12$cnt[1:len], by = data.frame(d = 0:(len-1)%%ws), FUN = mean)
  pr1 = data.frame(day = 1:len, cnt = 1:len)
  for(i in 1:len) {
    pr1$cnt[i] = n1$cnt[floor((i-1)/ws+1)] * a1$x[ (i-1)%%ws + 1 ]
  }
  
  t1 = data.frame(day = 1:len, cnt = 1:len)
  for(i in 1:len) {
    t1$cnt[i] = n1$cnt[floor((i-1)/ws)+1] / ws
  }
  s1 = data.frame(day = 1:len, cnt = 1:len)
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
  
  e1 = merge(ts5[1:len,],pr1,by="day")
  e1$err = e1$cnt.x - e1$cnt.y
  tsall1 = rbind(tsall1, data.frame(day = e1$day, cnt = e1$err, type = "random"))
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
    
    ts1 = aggregate(t1["null"],by=t1["day"],FUN=length)
    ts1 = data.frame(day = as.POSIXct(ts1$day,tz="UTC"), cnt = ts1$null)
    ts1start = as.numeric(strptime("2007-01-01","%Y-%m-%d",tz="UTC"))
    ts1end = as.numeric(strptime("2009-12-31","%Y-%m-%d",tz="UTC"))
    n1 = (ts1end - ts1start) / 86400
    
    ts1t = data.frame(day = as.POSIXct(ts1start + 0:n1*86400, origin = '1970-01-01', tz = 'UTC'))
    
    ts1 = merge(ts1t,ts1,by="day",all.x=TRUE)
    
    ts1[is.na(ts1$cnt),"cnt"] = 0
    ts1 = ts1[with(ts1,order(day)),]
    
    # Fourier-transform, power spectrum, figs. of the relative weight of most relevant
    #   periods
    fft1 = fft(ts1$cnt[1:1092])
    tmp12 = Re(fft1[1:546]*Conj(fft1[1:546]))
    tmp13 = data.frame(f = 0:545, w = tmp12)
    
#     p1 = ggplot(tmp13) + geom_line(aes(x=f,y=w),size=0.35) + scale_y_continuous(trans="log") +
#       scale_x_continuous(trans="log", breaks = c(3,12,36.4,156,312),
#                          labels = c("1 year","1/4 year","1 month","1 week","1/2 week"))
#     p1 = p1 + theme_bw(base_size=7) + theme(axis.text.y=element_blank()) + xlab("period") + 
#       ylab("weight (arbitrary units, log scale)")
#     head(tmp13[order(tmp13$w,decreasing = TRUE),],n = 10)
#     periods = data.frame(name=c("1 year","1/2 year","1/4 year","1 month","1 week","1/2 week",
#                                 "19.5 days","64.2 days"),
#                          f = c(3,6,12,36.4,156,312,56,17))
#     periods2 = sort(c(156,312,3,468,30,159,42,6,12))
#     1092/periods2
#     
#     periods = data.frame(f = sort(c(156,312,3,468,30,159,42,6,12)),
#           name = c("1 year","1/2 year","1/4 year","36.4 days","26 days","1 week",
#                    "6.87 days","1/2 week","1/3 week"))
    
    tmp13$w2 = tmp13$w / sum(tmp13$w)
    
    # write out 20 main components to check for later
#     tmp15 = head(tmp13[order(tmp13$w,decreasing = TRUE),],n = 20)
#     tmp15$d = 1092/tmp15$f
#     write.table(tmp15,paste0("fcomponents1_",c,"_",b,".dat"),quote=FALSE,row.names=FALSE,sep="\t")
    
    # plot first 11 components (expect 0,1,2) to see what to include later
    tmp15 = tmp13[ tmp13$f != 0, ]
    tmp15 = head(tmp15[order(tmp15$w,decreasing = TRUE),],n=11)
    tmp15$wmin = 1e-3
    
    p1 = ggplot(tmp15) + geom_linerange(aes(x=factor(f),ymin=wmin,ymax=w2,color=factor(f)),size=4)
    p1 = p1 + theme_bw(4) + scale_y_log10(limits=c(1e-3,1e-1),breaks=c(0.001,0.01,0.1))
    p1 = p1 + guides(color=FALSE) + scale_x_discrete(breaks = factor(tmp15$f),
              labels = format(1092/tmp15$f,digits=4,drop0trailing = TRUE,trim=TRUE))
    p1 = p1 + scale_color_manual(values=colors1) + xlab("period (days)") + ylab("rel. weight")
    p1 = p1 + theme(panel.grid = element_blank(), axis.ticks = element_line(size=0.25),
                    plot.background = element_blank())
    if(b=="R") p1 = p1 + ggtitle(paste0(cities2[cities2$ids == c,"names"],", residents"))
    else p1 = p1 + ggtitle(paste0(cities2[cities2$ids == c,"names"],", visitors"))
    fn0 = paste0("fourier_weights2_",c,"_",b)
    ggsave(paste0(fn0,".pdf"),p1,width=2.8,height=1.9)
    system(paste0("convert -density 300 ",fn0,".pdf -bordercolor white -border 0x0 -units PixelsPerInch -density 300 ",fn0,".png"))
 
     
# new run: all periods among the first five for any TS
#     periods = data.frame(f = sort(c(156,1,312,3,2,62,468,42,4,153,14,87,85,74,71,64,6,55,49,33,31,30,29,26,24,222,22,183,159,155,15,147,12)),
#         name = c("3 years", "1.5 year", "1 year","3/4 year","1/2 year","1/4 year",
#           "78 days (13 weeks)","72.8 days","49.63 days","45.5 days","42 days (6 weeks)",
#           "37.655 days","36.4 days (5.2 weeks)","35.2258 days","33.1 days","26 days",
#           "22.286 days","19.855 days","17.6129 days","17.0625 days","15.38 days",
#           "14.757 days","12.847 days","12.552 days","7.429 days","7.1373 days","7.045 days",
#           "1 week","6.87 days","5.9672 days","4.919 days","1/2 week","1/3 week"))
#     
#     
#     tmp14 = merge(tmp13,periods,by="f")
#     # + write out the "selected" components to adjust figures
#     # write.table(tmp14,paste0("fcomponents2_",c,"_",b,".dat"),quote=FALSE,row.names=FALSE,sep="\t")
#     
#     tmp14$wmin = 1e-5
#     p1 = ggplot(tmp14) + geom_linerange(aes(x=factor(f),ymin=wmin,ymax=w2,color=factor(f)),size=3)
#     p1 = p1 + theme_bw(4) + scale_y_log10(limits=c(1e-5,1e-1),breaks=c(1e-5,0.0001,0.001,0.01,0.1))
#     p1 = p1 + scale_x_discrete(breaks = factor(tmp14$f),labels=tmp14$name)
#     p1 = p1 + guides(color=FALSE) + theme(axis.text.x = element_text(angle=90,hjust=1))
#     p1 = p1 + scale_color_manual(values=c(colors1,colors1,colors1)) + xlab("period") + ylab("rel. weight")
#     p1 = p1 + theme(panel.grid = element_blank(), axis.ticks = element_line(size=0.25),
#                     plot.background = element_blank())
#     if(b=="R") p1 = p1 + ggtitle(paste0(cities2[cities2$ids == c,"names"],", residents"))
#     else p1 = p1 + ggtitle(paste0(cities2[cities2$ids == c,"names"],", visitors"))
#     fn0 = paste0("fourier_weights3_",c,"_",b)
#     ggsave(paste0(fn0,".pdf"),p1,width=5,height=2.5)
#     system(paste0("convert -density 300 ",fn0,".pdf -bordercolor white -border 0x0 -units PixelsPerInch -density 300 ",fn0,".png"))
  }
}


errors = data.frame(c = cities,rmerr = 1:10,rmstdev=1:10,raerr=1:10,rastdev=1:10,
                      tmerr = 1:10,tmstdev=1:10,taerr=1:10,tastdev=1:10)

# "predictability scores" / normalized residuals
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
    
    ts1 = aggregate(t1["null"],by=t1["day"],FUN=length)
    ts1 = data.frame(day = as.POSIXct(ts1$day,tz="UTC"), cnt = ts1$null)
    ts1start = as.numeric(strptime("2007-01-01","%Y-%m-%d",tz="UTC"))
    ts1end = as.numeric(strptime("2009-12-31","%Y-%m-%d",tz="UTC"))
    n1 = (ts1end - ts1start) / 86400
    
    ts1t = data.frame(day = as.POSIXct(ts1start + 0:n1*86400, origin = '1970-01-01', tz = 'UTC'))
    
    ts1 = merge(ts1t,ts1,by="day",all.x=TRUE)
    
    ts1[is.na(ts1$cnt),"cnt"] = 0
    ts1 = ts1[with(ts1,order(day)),]
    
    # 1. additive decomposition (multiplicative previously): 
    # measured = trend*seasonal + noise
    m11 = mdecompose(ts1$cnt[1:1092],period = 7)
    
    tmp1 = dcast(m11,day~type,value.var = "cnt")
    # sum(abs(tmp1$observed - tmp1$predicted - tmp1$random)) # 0, OK
    e1 = 2*abs(tmp1$random) / (tmp1$observed + tmp1$predicted)
    
    # 2. multiplicative decomposition:
    # measured = trend*seasonal*noise
    ts11 = ts(ts1$cnt[1:1092], frequency = 7, start = c(1,1))
    c11 = decompose(ts11,"multiplicative")
    e2 = abs(na.omit(c11$random)-1)
    
    i = which(errors$c == c)
    if(b == "R") {
      errors$raerr[i] = mean(e1)
      errors$rastdev[i] = sd(e1)
      errors$rmerr[i] = mean(e2)
      errors$rmstdev[i] = sd(e2)
    }
    else {
      errors$taerr[i] = mean(e1)
      errors$tastdev[i] = sd(e1)
      errors$tmerr[i] = mean(e2)
      errors$tmstdev[i] = sd(e2)
    }
  }
}


# relative noise as a function of total number of photos
errors$rcnt = c(641690,622243,332287,505379,113815,154223,305717,64334,270979,170845)
errors$tcnt = c(212069,199723,171821,157499,86558,92907,96282,87268,86851,86684)

p1 = ggplot(errors) + geom_line(aes(x=sqrt(rcnt),y=raerr),color="red") +
  geom_line(aes(x=sqrt(rcnt),y=rmerr),color="blue")


###############################################
# plot periogram together:
#   gray background + most important frequencies as colored dots
# 1. load all data, save power spectrum

frequencies = list()

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
    
    ts1 = aggregate(t1["null"],by=t1["day"],FUN=length)
    ts1 = data.frame(day = as.POSIXct(ts1$day,tz="UTC"), cnt = ts1$null)
    ts1start = as.numeric(strptime("2007-01-01","%Y-%m-%d",tz="UTC"))
    ts1end = as.numeric(strptime("2009-12-31","%Y-%m-%d",tz="UTC"))
    n1 = (ts1end - ts1start) / 86400
    
    ts1t = data.frame(day = as.POSIXct(ts1start + 0:n1*86400, origin = '1970-01-01', tz = 'UTC'))
    
    ts1 = merge(ts1t,ts1,by="day",all.x=TRUE)
    
    ts1[is.na(ts1$cnt),"cnt"] = 0
    ts1 = ts1[with(ts1,order(day)),]
    
    # Fourier-transform, power spectrum, figs. of the relative weight of most relevant
    #   periods
    fft1 = fft(ts1$cnt[1:1092])
    tmp12 = Re(fft1[1:546]*Conj(fft1[1:546]))
    tmp13 = data.frame(f = 0:545, w = tmp12)
    
    i = which(cities2$ids == c)
    if(b == "T") i = i + 10
    frequencies[[i]] = tmp13
  }
}


# example fig
p1 = ggplot(frequencies[[13]]) + geom_line(aes(x=f,y=w),size=0.35) + scale_y_continuous(trans="log") +
  scale_x_continuous(trans="log", breaks = c(3,12,36.4,156,312),
                     labels = c("1 year","1/4 year","1 month","1 week","1/2 week"))
p1 = p1 + theme_bw(base_size=7) + theme(axis.text.y=element_blank()) + xlab("period") + 
  ylab("weight (arbitrary units, log scale)")

# problem: DC components, omit this for all, normalize
for(i in 1:20) {
  frequencies[[i]] = frequencies[[i]][2:546,]
  frequencies[[i]]$wn = frequencies[[i]]$w / sum(frequencies[[i]]$w)
}


# residents
# collection of most important frequencies (first 10 for all?)
f10r1 = data.frame(f = integer())
for(i in 1:10) {
  f10r2 = data.frame(f = head(frequencies[[i]][order(frequencies[[i]]$w,decreasing = TRUE),],n = 20)$f)
  f10r1 = rbind(f10r1,f10r2)
}
f10r = data.frame(f = unique(f10r1$f))

frr = frequencies[[1]]
frr$i = 1
for(i in 2:10) {
  frr = rbind(frr,data.frame(frequencies[[i]],i=i))
}
p1 = ggplot(frr) + geom_line(aes(x=f,y=wn,group=i),size=0.35,color="grey")
p1 = p1 + scale_y_log10() + scale_x_continuous(trans="log", breaks = c(1,3,12,36.4,156,312),
            labels = c("3 years","1 year","1/4 year","1 month","1 week","1/2 week"))
p1 = p1 + theme_bw(8) + xlab("period") + ylab("relative weight")

f10r2 = merge(f10r,frequencies[[1]],by="f")
f10r2$i = 1
for(i in 2:10) {
  f10r2 = rbind(f10r2,data.frame(merge(f10r,frequencies[[i]],by="f"),i=i))
}
f10r2$i = factor(f10r2$i)
p2 = p1 + geom_point(data=f10r2, aes(x=f,y=wn,color=i)) + scale_color_manual(values=colors1,
        breaks=1:10,labels=cities2$names,name="city")
ggsave("fourier_all1.pdf",p2,width=6,height=4)
ggsave("fourier_all1.png",p2,width=6,height=4,dpi=300)

f20r = f10r
f20r2 = f10r2

# too many points (obviously)
# maybe display just the first ten? (previous: first 20)
f10r1 = data.frame(f = integer())
for(i in 1:10) {
  f10r2 = data.frame(f = head(frequencies[[i]][order(frequencies[[i]]$w,decreasing = TRUE),],n = 10)$f)
  f10r1 = rbind(f10r1,f10r2)
}
f10r = data.frame(f = unique(f10r1$f))

f10r2 = merge(f10r,frequencies[[1]],by="f")
f10r2$i = 1
for(i in 2:10) {
  f10r2 = rbind(f10r2,data.frame(merge(f10r,frequencies[[i]],by="f"),i=i))
}
f10r2$i = factor(f10r2$i)

p2 = p1 + geom_point(data=f10r2, aes(x=f,y=wn,color=i)) + scale_color_manual(values=colors1,
               breaks=1:10,labels=cities2$names,name="city")
ggsave("fourier_all10.pdf",p2,width=6,height=4)
ggsave("fourier_all10.png",p2,width=6,height=4,dpi=300)


# first 5
f5r1 = data.frame(f = integer())
for(i in 1:10) {
  f5r2 = data.frame(f = head(frequencies[[i]][order(frequencies[[i]]$w,decreasing = TRUE),],n = 5)$f)
  f5r1 = rbind(f5r1,f5r2)
}
f5r = data.frame(f = unique(f5r1$f))

f5r2 = merge(f5r,frequencies[[1]],by="f")
f5r2$i = 1
for(i in 2:10) {
  f5r2 = rbind(f5r2,data.frame(merge(f5r,frequencies[[i]],by="f"),i=i))
}
f5r2$i = factor(f5r2$i)

p2 = p1 + geom_point(data=f5r2, aes(x=f,y=wn,color=i)) + scale_color_manual(values=colors1,
             breaks=1:10,labels=cities2$names,name="city")
ggsave("fourier_all5.pdf",p2,width=6,height=4)
ggsave("fourier_all5.png",p2,width=6,height=4,dpi=300)


# do the same for tourists
frt = frequencies[[11]]
frt$i = 1
for(i in 12:20) {
  frt = rbind(frt,data.frame(frequencies[[i]],i=i-10))
}
p1 = ggplot(frt) + geom_line(aes(x=f,y=wn,group=i),size=0.35,color="grey")
p1 = p1 + scale_y_log10() + scale_x_continuous(trans="log", breaks = c(1,3,12,36.4,156,312),
                                               labels = c("3 years","1 year","1/4 year","1 month","1 week","1/2 week"))
p1 = p1 + theme_bw(8) + xlab("period") + ylab("relative weight")

# collection of most important frequencies (first 10 for all?)
f20t1 = data.frame(f = integer())
for(i in 11:20) {
  f20t2 = data.frame(f = head(frequencies[[i]][order(frequencies[[i]]$w,decreasing = TRUE),],n = 20)$f)
  f20t1 = rbind(f20t1,f20t2)
}
f20t = data.frame(f = unique(f20t1$f))

f20t2 = merge(f20t,frequencies[[11]],by="f")
f20t2$i = 1
for(i in 12:20) {
  f20t2 = rbind(f20t2,data.frame(merge(f20t,frequencies[[i]],by="f"),i=i-10))
}
f20t2$i = factor(f20t2$i)
p2 = p1 + geom_point(data=f20t2, aes(x=f,y=wn,color=i)) + scale_color_manual(values=colors1,
      breaks=1:10,labels=cities2$names,name="city")
ggsave("fourier_allt20.pdf",p2,width=6,height=4)
ggsave("fourier_allt20.png",p2,width=6,height=4,dpi=300)


# too many points (obviously)
# maybe display just the first ten? (previous: first 20)
f10t1 = data.frame(f = integer())
for(i in 11:20) {
  f10t2 = data.frame(f = head(frequencies[[i]][order(frequencies[[i]]$w,decreasing = TRUE),],n = 10)$f)
  f10t1 = rbind(f10t1,f10t2)
}
f10t = data.frame(f = unique(f10t1$f))

f10t2 = merge(f10t,frequencies[[11]],by="f")
f10t2$i = 1
for(i in 12:20) {
  f10t2 = rbind(f10t2,data.frame(merge(f10t,frequencies[[i]],by="f"),i=i-10))
}
f10t2$i = factor(f10t2$i)

p2 = p1 + geom_point(data=f10t2, aes(x=f,y=wn,color=i)) + scale_color_manual(values=colors1,
            breaks=1:10,labels=cities2$names,name="city")
ggsave("fourier_allt10.pdf",p2,width=6,height=4)
ggsave("fourier_allt10.png",p2,width=6,height=4,dpi=300)


# first 5
f5t1 = data.frame(f = integer())
for(i in 11:20) {
  f5t2 = data.frame(f = head(frequencies[[i]][order(frequencies[[i]]$w,decreasing = TRUE),],n = 5)$f)
  f5t1 = rbind(f5t1,f5t2)
}
f5t = data.frame(f = unique(f5t1$f))

f5t2 = merge(f5t,frequencies[[1]],by="f")
f5t2$i = 1
for(i in 12:20) {
  f5t2 = rbind(f5t2,data.frame(merge(f5t,frequencies[[i]],by="f"),i=i-10))
}
f5t2$i = factor(f5t2$i)

p2 = p1 + geom_point(data=f5t2, aes(x=f,y=wn,color=i)) + scale_color_manual(values=colors1,
        breaks=1:10,labels=cities2$names,name="city")
ggsave("fourier_allt5.pdf",p2,width=6,height=4)
ggsave("fourier_allt5.png",p2,width=6,height=4,dpi=300)


################################################
# re-do individual periodogram plots for all cities
for(i in 1:10) {
  c = cities2$ids[i]
  for(b in c("R","T")) {
    # 0. read file, aggregate, check for missing dates / times
    fn = paste0("../c_",c,"_",b,".txt")
    t1 = read.table(fn)
    # !! note: I replaced all the commas separating the coordinates with spaces in the files !!
    b2 = NA
    if(b == "R") {
      names(t1) = c("null","photoid","day","time","lat","lon")
      b2 = "residents"
    }
    else {
      names(t1) = c("null","photoid","day","time","lat","lon","country")
      b2 = "visitors"
    }
    
    ts1 = aggregate(t1["null"],by=t1["day"],FUN=length)
    ts1 = data.frame(day = as.POSIXct(ts1$day,tz="UTC"), cnt = ts1$null)
    ts1start = as.numeric(strptime("2007-01-01","%Y-%m-%d",tz="UTC"))
    ts1end = as.numeric(strptime("2009-12-31","%Y-%m-%d",tz="UTC"))
    n1 = (ts1end - ts1start) / 86400
    
    ts1t = data.frame(day = as.POSIXct(ts1start + 0:n1*86400, origin = '1970-01-01', tz = 'UTC'))
    
    ts1 = merge(ts1t,ts1,by="day",all.x=TRUE)
    
    ts1[is.na(ts1$cnt),"cnt"] = 0
    ts1 = ts1[with(ts1,order(day)),]
    
    # Fourier-transform, power spectrum, figs. of the relative weight of most relevant
    #   periods
    fft1 = fft(ts1$cnt[1:1092])
    tmp12 = Re(fft1[1:546]*Conj(fft1[1:546]))
    tmp13 = data.frame(f = 0:545, w = tmp12)
    tmp13$w2 = tmp13$w / sum(tmp13$w)
    
    # leave out DC component
    tmp13 = tmp13[tmp13$f != 0,]
    
    p1 = ggplot(tmp13) + geom_line(aes(x=f,y=w2),size=0.35) + scale_y_log10() +
      scale_x_continuous(trans="log", breaks = c(1,3,12,36.4,156,312),
            labels = c("3 years","1 year","1/4 year","1 month","1 week","1/2 week"))
    p1 = p1 + theme_bw(5) + xlab("period") + ylab("rel. weight") +
      ggtitle(paste0(cities2$names[i],", ",b2))
    ggsave(paste0("fourier_spectrum_",c,"_",b,".pdf"),p1,width=3,height=2)
    
  }
}



