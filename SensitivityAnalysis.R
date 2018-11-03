rm(list=ls())

#####################################
###### Variables to be changed ######
#####################################

# Function to be analyzed
func = function(x,y,z){
  return(sin(x)+7*sin(y)^2+0.1*z^4*sin(x))
}
# Fluctuation interval of input variables
bornes = data.frame(x=c(-pi,pi),y=c(-pi,pi),z=c(-pi,pi),row.names=c("min","max"))
# Number of samples of the experimental design
s = 2000
# Maximum order of sensitivity indices
max_order = 3
# Number of bootstrap samples
n_sample = 10000
# Tolerated error for the symmetric confidence interval of the indices by bootstrap
alpha = 0.10

#########################################
### Variables related to the analysis ###
#########################################

# Number of input variables
n = ncol(bornes)

# Detection of input errors
  # Checking the number of samples of the design of experiments
if(s <= 0){
  stop("The number of samples to be drawn must be strictly positive")
} else {
  s = as.integer(s)
  if(s!=round(s))
    stop("The number of samples to be drawn must be an integer")
}
  # Checking number of bootstrap samples
if(n_sample <= 0){
  stop("The number of bootstrap samples must be strictly positive")
} else {
  n_sample = as.integer(n_sample)
  if(n_sample!=round(n_sample))
    stop("The number of bootstrap samples must be an integer")
}
  # Checking maximum order of sensitivity indices
if(max_order!=round(max_order)){ 
  stop("The maximum order of sensitivity indices must be an integer")
} else {
  max_order = as.integer(max_order)
  if(max_order < 1 | max_order > n){
    stop("max_order must be an integer between 1 and ",n)
  } else {
    if(max_order %in% (n-1):(n-2)){
      max_order = n
      warning("The maximum order has been set at ",n," because it does not require any additional calculation\n")
    }
  }
}
  # Checking the function parameters
if(length(names(as.list(formals(func))))!=ncol(bornes)){
  stop("There must be as many columns in the data.frame bornes as there are parameters in the func function")
} else {
  if(!all(names(as.list(formals(func))) %in% names(bornes))){
    stop("The parameters of the func function must have the same name as the columns of the data.frame bornes")
  } else {
    if(!all(bornes[1,] < bornes[2,])){
      stop("The minimum of each variable must be less than its maximum")
    }
  }
}
  # Checking the alpha error
if(alpha<0|alpha>1){
  stop("alpha must be between 0 and 1")
}
    
# Downloading and loading the stringr library
libraries = c("stringr")
for(i in libraries){
  if(!(i %in% installed.packages()[,"Package"])) install.packages(i)
  library(i,character.only=T)
}
rm(libraries)

# List of matrices to consider
list_directories = c("a","b")
for(z in unique(c(1:max_order,n-1))){
  combinaisons = combn(1:n,z)
  for(i in 1:ncol(combinaisons)){
    x = NULL
    for(j in 1:z){
      if(j == 1)
        x = paste(x,'c',sep='')
      x = paste(x,combinaisons[j,i],sep='')
      if(j != z)
        x = paste(x,".",sep='')
    }
    list_directories = c(list_directories,x)
  }
  rm(combinaisons,i,j,z,x)
}
list_directories_without_b = list_directories[-2]

# LHS sampling function of a continuous uniform law
lhs_sampling_unifc = function(n,binf,bsup){
  v = NULL
  n = n+1
  for(q in (1:n)-1)
    v[q] = qunif(q/n+runif(1,0,1/n),binf,bsup)
  return(v)
}

# Converting a vector from time in seconds to time in D:H:M:S:MS
m_to_d = function(vec){
  n = length(vec)
  sortie = NULL
  for(i in 1:n){
    z = FALSE
    x = vec[i]
    str = ''
    d = floor(x / 1440)
    if(z | d > 0){
      x = x - d * 1440
      z = TRUE
      str = paste(str,d,'d',sep='')
    }
    h = floor(x / 60)
    if(z | h > 0){
      x = x - h * 60
      if(z & h < 10) h = paste(0,h,sep='')
      if(!z) z = T
      str = paste(str,h,'h',sep='')
    }
    m = floor(x)
    if(z | m > 0){
      x = x - m
      if(z & m < 10) m = paste(0,m,sep='')
      if(!z) z = T
      str = paste(str,m,'m',sep='')
    }
    s = floor(x * 60)
    if(z | s > 0){
      x = x - s / 60
      if(z & s < 10) s = paste(0,s,sep='')
      if(!z) z = T
      str = paste(str,s,'s',sep='')
    }
    ms = floor(x * 6000)
    if(z | ms >= 0){
      x = x - ms / 6000
      if(z & ms < 10) ms = paste(0,ms,sep='')
      if(!z) z = T
      str = paste(str,ms,'ms',sep='')
    }
    sortie[i] = str
    rm(str,x,d,h,m,s,ms)
  }
  return(sortie)
}

##########################################
### Launching the sensitivity analysis ###
##########################################

v_t = NULL # Vector of the cumulative times of the script progress
v_t[1] = 0

# Initialization
print("- Initialization")
dir.create("SensitivityAnalysis")
setwd(paste(getwd(),'/SensitivityAnalysis',sep='')) # Creating a working directory
t_start_script = Sys.time()
v_t[2] = as.numeric(difftime(Sys.time(),t_start_script,units="mins"))

# Creation of matrices A and B
print("- Creation of matrices A and B")
if(!dir.exists("matrices")) dir.create("matrices")
for(z in c("a.csv","b.csv")){
  d = matrix(nrow=s,ncol=n)
  for(i in 1:n){ # Sampling by LHS of each variable
    l_param = bornes[,i]
    d[,i] = sample(lhs_sampling_unifc(s,l_param[1],l_param[2]))
  }
  d = as.data.frame(d)
  names(d) = names(bornes)
  write.csv(d,paste('matrices/',z,sep=''),row.names=F) # Recording of experimental designs in .csv
  rm(d)
}
rm(z,l_param,i)
v_t[3] = as.numeric(difftime(Sys.time(),t_start_script,units="mins"))

# Creation of matrices C
print("- Creation of matrices C")
a = read.csv(paste("matrices/","a.csv",sep=''))
b = read.csv(paste("matrices/","b.csv",sep=''))
for(z in unique(c(1:max_order,n-1))){
  c = combn(1:n,z) # Creating factor combinations
  for(i in 1:ncol(c)){
    tmp = b
    tmp[,c[,i]] = a[,c[,i]]
    str = NULL
    for(j in 1:z){
      if(j == 1)
        str = paste(str,'c',sep='')
      str = paste(str,c[j,i],sep='')
      if(j != z)
        str = paste(str,".",sep='')
    }
    write.csv(tmp,paste("matrices/",str,".csv",sep=''),row.names=F)
  }
  rm(c,tmp,str,i,j)
}
rm(a,b,z)
v_t[4] = as.numeric(difftime(Sys.time(),t_start_script,units="mins"))

# Collection of the outputs Y
print("- Collection of the outputs Y")
if(!dir.exists("y")) dir.create("y")
for(i in list_directories_without_b){
  v = NULL
  for(j in 1:s)
   v[j] = do.call(func,unname(read.csv(paste("matrices/",i,".csv",sep=''))[j,]))
  write.csv(v,paste("y/","y",i,".csv",sep=''),row.names=F)
}
rm(i,j,v)
v_t[5] = as.numeric(difftime(Sys.time(),t_start_script,units="mins"))

# Calculation of sensitivity indices
print("- Calculation of sensitivity indices")
for(i in list_directories_without_b)
  assign(paste("y",i,sep=''),as.numeric(read.csv(paste('y/y',i,".csv",sep=''))[,1]))
rm(i)
eya = mean(ya)
vya = mean(ya^2)-eya^2
str = ''
for(z in 1:max_order){ # Calculation of z-th order sensitivity indices
  v_ind = NULL
  vars = apply(t(combn(1:n,z)),1,function(x) paste(x,collapse='.'))
  for(i in vars){
    v = get(paste("yc",i,sep=''))
    ind = (mean(ya*v)-eya*mean(v))/vya
    list_mat = NULL
    if(z > 1){
      for(j in 1:(z-1))
        list_mat[length(list_mat)+1] = list(apply(t(combn(as.integer(unlist(str_split(i,'[.]'))),j)),1,function(x) paste(x,collapse='.')))
      list_mat = unlist(list_mat)
      for(j in list_mat)
        ind = ind - get(paste("s",j,sep="")) 
    }
    assign(paste("s",i,sep=''),ind)
    v_ind[length(v_ind)+1] = ind
    names(v_ind)[length(v_ind)] = i
  }
  rm(i,ind,vars,v,list_mat)
  str = str_c(str,paste("Sensivity indice",ifelse(z==max_order,'','s')," of ",z,ifelse(z%%10==1,"st",ifelse(z==2,"nd",ifelse(z==3,"rd","th")))," order :\n",sep=''))
  #v_sorted = sort(v_ind,decreasing=T)
  v_sorted = v_ind
  lmax = max(unlist(lapply(lapply(lapply(strsplit(names(v_sorted),'[.]'),as.integer),function(x) paste(names(bornes)[x],collapse=',')),nchar)))
  names_mat = names(v_sorted)
  for(i in 1:length(v_sorted)){
    name = paste(names(bornes)[as.integer(unlist(strsplit(names_mat[i],"[.]")))],collapse=",")
    ncar = nchar(name)
    str = str_c(str,name,paste(rep(' ',lmax-ncar),collapse='')," : ",ifelse(v_sorted[i]<0,'',' '),signif(unname(v_sorted[i]),4),'\n',sep='')
  }
  rm(i,v_sorted,ncar,name,names_mat,lmax)
  str = str_c(str,'\n')
  assign(paste("sum_",z,sep=''),sum(v_ind))
}
rm(z,v_ind,j)
for(z in 1:max_order)
  str = str_c(str,"Sum of sensitivity indices of ",z,ifelse(z%%10==1,"st",ifelse(z==2,"nd",ifelse(z==3,"rd","th")))," order : ",ifelse(get(paste("sum_",z,sep=''))<0,'',' '),signif(get(paste("sum_",z,sep='')),4),'\n',sep='')
rm(z)
str = str_c(str,'\n')
# Calculation of total sensitivity indices
str = str_c(str,'Calculation by complementarity','\n')
for(z in 1:n){
  v = get(paste("yc",paste((1:n)[-z],collapse='.'),sep=''))
  ind = 1-(mean(ya*v)-eya*mean(v))/vya
  assign(paste("st",z,sep=''),ind)
  str = str_c(str,"Total sensitivity index of ",names(bornes)[z],
              paste(rep(" ",max(apply(t(t(names(bornes))),1,nchar))-apply(t(t(names(bornes))),1,nchar)[z]),collapse='')," : ",ifelse(ind<0,'',' '),signif(ind,4),'\n',sep='')
}
rm(ind,z,v)
str = str_c(str,'\n')
if(max_order == n){
  str = str_c(str,'Calculation by sum','\n')
  v=NULL
  for(z in 1:n){
    ss=0
    for(i in ls(pattern=paste("^s[0-9.*]*[",z,"][0-9.*]*",sep=''))) 
      ss=ss+get(i)
    v[z] = ss
    assign(paste("sts",z,sep=''),ss)
  }
  v = signif(v,4)
  for(z in 1:n)
    str = str_c(str,"Total sensitivity index of ",names(bornes)[z],paste(rep(' ',max(unlist(lapply(t(t(names(bornes))),nchar)))-nchar(names(bornes)[z])),collapse='')," : ",ifelse(v[z]<0,'',' '),v[z],'\n')
}
rm(z,v)
str = str_c(str,'\n')
v_t[6] = as.numeric(difftime(Sys.time(),t_start_script,units="mins"))

# Estimation of confidence intervals by bootstrap
print("- Estimation of confidence intervals by bootstrap")
list_directories_c = paste('c',1:n,sep='')
for(i in list_directories_c){
  m = NULL
  mt = NULL
  for(j in 1:n_sample){# Resampling j
    ind = sample(1:s,s,replace=T) # Resampled indices
    y = get(paste("y",i,sep=''))[ind] # Ci
    ym = get(paste("yc",paste((1:n)[-as.numeric(substr(i,2,nchar(i)))],collapse='.'),sep=''))[ind] # C-i
    ya2 = ya[ind] # A
    rm(ind)
    # Calculation of new indicators
    eya2 = mean(ya2) # Mean of A
    vya2 = mean(ya2^2)-eya2^2 # Variance of A
    m[j] = (mean(ya2*y)-eya2*mean(y))/vya2 # Si
    mt[j] = 1-(mean(ya2*ym)-eya2*mean(ym))/vya2 # STi
    rm(ya2,eya2,vya2,y,ym)
  }
  assign(paste("m",i,sep=''),mean(m))
  assign(paste("sd",i,sep=''),sd(m))
  assign(paste("mt",i,sep=''),mean(mt))
  assign(paste("sdt",i,sep=''),sd(mt))
}
rm(m,mt,i,j)
sd_vec = NULL;sdt_vec = NULL;m_vec = NULL;mt_vec = NULL
for(i in 1:n){
  sd_vec[i] = get(paste("sdc",i,sep=''))
  sdt_vec[i] = get(paste("sdtc",i,sep=''))
  m_vec[i] = get(paste("mc",i,sep=''))
  mt_vec[i] = get(paste("mtc",i,sep=''))
}
rm(i,eya,vya)
rm(list=paste(c("sdc","sdtc","mc","mtc"),rep(1:n,4),sep=''))
v_t[7] = as.numeric(difftime(Sys.time(),t_start_script,units="mins"))

# Creation of the output graph
print("- Creation of the output graph")
v_1 = NULL
v_tot = NULL
for(i in 1:n) {
  v_1[i] = get(paste("s",i,sep=''))
  v_tot[i] = get(paste("st",i,sep=''))
}
c = qt(1-alpha/2,s-1)
di = rbind(v_1,v_tot-v_1)
val_s1 = round(v_1*100,2)
val_st = round(v_tot*100,2)
d = di
tmp = d
d = matrix(0,nrow=length(di),ncol=ncol(di))
d[cbind(1:length(di),rep(1:ncol(di),each=2))] = tmp
png(filename=paste("output.jpg",sep=''),width=600,height=600)
default = par(no.readonly=T)
layout(matrix(c(1,2),ncol=1,byrow=TRUE),heights=c(20,2))
barp = barplot(d,ylim=c(0,min(max(v_tot,v_1)*1.35,1)),main="Sensitivity indices",sub=paste((1-alpha)*100,"% confidence intervals",sep=''),col=unlist(ifelse(di[2,]<0,list(c("gray80","gray40")),list(c("gray40","gray80")))))
abline(h=0)
text(barp,par("usr")[3],labels=names(bornes),adj=c(0,3),cex=.9,xpd=T)
text(barp,-0.035,labels=paste(val_s1,"%",sep=''),adj=c(0.5,-0.5),cex=0.8,xpd=T,font=2)
text(barp,ifelse(apply(rbind(v_1,v_tot),2,max)>0,apply(rbind(v_1,v_tot),2,max),0),labels=paste(val_st,"%",sep=''),adj=c(0.5,-0.3),cex=0.8,xpd=T,font=2)
segments(barp,mt_vec-sdt_vec*c,barp,mt_vec+sdt_vec*c,lwd=1.5,col='red',xpd=T)
arrows(barp,mt_vec-sdt_vec*c,barp,mt_vec+sdt_vec*c,lwd=1.5,angle=90,code=3,length=0.07,col='red',xpd=T)
segments(barp,m_vec-sd_vec*c,barp,m_vec+sd_vec*c,lwd=1.5,col='green',xpd=T)
arrows(barp,m_vec-sd_vec*c,barp,m_vec+sd_vec*c,lwd=1.5,angle=90,code=3,length=0.07,col='green',xpd=T)
par(mai=c(0,0,0,0))
plot.new()
legend("center",ncol=2,legend=c("Total indices","First-order indices","CI of total indices","CI of first-order indices"),col=c("gray80","gray40","red","green"),lty=c(NA,NA,1,1),pch=c(15,15,NA,NA),xpd=T,bty="n")
par(default)
dev.off()
v_t[8] = as.numeric(difftime(Sys.time(),t_start_script,units="mins"))

# Storage of the info.txt file
print("- Storage of the info.txt file")
sink("infos.txt")
cat("------------------------------------------------------------\n")
cat("------------------- Sensitivity analysis -------------------\n")
cat("------------------------------------------------------------\n")
cat('\n')
cat("------------------- Function parameters --------------------\n")
cat('\n')
print(func)
cat('\n')
cat("------------- Sensitivity analysis parameters --------------\n")
cat('\n')
cat("Number of samples        :",s,"\n")
cat("Maximum order of indices :",max_order,"\n")
cat("Number of factors        :",n,"\n")
for(i in 1:n)
    cat("-> ",names(bornes)[i],", sampled by LHS according to uniform continuous law of parameters min = ",bornes[1,i]," and max = ",bornes[2,i],"\n",sep='')
cat('\n')
cat("----------- Results of the sensitivity analysis ------------\n")
cat('\n')
cat(str)
cat("------------- Cumulative script execution time -------------\n")
cat('\n')
v_t[9] = as.numeric(difftime(Sys.time(),t_start_script,units="mins"))
v_t = m_to_d(v_t)
names(v_t) = c("Initialization","Creation of matrices A and B",
               "Creation of matrices C","Collection of the outputs Y",
               "Calculation of sensitivity indices",
               "Estimation of confidence intervals by bootstrap",
               "Creation of the output graph","Storage of the info.txt file","End")
lbeforemax = max(unname(sapply(names(v_t),nchar)))
laftermax = max(unname(sapply(v_t,nchar)))
for(i in 1:length(v_t))
  cat(names(v_t)[i],rep(' ',lbeforemax-nchar(names(v_t)[i])+1),':',rep(' ',laftermax-nchar(v_t[i])+1),unname(v_t[i]),'\n',sep='')
sink()
setwd("..")