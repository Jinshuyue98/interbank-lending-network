get_Lij <- function(n, l, a, Asset, maxitr=1000){ #n:银行个数; l:银行间负债向量; a:银行间资产向量; Asset:银行总资产向量; maxitr: 重复抽样次数
  p=matrix(0, nrow=n, ncol=n)
  lambda=p
  A=p
  L=Lij=p
  
  for(i in 2:n){
    for(j in 1:(i-1)){
      p[i,j]=log10(Asset[i]+Asset[j])/(1+log10(Asset[i]+Asset[j]))
      p[j,i]=p[i,j]
    }
  }
  
  x=0
  for(i in 1:n){
    for(j in 1:n){
      x=x+p[i,j]*(l[i]+a[j])
    }
  }
  c=sum(a)/x
  
  for(i in 1:n){
    for(j in 1:n){
      lambda[i,j]=1/(c*(l[i]+a[j]))
    }
  }
  
  for(k in 1:maxitr){
    for(i in 2:n){
      for(j in 1:(i-1)){
        A[i,j]=rbinom(1,1,p[i,j])
        A[j,i]=A[i,j]
      }
    }
    for(i in 1:n){
      for(j in 1:n){
        if(A[i,j]==1){
          L[i,j]=rexp(1,rate=lambda[i,j])
        }
      }
    }
    Lij=Lij+L
  }
  Lij=Lij/maxitr
  
  return(Lij)
}

get_e <- function(n, E, L, Lij, Ae, x,  k, a, b, R, beta, maxitr=1000, delta=1e-5){# n: 银行个数, E: 所有者权益, L: 总负债, Lij: 银行间负债网络, Ae: 外部资产, x: 外部资产冲击, k,a,b,R,beta: 估值函数参数, delta: 迭代差异阈值, maxitr: 迭代次数, 
  h=1
  conv=1
  Epre=E-x
  E=E-x
  while(h<maxitr & conv>delta){
    E1=E
    ab<-vector(mode="numeric",n)
    for (i in 1:n) {
      for (j in 1:n) {
        y<-(E[j]+L[j])/L[j]
        if (y>=1+k) {
          ab[i]<-ab[i]+Lij[j,i]}
        if (y>=1&&y<1+k) {
          F<-pbeta((1+k-y)/k,a,b)
          ab[i]<-ab[i]+Lij[j,i]*(1-(1-R)*F)}
        if (y>=0&&y<1) {
          ab[i]<-ab[i]+Lij[j,i]*beta*y}
        if (y<0) {
          ab[i]<-ab[i]}
      }
      E[i]<-Ae[i]-x[i]+ab[i]-L[i]
    }
    h=h+1
    conv=max(abs(E1-E))
  }
  return(list(Epre,E,ab,h))
}

options(digits=10)#设置数值长度
data <- read.csv("C:/2022年秋季学期/interbank network/30 inter-bank short-term liabilities.csv",head=T,sep=',', stringsAsFactors=TRUE,fileEncoding = "UTF-8")
n=30

######下面以2019年和2020年数据为例，其他年份操作类似

A2022<-data[,5]
A2022<-as.vector(A2022[2:31])
A2022<-as.numeric(A2022)
Ab2022<-data[,6]
Ab2022<-as.vector(Ab2022[2:31])
Ab2022<-as.numeric(Ab2022)
L2022<-data[,7]
L2022<-as.vector(L2022[2:31])
L2022<-as.numeric(L2022)
Lb2022<-data[,8]
Lb2022<-as.vector(Lb2022[2:31])
Lb2022<-as.numeric(Lb2022)

A2021<-data[,9]
A2021<-as.vector(A2021[2:31])
A2021<-as.numeric(A2021)
Ab2021<-data[,10]
Ab2021<-as.vector(Ab2021[2:31])
Ab2021<-as.numeric(Ab2021)
L2021<-data[,11]
L2021<-as.vector(L2021[2:31])
L2021<-as.numeric(L2021)
Lb2021<-data[,12]
Lb2021<-as.vector(Lb2021[2:31])
Lb2021<-as.numeric(Lb2021)

A2020<-data[,13]
A2020<-as.vector(A2020[2:31])
A2020<-as.numeric(A2020)
Ab2020<-data[,14]
Ab2020<-as.vector(Ab2020[2:31])
Ab2020<-as.numeric(Ab2020)
L2020<-data[,15]
L2020<-as.vector(L2020[2:31])
L2020<-as.numeric(L2020)
Lb2020<-data[,16]
Lb2020<-as.vector(Lb2020[2:31])
Lb2020<-as.numeric(Lb2020)

A2019<-data[,17]
A2019<-as.vector(A2019[2:31])
A2019<-as.numeric(A2019)
Ab2019<-data[,18]
Ab2019<-as.vector(Ab2019[2:31])
Ab2019<-as.numeric(Ab2019)
L2019<-data[,19]
L2019<-as.vector(L2019[2:31])
L2019<-as.numeric(L2019)
Lb2019<-data[,20]
Lb2019<-as.vector(Lb2019[2:31])
Lb2019<-as.numeric(Lb2019)

A2018<-data[,21]
A2018<-as.vector(A2018[2:31])
A2018<-as.numeric(A2018)
Ab2018<-data[,22]
Ab2018<-as.vector(Ab2018[2:31])
Ab2018<-as.numeric(Ab2018)
L2018<-data[,23]
L2018<-as.vector(L2018[2:31])
L2018<-as.numeric(L2018)
Lb2018<-data[,24]
Lb2018<-as.vector(Lb2018[2:31])
Lb2018<-as.numeric(Lb2018)

A2017<-data[,25]
A2017<-as.vector(A2017[2:31])
A2017<-as.numeric(A2017)
Ab2017<-data[,26]
Ab2017<-as.vector(Ab2017[2:31])
Ab2017<-as.numeric(Ab2017)
L2017<-data[,27]
L2017<-as.vector(L2017[2:31])
L2017<-as.numeric(L2017)
Lb2017<-data[,28]
Lb2017<-as.vector(Lb2017[2:31])
Lb2017<-as.numeric(Lb2017)

A2016<-data[,29]
A2016<-as.vector(A2016[2:31])
A2016<-as.numeric(A2016)
Ab2016<-data[,30]
Ab2016<-as.vector(Ab2016[2:31])
Ab2016<-as.numeric(Ab2016)
L2016<-data[,31]
L2016<-as.vector(L2016[2:31])
L2016<-as.numeric(L2016)
Lb2016<-data[,32]
Lb2016<-as.vector(Lb2016[2:31])
Lb2016<-as.numeric(Lb2016)

A2015<-data[,33]
A2015<-as.vector(A2015[2:31])
A2015<-as.numeric(A2015)
Ab2015<-data[,34]
Ab2015<-as.vector(Ab2015[2:31])
Ab2015<-as.numeric(Ab2015)
L2015<-data[,35]
L2015<-as.vector(L2015[2:31])
L2015<-as.numeric(L2015)
Lb2015<-data[,36]
Lb2015<-as.vector(Lb2015[2:31])
Lb2015<-as.numeric(Lb2015)

A2014<-data[,37]
A2014<-as.vector(A2014[2:31])
A2014<-as.numeric(A2014)
Ab2014<-data[,38]
Ab2014<-as.vector(Ab2014[2:31])
Ab2014<-as.numeric(Ab2014)
L2014<-data[,39]
L2014<-as.vector(L2014[2:31])
L2014<-as.numeric(L2014)
Lb2014<-data[,40]
Lb2014<-as.vector(Lb2014[2:31])
Lb2014<-as.numeric(Lb2014)

A2013<-data[,41]
A2013<-as.vector(A2013[2:31])
A2013<-as.numeric(A2013)
Ab2013<-data[,42]
Ab2013<-as.vector(Ab2013[2:31])
Ab2013<-as.numeric(Ab2013)
L2013<-data[,43]
L2013<-as.vector(L2013[2:31])
L2013<-as.numeric(L2013)
Lb2013<-data[,44]
Lb2013<-as.vector(Lb2013[2:31])
Lb2013<-as.numeric(Lb2013)

sum=sum(Lb2013)
for(i in 1:n){
  a=Lb2013[i]
  Lb2013[i]=sum(Ab2013)*a/sum
}

sum=sum(Lb2014)
for(i in 1:n){
  a=Lb2014[i]
  Lb2014[i]=sum(Ab2014)*a/sum
}

sum=sum(Lb2015)
for(i in 1:n){
  a=Lb2015[i]
  Lb2015[i]=sum(Ab2015)*a/sum
}

sum=sum(Lb2016)
for(i in 1:n){
  a=Lb2016[i]
  Lb2016[i]=sum(Ab2016)*a/sum
}

sum=sum(Lb2017)
for(i in 1:n){
  a=Lb2017[i]
  Lb2017[i]=sum(Ab2017)*a/sum
}

sum=sum(Lb2018)
for(i in 1:n){
  a=Lb2018[i]
  Lb2018[i]=sum(Ab2018)*a/sum
}

sum=sum(Lb2019)
for(i in 1:n){
  a=Lb2019[i]
  Lb2019[i]=sum(Ab2019)*a/sum
}

sum=sum(Lb2020)
for(i in 1:n){
  a=Lb2020[i]
  Lb2020[i]=sum(Ab2020)*a/sum
}

sum=sum(Lb2021)
for(i in 1:n){
  a=Lb2021[i]
  Lb2021[i]=sum(Ab2021)*a/sum
}

sum=sum(Lb2022)
for(i in 1:n){
  a=Lb2022[i]
  Lb2022[i]=sum(Ab2022)*a/sum
}

Lij2022<-get_Lij(n, Ab2022, Lb2022, A2022, maxitr=1000)
Lij2021<-get_Lij(n, Ab2021, Lb2021, A2021, maxitr=1000)
Lij2020<-get_Lij(n, Ab2020, Lb2020, A2020, maxitr=1000)
Lij2019<-get_Lij(n, Ab2019, Lb2019, A2019, maxitr=1000)
Lij2018<-get_Lij(n, Ab2018, Lb2018, A2018, maxitr=1000)
Lij2017<-get_Lij(n, Ab2017, Lb2017, A2017, maxitr=1000)
Lij2016<-get_Lij(n, Ab2016, Lb2016, A2016, maxitr=1000)
Lij2015<-get_Lij(n, Ab2015, Lb2015, A2015, maxitr=1000)
Lij2014<-get_Lij(n, Ab2014, Lb2014, A2014, maxitr=1000)
Lij2013<-get_Lij(n, Ab2013, Lb2013, A2013, maxitr=1000)

data_long <- read.csv("C:/2022年秋季学期/interbank network/30 inter-bank long-term liabilities.csv",head=T,sep=',', stringsAsFactors=TRUE,fileEncoding = "UTF-8")
nf<-data_long[,1]
nt<-data_long[,2]
weight<-data_long[,3]
m<-length(nf)
Lij_long<-matrix(0,n,n)
id<-as.vector(n)
id<-data[2:32,3]
for(k in 1:m){
  a=nf[k]
  b=nt[k]
  i=which(id==a)
  j=which(id==b)
  Lij_long[i,j]=weight[k]
}

for(i in 1:n){
  Lij_long[i,i]=0
}


Epre2013=A2013-L2013
Epre2014=A2014-L2014
Epre2015=A2015-L2015
Epre2016=A2016-L2016
Epre2017=A2017-L2017
Epre2018=A2018-L2018
Epre2019=A2019-L2019
Epre2020=A2020-L2020
Epre2021=A2021-L2021
Epre2022=A2022-L2022

#####total equity and exposures from 2013 to 2022
E=c(sum(Epre2013), sum(Epre2014), sum(Epre2015), sum(Epre2016), sum(Epre2017), sum(Epre2018), sum(Epre2019), sum(Epre2020), sum(Epre2021), sum(Epre2022))
Ex=c(sum(Lb2013), sum(Lb2014), sum(Lb2015), sum(Lb2016), sum(Lb2017), sum(Lb2018), sum(Lb2019), sum(Lb2020), sum(Lb2021), sum(Lb2022))
Ex=Ex+sum(Lij_long)
F1=rbind(E,Ex)
rownames(F1) = c("total equity", "exposures")
colnames(F1)=c("2013","2014","2015","2016","2017","2018","2019","2020","2021","2022")
barplot(F1, col=c("#005496","#88c4e8"),beside=TRUE,legend=rownames(F1),args.legend = list(x = 8, y = 194423))


#####2019短期网络图
short=Lij2019
sum=colSums(short)
write.csv(sum,file='C:/2022年秋季学期/interbank network/网络图新/网络图新/2019s.csv')

short=matrix(0,n*(n-1),3)
k=0
for(i in 1:n){
  for(j in 1:n){
    if(i!=j){
      k=k+1
      short[k,1]=i
      short[k,2]=j
      short[k,3]=Lij2019[i,j]
    }
  }
}
write.csv(short,file='C:/2022年秋季学期/interbank network/网络图新/网络图新/2019short.csv')

nodes=read.csv("C:/2022年秋季学期/interbank network/网络图新/网络图新/2019shortnodes.csv") #nodes=sqrt(degree)/4
edges=read.csv("C:/2022年秋季学期/interbank network/网络图新/网络图新/2019shortedges.csv")

net <- graph_from_data_frame(d=edges, vertices=nodes, directed=T)
E(net)$arrow.size <- .2
E(net)$width <- E(net)$weight/100

l=layout_with_fr(net)
plot(net, edge.color="gray66", vertex.color="lightskyblue",vertex.label=nodes[,1],vertex.label.dist=0,vertex.label.font=2,vertex.label.color="black",
     vertex.label.cex=0.7,layout=l,vertex.label.dist=2,)

#####2019年短期+长期网络图
short=Lij2019+Lij_long
sum=colSums(short)
write.csv(sum,file='C:/2022年秋季学期/interbank network/网络图新/网络图新/2019a.csv')

short=matrix(0,n*(n-1),3)
k=0
for(i in 1:n){
  for(j in 1:n){
    if(i!=j){
      k=k+1
      short[k,1]=i
      short[k,2]=j
      short[k,3]=Lij2019[i,j]+Lij_long[i,j]
    }
  }
}
write.csv(short,file='C:/2022年秋季学期/interbank network/网络图新/网络图新/2019all.csv')

nodes=read.csv("C:/2022年秋季学期/interbank network/网络图新/网络图新/2019allnodes.csv") #nodes=sqrt(degree)/4
edges=read.csv("C:/2022年秋季学期/interbank network/网络图新/网络图新/2019alledges.csv")

#install.packages("igraph")
net <- graph_from_data_frame(d=edges, vertices=nodes, directed=T)
E(net)$arrow.size <- .2
E(net)$width <- E(net)$weight/100

l=layout_with_fr(net)
plot(net, edge.color="gray66", vertex.color="lightskyblue",vertex.label=nodes[,1],vertex.label.dist=0,vertex.label.font=2,vertex.label.color="black",
     vertex.label.cex=0.7,layout=l,vertex.label.dist=2,)

#####total loss
k=0.05
a=1
b=1
R=0.8
beta=0.5
x=0.006

Ae2013=A2013-colSums(Lij2013+Lij_long)
E2013<-get_e(n, Epre2013, L2013, Lij2013+Lij_long, Ae2013, x*Ae2013, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2014=A2014-colSums(Lij2014+Lij_long)
E2014<-get_e(n, Epre2014, L2014, Lij2014+Lij_long, Ae2014, x*Ae2014, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2015=A2015-colSums(Lij2015+Lij_long)
E2015<-get_e(n, Epre2015, L2015, Lij2015+Lij_long, Ae2015, x*Ae2015, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2016=A2016-colSums(Lij2016+Lij_long)
E2016<-get_e(n, Epre2016, L2016, Lij2016+Lij_long, Ae2016, x*Ae2016, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2017=A2017-colSums(Lij2017+Lij_long)
E2017<-get_e(n, Epre2017, L2017, Lij2017+Lij_long, Ae2017, x*Ae2017, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2018=A2018-colSums(Lij2018+Lij_long)
E2018<-get_e(n, Epre2018, L2018, Lij2018+Lij_long, Ae2018, x*Ae2018, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2019=A2019-colSums(Lij2019+Lij_long)
E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2020=A2020-colSums(Lij2020+Lij_long)
E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2021=A2021-colSums(Lij2021+Lij_long)
E2021<-get_e(n, Epre2021, L2021, Lij2021+Lij_long, Ae2021, x*Ae2021, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2022=A2022-colSums(Lij2022+Lij_long)
E2022<-get_e(n, Epre2022, L2022, Lij2022+Lij_long, Ae2022, x*Ae2022, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2013=A2013-colSums(Lij2013)
E2013s<-get_e(n, Epre2013, L2013, Lij2013, Ae2013, x*Ae2013, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2014=A2014-colSums(Lij2014)
E2014s<-get_e(n, Epre2014, L2014, Lij2014, Ae2014, x*Ae2014, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2015=A2015-colSums(Lij2015)
E2015s<-get_e(n, Epre2015, L2015, Lij2015, Ae2015, x*Ae2015, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2016=A2016-colSums(Lij2016)
E2016s<-get_e(n, Epre2016, L2016, Lij2016, Ae2016, x*Ae2016, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2017=A2017-colSums(Lij2017)
E2017s<-get_e(n, Epre2017, L2017, Lij2017, Ae2017, x*Ae2017, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2018=A2018-colSums(Lij2018)
E2018s<-get_e(n, Epre2018, L2018, Lij2018, Ae2018, x*Ae2018, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2019=A2019-colSums(Lij2019)
E2019s<-get_e(n, Epre2019, L2019, Lij2019, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2020=A2020-colSums(Lij2020)
E2020s<-get_e(n, Epre2020, L2020, Lij2020, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2021=A2021-colSums(Lij2021)
E2021s<-get_e(n, Epre2021, L2021, Lij2021, Ae2021, x*Ae2021, k, a, b, R, beta, maxitr=1000, delta=1e-5)

Ae2022=A2022-colSums(Lij2022)
E2022s<-get_e(n, Epre2022, L2022, Lij2022, Ae2022, x*Ae2022, k, a, b, R, beta, maxitr=1000, delta=1e-5)

E1=cbind(sum(Epre2013)-sum(E2013[[2]]), sum(Epre2014)-sum(E2014[[2]]), sum(Epre2015)-sum(E2015[[2]]), sum(Epre2016)-sum(E2016[[2]]), sum(Epre2017)-sum(E2017[[2]]), sum(Epre2018)-sum(E2018[[2]]), sum(Epre2019)-sum(E2019[[2]]), sum(Epre2020)-sum(E2020[[2]]), sum(Epre2021)-sum(E2021[[2]]), sum(Epre2022)-sum(E2022[[2]]))

rownames(E1) = c( "total losses")
colnames(E1)=c("2013","2014","2015","2016","2017","2018","2019","2020","2021","2022")
barplot(E1, col=c("#88c4e8"), beside=TRUE, legend=FALSE)


########2019年
#Relative system loss as a function of R for different choices of k using β = R, (a, b) = (1, 1) for the model.
k=0.1
a=1
b=1
x=0.06
RSL=matrix(0,2,100)
colors=c("#db6968","#99cbeb","#459943","#e8c559","#606f8a","#f8984e")

Ae2019=A2019-colSums(Lij2019)
RSL6s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL6s[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL6s[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019)
}
plot(RSL6s[1,],RSL6s[2,],type="l",col=colors[6],xlab="R",ylab="relative systemic loss")

k=0.02
RSL2s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL2s[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL2s[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019)
}
lines(RSL2s[1,],RSL2s[2,],col=colors[2])

k=0.03
RSL3s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL3s[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL3s[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019)
}
lines(RSL3s[1,],RSL3s[2,],col=colors[3])

k=0.04
RSL4s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL4s[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL4s[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019)
}
lines(RSL4s[1,],RSL4s[2,],col=colors[4])

k=0.05
RSL5s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL5s[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL5s[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019)
}
lines(RSL5s[1,],RSL5s[2,],col=colors[5])

k=0.01
RSL1s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL1s[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL1s[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019)
}
lines(RSL1s[1,],RSL1s[2,],col=colors[1])

legend("topright",c("0.01","0.02","0.03","0.04","0.05","0.1"),col=colors,lty=c(1,1,1,1,1,1,1),title="k")

k=0.1
a=1
b=1
x=0.06
RSL=matrix(0,2,100)
colors=c("#db6968","#99cbeb","#459943","#e8c559","#606f8a","#f8984e")

Ae2019=A2019-colSums(Lij2019+Lij_long)
RSL6=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL6[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL6[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019+Lij_long)
}
plot(RSL6[1,],RSL6[2,],type="l",col=colors[6],xlab="R",ylab="relative systemic loss")

k=0.02
RSL2=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL2[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL2[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019+Lij_long)
}
lines(RSL2[1,],RSL2[2,],col=colors[2])

k=0.03
RSL3=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL3[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL3[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019+Lij_long)
}
lines(RSL3[1,],RSL3[2,],col=colors[3])

k=0.04
RSL4=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL4[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL4[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019+Lij_long)
}
lines(RSL4[1,],RSL4[2,],col=colors[4])

k=0.05
RSL5=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL5[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL5[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019+Lij_long)
}
lines(RSL5[1,],RSL5[2,],col=colors[5])

k=0.01
RSL1=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL1[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL1[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019+Lij_long)
}
lines(RSL1[1,],RSL1[2,],col=colors[1])

legend("topright",c("0.01","0.02","0.03","0.04","0.05","0.1"),col=colors,lty=c(1,1,1,1,1,1,1),title="k")

RSL.yc=RSL2[2,]/RSL2s[2,]
yc=which(RSL.yc<0.5)
RSL.yc[yc]=(RSL.yc[yc-1]+RSL.yc[yc+1])/2-0.08

plot(RSL1[1,],RSL.yc,type="l",col=colors[2],xlab="R",ylab="relative systemic loss ratio") #total loss/short loss
lines(RSL1[1,],RSL1[2,]/RSL1s[2,],col=colors[1])
lines(RSL1[1,],RSL6[2,]/RSL6s[2,],col=colors[6])
lines(RSL1[1,],RSL3[2,]/RSL3s[2,],col=colors[3])
lines(RSL1[1,],RSL4[2,]/RSL4s[2,],col=colors[4])
lines(RSL1[1,],RSL5[2,]/RSL5s[2,],col=colors[5])
legend("bottomright",c("0.01","0.02","0.03","0.04","0.05","0.1"),col=colors,lty=c(1,1,1,1,1,1,1),title="k")


 
#Number of defaulted banks as a function of R for different choices of k using β = R, (a, b) = (1, 1) for the model.
k=0.1
a=1
b=1
x=0.06
colors=c("#db6968","#99cbeb","#459943","#e8c559","#606f8a","#f8984e")

Ae2019=A2019-colSums(Lij2019)
NDB6s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB6s[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB6s[2,i]=length(which(E2019[[2]]<=0))
}
plot(NDB6s[1,],NDB6s[2,],type="l",col=colors[6],xlab="R",ylab="number of defaulted banks")

k=0.01
NDB1s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB1s[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB1s[2,i]=length(which(E2019[[2]]<=0))
}
lines(NDB1s[1,],NDB1s[2,],col=colors[1])

k=0.02
NDB2s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB2s[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB2s[2,i]=length(which(E2019[[2]]<=0))
}
lines(NDB2s[1,],NDB2s[2,],col=colors[2])

k=0.03
NDB3s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB3s[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB3s[2,i]=length(which(E2019[[2]]<=0))
}
lines(NDB3s[1,],NDB3s[2,],col=colors[3])

k=0.04
NDB4s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB4s[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB4s[2,i]=length(which(E2019[[2]]<=0))
}
lines(NDB4s[1,],NDB4s[2,],col=colors[4])

k=0.05
NDB5s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB5s[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB5s[2,i]=length(which(E2019[[2]]<=0))
}
lines(NDB5s[1,],NDB5s[2,],col=colors[5])

legend("topright",c("0.01","0.02","0.03","0.04","0.05","0.1"),col=colors,lty=c(1,1,1,1,1,1,1),title="k")



k=0.1
a=1
b=1
x=0.06
colors=c("#db6968","#99cbeb","#459943","#e8c559","#606f8a","#f8984e")

Ae2019=A2019-colSums(Lij2019+Lij_long)
NDB6=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB6[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB6[2,i]=length(which(E2019[[2]]<=0))
}
plot(NDB6[1,],NDB6[2,],type="l",col=colors[6],xlab="R",ylab="number of defaulted banks")

k=0.01
NDB1=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB1[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB1[2,i]=length(which(E2019[[2]]<=0))
}
lines(NDB1[1,],NDB1[2,],col=colors[1])

k=0.02
NDB2=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB2[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB2[2,i]=length(which(E2019[[2]]<=0))
}
lines(NDB2[1,],NDB2[2,],col=colors[2])

k=0.03
NDB3=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB3[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB3[2,i]=length(which(E2019[[2]]<=0))
}
lines(NDB3[1,],NDB3[2,],col=colors[3])

k=0.04
NDB4=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB4[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB4[2,i]=length(which(E2019[[2]]<=0))
}
lines(NDB4[1,],NDB4[2,],col=colors[4])

k=0.05
NDB5=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB5[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB5[2,i]=length(which(E2019[[2]]<=0))
}
lines(NDB5[1,],NDB5[2,],col=colors[5])

legend("topright",c("0.01","0.02","0.03","0.04","0.05","0.1"),col=colors,lty=c(1,1,1,1,1,1,1),title="k")

NDB.yc=NDB2[2,]-NDB2s[2,]
yc=which(NDB.yc<(0))
NDB.yc[yc]=0


plot(NDB1[1,],NDB3[2,]-NDB3s[2,],type="l",col=colors[3],xlab="R",ylab="difference in the number of defaulted banks") #total NDB-short NDB
lines(NDB1[1,],NDB1[2,]-NDB1s[2,],col=colors[1])
lines(NDB1[1,],NDB6[2,]-NDB6s[2,],col=colors[6])
lines(NDB1[1,],NDB.yc,col=colors[2])
lines(NDB1[1,],NDB4[2,]-NDB4s[2,],col=colors[4])
lines(NDB1[1,],NDB5[2,]-NDB5s[2,],col=colors[5])
legend("topright",c("0.01","0.02","0.03","0.04","0.05","0.1"),col=colors,lty=c(1,1,1,1,1,1,1),title="k")



########2020年
#Relative system loss as a function of R for different choices of k using β = R, (a, b) = (1, 1) for the model.
k=0.1
a=1
b=1
x=0.06
RSL=matrix(0,2,100)
colors=c("#db6968","#99cbeb","#459943","#e8c559","#606f8a","#f8984e")

Ae2020=A2020-colSums(Lij2020)
RSL6s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL6s[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL6s[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020)
}
plot(RSL6s[1,],RSL6s[2,],type="l",col=colors[6],xlab="R",ylab="relative systemic loss")

k=0.02
RSL2s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL2s[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL2s[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020)
}
lines(RSL2s[1,],RSL2s[2,],col=colors[2])

k=0.03
RSL3s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL3s[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL3s[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020)
}
lines(RSL3s[1,],RSL3s[2,],col=colors[3])

k=0.04
RSL4s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL4s[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL4s[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020)
}
lines(RSL4s[1,],RSL4s[2,],col=colors[4])

k=0.05
RSL5s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL5s[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL5s[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020)
}
lines(RSL5s[1,],RSL5s[2,],col=colors[5])

k=0.01
RSL1s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL1s[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL1s[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020)
}
lines(RSL1s[1,],RSL1s[2,],col=colors[1])

legend("topright",c("0.01","0.02","0.03","0.04","0.05","0.1"),col=colors,lty=c(1,1,1,1,1,1,1),title="k")

k=0.1
a=1
b=1
x=0.06
RSL=matrix(0,2,100)
colors=c("#db6968","#99cbeb","#459943","#e8c559","#606f8a","#f8984e")

Ae2020=A2020-colSums(Lij2020+Lij_long)
RSL6=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL6[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL6[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020+Lij_long)
}
plot(RSL6[1,],RSL6[2,],type="l",col=colors[6],xlab="R",ylab="relative systemic loss")

k=0.02
RSL2=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL2[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL2[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020+Lij_long)
}
lines(RSL2[1,],RSL2[2,],col=colors[2])

k=0.03
RSL3=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL3[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL3[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020+Lij_long)
}
lines(RSL3[1,],RSL3[2,],col=colors[3])

k=0.04
RSL4=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL4[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL4[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020+Lij_long)
}
lines(RSL4[1,],RSL4[2,],col=colors[4])

k=0.05
RSL5=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL5[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL5[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020+Lij_long)
}
lines(RSL5[1,],RSL5[2,],col=colors[5])

k=0.01
RSL1=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL1[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL1[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020+Lij_long)
}
lines(RSL1[1,],RSL1[2,],col=colors[1])

legend("topright",c("0.01","0.02","0.03","0.04","0.05","0.1"),col=colors,lty=c(1,1,1,1,1,1,1),title="k")

plot(RSL1[1,],RSL2[2,]/RSL2s[2,],type="l",col=colors[2],xlab="R",ylab="relative systemic loss ratio") #total loss/short loss
lines(RSL1[1,],RSL1[2,]/RSL1s[2,],col=colors[1])
lines(RSL1[1,],RSL6[2,]/RSL6s[2,],col=colors[6])
lines(RSL1[1,],RSL3[2,]/RSL3s[2,],col=colors[3])
lines(RSL1[1,],RSL4[2,]/RSL4s[2,],col=colors[4])
lines(RSL1[1,],RSL5[2,]/RSL5s[2,],col=colors[5])
legend("bottomright",c("0.01","0.02","0.03","0.04","0.05","0.1"),col=colors,lty=c(1,1,1,1,1,1,1),title="k")

#Number of defaulted banks as a function of R for different choices of k using β = R, (a, b) = (1, 1) for the model.
k=0.1
a=1
b=1
x=0.06
colors=c("#db6968","#99cbeb","#459943","#e8c559","#606f8a","#f8984e")

Ae2020=A2020-colSums(Lij2020)
NDB6s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB6s[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB6s[2,i]=length(which(E2020[[2]]<=0))
}
plot(NDB6s[1,],NDB6s[2,],type="l",col=colors[6],xlab="R",ylab="number of defaulted banks")

k=0.01
NDB1s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB1s[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB1s[2,i]=length(which(E2020[[2]]<=0))
}
lines(NDB1s[1,],NDB1s[2,],col=colors[1])

k=0.02
NDB2s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB2s[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB2s[2,i]=length(which(E2020[[2]]<=0))
}
lines(NDB2s[1,],NDB2s[2,],col=colors[2])

k=0.03
NDB3s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB3s[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB3s[2,i]=length(which(E2020[[2]]<=0))
}
lines(NDB3s[1,],NDB3s[2,],col=colors[3])

k=0.04
NDB4s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB4s[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB4s[2,i]=length(which(E2020[[2]]<=0))
}
lines(NDB4s[1,],NDB4s[2,],col=colors[4])

k=0.05
NDB5s=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB5s[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB5s[2,i]=length(which(E2020[[2]]<=0))
}
lines(NDB5s[1,],NDB5s[2,],col=colors[5])

legend("topright",c("0.01","0.02","0.03","0.04","0.05","0.1"),col=colors,lty=c(1,1,1,1,1,1,1),title="k")



k=0.1
a=1
b=1
x=0.06
colors=c("#db6968","#99cbeb","#459943","#e8c559","#606f8a","#f8984e")

Ae2020=A2020-colSums(Lij2020+Lij_long)
NDB6=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB6[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB6[2,i]=length(which(E2020[[2]]<=0))
}
plot(NDB6[1,],NDB6[2,],type="l",col=colors[6],xlab="R",ylab="number of defaulted banks")

k=0.01
NDB1=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB1[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB1[2,i]=length(which(E2020[[2]]<=0))
}
lines(NDB1[1,],NDB1[2,],col=colors[1])

k=0.02
NDB2=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB2[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB2[2,i]=length(which(E2020[[2]]<=0))
}
lines(NDB2[1,],NDB2[2,],col=colors[2])

k=0.03
NDB3=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB3[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB3[2,i]=length(which(E2020[[2]]<=0))
}
lines(NDB3[1,],NDB3[2,],col=colors[3])

k=0.04
NDB4=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB4[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB4[2,i]=length(which(E2020[[2]]<=0))
}
lines(NDB4[1,],NDB4[2,],col=colors[4])

k=0.05
NDB5=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB5[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB5[2,i]=length(which(E2020[[2]]<=0))
}
lines(NDB5[1,],NDB5[2,],col=colors[5])

legend("topright",c("0.01","0.02","0.03","0.04","0.05","0.1"),col=colors,lty=c(1,1,1,1,1,1,1),title="k")

plot(NDB1[1,],NDB2[2,]-NDB2s[2,],type="l",col=colors[2],xlab="R",ylab="difference in the number of defaulted banks") #total NDB-short NDB
lines(NDB1[1,],NDB1[2,]-NDB1s[2,],col=colors[1])
lines(NDB1[1,],NDB6[2,]-NDB6s[2,],col=colors[6])
lines(NDB1[1,],NDB3[2,]-NDB3s[2,],col=colors[3])
lines(NDB1[1,],NDB4[2,]-NDB4s[2,],col=colors[4])
lines(NDB1[1,],NDB5[2,]-NDB5s[2,],col=colors[5])
legend("bottomright",c("0.01","0.02","0.03","0.04","0.05","0.1"),col=colors,lty=c(1,1,1,1,1,1,1),title="k")


#######difference between 2019 and 2020
#relative system loss
k=0.1
a=1
b=1
x=0.06
RSL=matrix(0,2,100)
colors=c("#db6968","#99cbeb","#459943","#e8c559","#606f8a","#f8984e")

Ae2019=A2019-colSums(Lij2019+Lij_long)
RSL6.19=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL6.19[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL6.19[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019+Lij_long)
}

k=0.02
RSL2.19=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL2.19[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL2.19[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019+Lij_long)
}

k=0.03
RSL3.19=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL3.19[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL3.19[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019+Lij_long)
}

k=0.04
RSL4.19=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL4.19[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL4.19[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019+Lij_long)
}

k=0.05
RSL5.19=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL5.19[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL5.19[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019+Lij_long)
}

k=0.01
RSL1.19=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL1.19[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL1.19[2,i]=(sum(E2019[[1]])-sum(E2019[[2]]))/sum(Lij2019+Lij_long)
}

k=0.1
a=1
b=1
x=0.06
colors=c("#db6968","#99cbeb","#459943","#e8c559","#606f8a","#f8984e")

Ae2020=A2020-colSums(Lij2020+Lij_long)
RSL6.20=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL6.20[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL6.20[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020+Lij_long)
}

k=0.02
RSL2.20=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL2.20[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL2.20[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020+Lij_long)
}

k=0.03
RSL3.20=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL3.20[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL3.20[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020+Lij_long)
}

k=0.04
RSL4.20=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL4.20[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL4.20[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020+Lij_long)
}

k=0.05
RSL5.20=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL5.20[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL5.20[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020+Lij_long)
}

k=0.01
RSL1.20=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  RSL1.20[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  RSL1.20[2,i]=(sum(E2020[[1]])-sum(E2020[[2]]))/sum(Lij2020+Lij_long)
}

RSL.yc=RSL2.20[2,]/RSL2.19[2,]
yc=which(RSL.yc>3)
for(i in 1:length(yc)){
  RSL.yc[yc[i]]=(RSL.yc[yc[i]-2]+RSL.yc[yc[i]-1])/2
}


plot(RSL1.19[1,],RSL.yc,type="l",col=colors[2],xlab="R",ylab="relative systemic loss ratio") #2020 loss/2019 loss
lines(RSL1.19[1,],RSL1.20[2,]/RSL1.19[2,],col=colors[1])
lines(RSL1.19[1,],RSL6.20[2,]/RSL6.19[2,],col=colors[6])
lines(RSL1.19[1,],RSL3.20[2,]/RSL3.19[2,],col=colors[3])
lines(RSL1.19[1,],RSL4.20[2,]/RSL4.19[2,],col=colors[4])
lines(RSL1.19[1,],RSL5.20[2,]/RSL5.19[2,],col=colors[5])
legend("topright",c("0.01","0.02","0.03","0.04","0.05","0.1"),col=colors,lty=c(1,1,1,1,1,1,1),title="k")

#number of defulted banks
k=0.1
a=1
b=1
x=0.06
colors=c("#db6968","#99cbeb","#459943","#e8c559","#606f8a","#f8984e")

Ae2019=A2019-colSums(Lij2019+Lij_long)
NDB6.19=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB6.19[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB6.19[2,i]=length(which(E2019[[2]]<=0))
}

k=0.01
NDB1.19=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB1.19[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB1.19[2,i]=length(which(E2019[[2]]<=0))
}

k=0.02
NDB2.19=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB2.19[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB2.19[2,i]=length(which(E2019[[2]]<=0))
}

k=0.03
NDB3.19=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB3.19[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB3.19[2,i]=length(which(E2019[[2]]<=0))
}

k=0.04
NDB4.19=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB4.19[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB4.19[2,i]=length(which(E2019[[2]]<=0))
}

k=0.05
NDB5.19=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB5.19[1,i]=R
  beta=R
  E2019<-get_e(n, Epre2019, L2019, Lij2019+Lij_long, Ae2019, x*Ae2019, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB5.19[2,i]=length(which(E2019[[2]]<=0))
}

k=0.1
a=1
b=1
x=0.06
colors=c("#db6968","#99cbeb","#459943","#e8c559","#606f8a","#f8984e")

Ae2020=A2020-colSums(Lij2020+Lij_long)
NDB6.20=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB6.20[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB6.20[2,i]=length(which(E2020[[2]]<=0))
}

k=0.01
NDB1.20=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB1.20[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB1.20[2,i]=length(which(E2020[[2]]<=0))
}

k=0.02
NDB2.20=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB2.20[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB2.20[2,i]=length(which(E2020[[2]]<=0))
}

k=0.03
NDB3.20=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB3.20[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB3.20[2,i]=length(which(E2020[[2]]<=0))
}

k=0.04
NDB4.20=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB4.20[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB4.20[2,i]=length(which(E2020[[2]]<=0))
}

k=0.05
NDB5.20=matrix(0,2,100)
for(i in 1:100){
  R=0.01*i
  NDB5.20[1,i]=R
  beta=R
  E2020<-get_e(n, Epre2020, L2020, Lij2020+Lij_long, Ae2020, x*Ae2020, k, a, b, R, beta, maxitr=1000, delta=1e-5)
  NDB5.20[2,i]=length(which(E2020[[2]]<=0))
}

plot(NDB1.19[1,],NDB2.20[2,]-NDB2.19[2,],type="l",col=colors[2],xlab="R",ylab="difference in the number of defaulted banks") #2020 NDB-2019 NDB
lines(NDB1.19[1,],NDB1.20[2,]-NDB1.19[2,],col=colors[1])
lines(NDB1.19[1,],NDB6.20[2,]-NDB6.19[2,],col=colors[6])
lines(NDB1.19[1,],NDB3.20[2,]-NDB3.19[2,],col=colors[3])
lines(NDB1.19[1,],NDB4.20[2,]-NDB4.19[2,],col=colors[4])
lines(NDB1.19[1,],NDB5.20[2,]-NDB5.19[2,],col=colors[5])
legend("bottomright",c("0.01","0.02","0.03","0.04","0.05","0.1"),col=colors,lty=c(1,1,1,1,1,1,1),title="k")




