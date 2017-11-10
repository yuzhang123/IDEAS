plotHeat<-function(fpref, suffix, fout=NULL, inputlab = NULL, method="complete", orient=0)
{
    if(length(fout)>0)
    {   pdf(file=paste(fout,".pdf",sep=""), height=8-4*orient, width=8+10*orient, onefile=TRUE, family='Helvetica', paper='a4r', pointsize=12);
    }
    source("../source/plot.R");
    rt=readPara1(fpref, suffix=suffix, 1, 14);
    m=t(array(rt$m-0*0.7445,dim=c(14,length(rt$p))));
    rownames(m)=1:dim(m)[1]-1;  
#m=m[rt$p>1e-6,];

    rt1=readPara1("chrHMM.complete", 1, 14, reorder=F);
    m1=t(array(rt1$m,dim=c(14,length(rt1$p))));
    rownames(m1)=paste("[",read.table("chrHMM.labmap")[,2],"]");

    rt2=readPara1("encodeComplete.All", suffix="", 1,14);#, reorder=T);
    m2=t(array(rt2$m, dim=c(14,length(rt2$p))));
    #m2=m2[1:25,];
    rownames(m2)=paste("(", read.table("mineComplete.labmap")[,2], ")");

print(dim(m));print(dim(m1));print(dim(m2));
    data=rbind(m, m1, m2);
    cc = c("red","white","blue")[c(rep(1,dim(m)[1]),rep(2,dim(m1)[1]),rep(3,dim(m2)[1]))];
    if(length(inputlab) > 0) 
    {   data=rbind(m, m2); 
        cc = #rep("red",dim(m)[1]);
            c("red","black")[c(rep(1,dim(m)[1]),rep(2,dim(m2)[1]))]
        rownames(data)[1:dim(m)[1]] = as.character(inputlab);
    } 
    colnames(data)=read.table("encode6.labs")[,1];

    rang=range(data);
    colors=0:100/100*(rang[2]-rang[1])+rang[1];
    my_palette=colorRampPalette(c("dark blue","cyan","green","yellow","red","white"))(n=100);
    source("heatmap3.R");
    myhclust<-function(x){return(hclust(x,method=method))};
    if(orient==0)
    {   heatmap.3(data, Colv=NA, RowSideColors=cc, cellnote=round(data*10)/10, notecol="black", hclust=myhclust, col=my_palette, breaks=colors);
    }
    else
    { heatmap.3(t(data), Rowv=NA, ColSideColors=cc, cellnote=t(round(data*10)/10), notecol="black", hclust=myhclust);#,colsep=1:dim(data)[1],sepwidth=c(0.001,0), col=my_palette, breaks=colors);
    }

    if(length(fout)>0) dev.off();
    return(data);
}

plotP300Signal<-function(mstates, cstates=24, sstates=10:12, showcol = 6, log10=F, datasource=paste0(rep(paste0("P300/",c("Gm12878","H1hesc","Hepg2"),".p300.rep"),each=2),c(1,2)), pref="enhancer_p300_list", remove0 = -1)
{
    X = NULL;
    Lm = Lc = Ls = NULL;
    for(i in 1:length(datasource))
    {   x = read.table(datasource[i]);
        chr=NULL;
        if(length(grep(":",x[1,1]))!=0)
        {   t=t(array(unlist(strsplit(as.character(x[,1]),":")),dim=c(2,dim(x)[1])));
            t2=t(array(unlist(strsplit(t[,2],"-")),dim=c(2,dim(x)[1])));
            chr = as.numeric(substring(t[,1],4));
            chr[which(t[,1]=="chrX")] = 23;
            chr[which(t[,1]=="chrY")] = 24;
            x = cbind(1:dim(x)[1],chr,as.numeric(t2[,1]),as.numeric(t2[,2]),as.numeric(x[,2]));
            nx = NULL;
            for(j in 1:24)
            {   t=which(x[,2]==j);
                if(length(t)>0)
                {   nx=rbind(nx,x[t[order(x[t,3])],]);
                }
            }
            x=nx;
        }
        else
        {   chr = as.integer(substr(x[,2],4,100));
            chr[which(x[,2]=="chrX")] = 23;
            chr[which(x[,2]=="chrY")] = 24;
        }
        nx = NULL;
        for(k in 1:24)
        {   t = which(chr == k);
            if(length(t) == 0) next;
            nx = rbind(nx, x[t[order(x[t,3])],]);
        }
        x = nx;
        X = rbind(X, x);
print(i);
print(dim(x)[1]);
        y = read.table(paste(pref,"m",i-1, sep="."));
print(dim(y)[1]);
        sel = rep(0, dim(y)[1]);
        for(k in mstates) sel[which(y[,2+k] > 0)] = 1;
        Lm = c(Lm, sel);
        #Lm = c(Lm, as.integer(apply(as.matrix(y[,2+mstates]),1,sum)>0));
        sm = apply(y[,-1], 1, sum);

        y = read.table(paste(pref,"c",i-1, sep="."));
print(dim(y)[1]);
        ss = cstates;
        sel = rep(0, dim(y)[1]);
        for(k in ss) sel[which(y[,2+k] > 0)] = 1;
        Lc = c(Lc, sel);
        #Lc = c(Lc, as.integer((y[,2+24]+y[,2+3]>0)));
        #Lc = c(Lc, as.integer(apply(as.matrix(y[,2+ss]),1,sum)>0));

        y = read.table(paste(pref,"s",i-1, sep="."));
print(dim(y)[1]);
        ss = sstates;
        sel = rep(0, dim(y)[1]);
        for(k in ss) sel[which(y[,2+k] > 0)] = 1;
        Ls = c(Ls, sel);
        #Ls = c(Ls, as.integer((y[,2+10]+y[,2+11]+y[,2+12]+y[,21]+y[,22]+y[,23]>0)));
        #Ls = c(Ls, as.integer(apply(as.matrix(y[,2+ss]),1,sum)>0));
    }

    par(mfrow=c(1,2));

    id = showcol;
    print(id);print(dim(X));

    if(remove0>=0)
    {   t=which(X[,id]<remove0);
        X=X[-t,];
        Lm=Lm[-t];Lc=Lc[-t];Ls=Ls[-t];
    }

    t11 = which(Lm > 0 & Lc > 0);
    t10 = which(Lm > 0 & Lc == 0);
    t01 = which(Lm == 0 & Lc > 0);
    
    if(log10==T) X[,id]=log10(X[,id]+min(X[X[,id]>0,id]/1));

    boxplot(X[t11,id], X[t10,id], X[t01,id], col=c("white", "red", "blue"), labs=c("Common", "IDEAS Only", "ChrHMM Only"), varwidth=T);
    print(c(length(t11),length(t10),length(t01)));

    t11 = which(Lm > 0 & Ls > 0);
    t10 = which(Lm > 0 & Ls == 0);
    t01 = which(Lm == 0 & Ls > 0);
    boxplot(X[t11,id], X[t10,id], X[t01,id], col=c("white", "red", "blue"), labs=c("Common", "IDEAS Only", "Segway Only"), varwidth=T);
    print(c(length(t11),length(t10),length(t01)));

    return(cbind(X,Lm, Lc, Ls));
}

plotPair<-function(fpair, labs = NULL, enrich = T, celltext = F, cut = 100, rang = NULL)
{   x = read.table(fpair);
    x = as.matrix(x);
 
    mm=apply(x,1,sum);
        nn=apply(x,2,sum);
    x=x[mm>=cut,nn>=cut];
    mm=mm[mm>=cut];
    nn=nn[nn>=cut];
        ee=mm%*%t(nn)/sum(x);
        if(enrich==TRUE)
        {
                rr=log(abs(x-ee)/(sqrt(ee+1e-10)*sqrt(sum(x)*0+1))+1)*sign(x-ee);
        }
        else
        {       rr=log((x+cut)/(ee+cut));
        }

    source("heatmap3.R");
        crmstates=labs;
        if(length(labs)==0)
        { crmstates = 1:dim(rr)[1]-1; }
print(dim(rr));
print(crmstates);
        rownames(rr)=colnames(rr)=crmstates;

    if(length(rang) == 0)
    {
        rang=range(rr);
        rang[2]=max(rang[2],abs(rang[1]));
        rang[1]=-max(abs(rang[1]),rang[2]);
    }
    print(rang);
    colors=0:100/100*(rang[2]-rang[1])+rang[1];
    my_palette=colorRampPalette(c("white","blue","black","red","yellow"))(n=100);

    if(celltext == T)
    {   cx = round(rr * 1)/1;
        o=heatmap.3(rr, cellnote = cx, notecol = array(c("white","black")[2+(sign(rr)-1)/2],dim=dim(rr)), breaks=colors, col=my_palette);
    }
    else
    {   heatmap.3(rr, breaks=colors, col=my_palette);
    }
    return(rr);
}

plotP300Pie<-function(fpref = "enhancer_p300")
{   names=c("mineCompleten4","chrHMM","segway");
    suffix=c("mineCompleten4","chrHMMf","segwayf");
    
    par(mfrow=c(1,3));
    for(i in 1:3)
    {   x=read.table(paste0(fpref, ".",suffix[i]));
        if(length(which(is.na(x[,6])==T))>0)
        {   x[which(is.na(x[,6])==T),6]=0;
        }
        lab=read.table(paste0(names[i],".labmap"));

        if(i==1) 
        { mycol=c("lightyellow", "lightblue", "mistyrose", "lightcyan", "lavender", "cornsilk"); }
        else if(i==2) 
        { mycol=c("mistyrose", "lightyellow", "lightcyan", "lightblue", "lavender", "cornsilk"); }
        else if(i==3) 
        { mycol=c("lightyellow", "mistyrose", "lightblue", "lightcyan", "lavender", "cornsilk"); }
        pie(x[,6],labels=paste0(lab[x[which(x[,6]>2),1]+1,2]," ",as.integer(x[x[,6]>2,6]),"%"),clockwise=T, main=names[i], col=mycol, init.angle=0);
    }
}

plotStatePairEnrich<-function(rt,method="complete", celltext=T)
{   rt=as.matrix(rt);
    rang=range(rt);
    colors=0:100/100*(rang[2]-rang[1])+rang[1];
    my_palette=colorRampPalette(c("white","blue","black","red","yellow"))(n=100);
    
    source("heatmap3.R");
    myhclust<-function(x){return(hclust(x,method=method))};
    if(celltext==T)
    {   cx = round(rt);
        heatmap.3(rt,symm=T,hclustfun = myhclust, breaks=colors,col=my_palette, cellnote = cx, notecol = (array(c("white","black")[2+(sign(rt)-1)/2],dim=dim(rt))), drawlayout=T);
    }
    else
    {   heatmap.3(rt,symm=T,hclustfun = myhclust, breaks=colors,col=my_palette, drawlayout=T);
    }
}

plotHeterogeneity<-function(fpair, enrich=T, cut=100)
{   
    par(mfrow=c(1,3));
    x=read.table(fpair);
    x=as.matrix(x);
    t=which(apply(x,1,sum)>=cut);
    x=x[t,t];
    m=apply(x,1,sum);
    p=median(diag(x)/m);
p=sum(diag(x))/sum(m);
    h=log((m-diag(x) + cut)/(m*(1-p)+cut));
if(enrich) { e=m*(1-p);h=sign(m-diag(x)-e)*log(abs(m-diag(x)-e)/(sqrt(e)+1e-10)); }
    r=range(h*1.2);
    barplot(h[order(h)],names.arg=(1:dim(x)[1]-1)[order(h)], ylim=r, col=rainbow(dim(x)[1]));

    x=read.table("chrHMM.pair");
    x=as.matrix(x);
    t=which(apply(x,1,sum)>=cut);
    x=x[t,t];
    m=apply(x,1,sum);
    p=median(diag(x)/m);
p=sum(diag(x))/sum(m);
    h=log((m-diag(x) + cut)/(m*(1-p)+cut));
if(enrich) { e=m*(1-p);h=sign(m-diag(x)-e)*log(abs(m-diag(x)-e)/(sqrt(e)+1e-10)); }
    barplot(h[order(h)],names.arg=(1:dim(x)[1]-1)[order(h)],ylim=r, col=rainbow(dim(x)[1]));

    x=read.table("segway.pair");
    x=as.matrix(x);
    t=which(apply(x,1,sum)>=cut);
    x=x[t,t];
    m=apply(x,1,sum);
    p=median(diag(x)/m);
p=sum(diag(x))/sum(m);
    h=log((m-diag(x) + cut)/(m*(1-p)+cut));
if(enrich) { e=m*(1-p);h=sign(m-diag(x)-e)*log(abs(m-diag(x)-e)/(sqrt(e)+1e-10)); }
    barplot(h[order(h)],names.arg=(1:dim(x)[1]-1)[order(h)],ylim=r, col=rainbow(dim(x)[1]));
}

plotPrecisionRecall<-function(pref, id = 1:6, i1 = 7, i2 = 9, xmax=-1, ymax=-1)
{   k = 2;  
    for(i in c(id, 100, 101))
    {   
        if(i < 100) { x = read.table(paste(pref, ".mineComplete",i,sep="")); }
        else if(i == 100) { x = read.table(paste(pref, ".chrHMMf",sep="")); }
        else if(i == 101) { x = read.table(paste(pref, ".segwayf",sep="")); }
        
        if(ymax < 0)
        {   ymax = 100;
            if(i2 == 9) { ymax = max(x[,i2]) * 1.5; }
        }
        if(xmax < 0) xmax=min(max(x[,i1]) + 10,100);
        if(i == id[1]) { plot(x[,i1], x[,i2],xlim=c(0,xmax),ylim=c(0,ymax),type="o", col=k); }
        else { lines(x[,i1],x[,i2],type="o",col=k,pch=k); }
        k=k+1;
    }
}

plotPrecisionRecallFixState<-function(pref, id = "mineCompleten4", i1 = 2, i2 = 4, xmax=-1, ymax=-1, smooth=F, cex=2, mode=0)
{   k = 2;  
    id = c(id, "chrHMMf", "segwayf", "enhfinder");
    n = 0;
    X = NULL;
    for(i in id)
    {   x = read.table(paste0(pref, ".fixstate.", i));
        n = dim(x)[1];
        if(smooth==T & dim(x)[1]>2)
        {   x[,i2]=predict(loess(x[,i2]~x[,i1]),x[,i1]);    
        }
        X = rbind(X, x);
        
        if(ymax < 0)
        {   ymax = min(100, max(x[which(is.na(x[,i2])==F),i2]) + 20);
            if(i2 == 9) { ymax = max(x[which(is.na(x[,i2])==F),i2]) * 1.5; }
        }
        if(xmax < 0) xmax=min(max(x[,i1]) + 10, 100);
        l=1:dim(x)[1];
        if(length(which(is.na(x[,i2])==F))>10) l=seq(1,length(l),2);
        if(i == id[1]) { plot(x[l,i1], x[l,i2],xlim=c(0,xmax),ylim=c(0,ymax),type="o", col=k, pch=k-1, cex=cex); }
        else { lines(x[l,i1],x[l,i2],type="o",col=k,pch=k-1,cex=cex); }
        k=k+1;
        if(mode==1 & i == "segwayf") break;
        if(i == "segwayf") k=k+2;
    }
}

plotFixStatePR_Enh<-function(i1=2, i2=4,smooth=T, xmax = -1, ymax = -1, cex=2, mode=0)
{   

if(FALSE)
{
    x=read.table("enhancer_vista.mineCompleten4");
    y=read.table("enhancer_vista.chrHMMf");
    z=read.table("enhancer_vista.10.segwayf");
u=read.table("enhancer_vista.enhfinder");
    
    par(mfrow=c(1,3+mode));
    par(mar=c(2,2,2,2));
    i=which(x[,1]==19);
    tx=xmax;ty=ymax;
    if(tx < 0) tx = min(100, x[i,i1] + 20);
    if(ty < 0) ty = min(u[u[,1]==0,i2] + 20, 100);
    k=2;
    plot(x[i,i1], x[i,i2], xlim=c(0,tx), ylim=c(0,ty), type="o", col=k,pch=k-1, cex=cex);
    k=k+1;
    i=which(y[,1]==24);
    points(y[i,i1], y[i,i2], col=k,pch=k-1, cex=cex);
    k=k+3;
    i=which(z[,1]==10);
    points(z[i,i1],z[i,i2],col=k,pch=k-1, cex=cex); 
k=k+1;
i=which(u[,1]==0);
points(u[i,i1],u[i,i2],col=k,pch=k-1,cex=cex);
        
    plotPrecisionRecallFixState("enhancer_cagetest",i1=i1,i2=i2,smooth=smooth,xmax=xmax,ymax=ymax,cex=cex);
    plotPrecisionRecallFixState("enhancer_cageusage",i1=i1,i2=i2,smooth=smooth,xmax=xmax,ymax=ymax,cex=cex);
    if(mode==1) plotPrecisionRecallFixState("enhancer_efinder",i1=i1,i2=i2,smooth=smooth,xmax=xmax,ymax=ymax,cex=cex);
}

    par(mfrow=c(1,4));
    par(mar=c(2,2,2,2));
    plotPrecisionRecallFixState("enhancer_vista",i1=i1,i2=i2,smooth=smooth,xmax=xmax,ymax=ymax,cex=cex,mode=1);
    plotPrecisionRecallFixState("enhancer_cagetest",i1=i1,i2=i2,smooth=smooth,xmax=xmax,ymax=ymax,cex=cex,mode=0);
    plotPrecisionRecallFixState("enhancer_cageusage",i1=i1,i2=i2,smooth=smooth,xmax=xmax,ymax=ymax,cex=cex,mode=0);
    plotPrecisionRecallFixState("enhancer_efinder",i1=i1,i2=i2,smooth=smooth,xmax=xmax,ymax=ymax,cex=cex,mode=1);
}


plotTssTdd<-function(pref, id = 10, rm=NULL, type = 1, smooth=0, pc1=1,pc2=2, K=5,method="ward")
{   x=read.table(paste(pref, ".",id, sep=""));
    y=read.table(paste(pref, ".chrHMMf", sep=""));
    z=read.table(paste(pref, ".segwayf", sep=""));
    x=as.matrix(x);
    y=as.matrix(y);
    z=as.matrix(z);

    mode = 1;
    if(substr(pref,1,3)=="tdd")
    {   mode = 2;
    }
        
    if(mode == 1)
    {   if(length(rm)==0) { rm = which(x[-1,1]<100); if(length(rm)>0){rm=rm+1;}}
        if(length(rm)>0) x=x[-rm,];
        mx=apply(x[-1,-1]*x[-1,1],2,sum);
        x = t(x[-1, -1]);
        my=apply(y[-1,-1]*y[-1,1],2,sum);
        y = t(y[-1, -1]);
        mz=apply(z[-1,-1]*z[-1,1],2,sum);
        z = t(z[-1, -1]);
    }
    else
    {   if(length(rm)==0) { rm = which(apply(x,2,sum)<100); }
        if(length(rm)>0) { x=x[,-rm]; }
        mx=x[,dim(x)[2]];
        x = x[,-dim(x)[2]];
        my=y[,dim(y)[2]];
        y = y[,-dim(y)[2]];
        mz=z[,dim(z)[2]];
        z = z[,-dim(z)[2]];
    }

    nnn=dim(x)[1];
    if(smooth>0)
    {   for(i in 1:dim(x)[2]) x[,i]=density(1:nnn,weights=(x[,i] + 1)/sum(x[,i] + 1),bw=smooth,n=nnn,from=1,to=nnn)$y; 
        for(i in 1:dim(y)[2]) y[,i]=density(1:nnn,weights=(y[,i] + 1)/sum(y[,i] + 1),bw=smooth,n=nnn,from=1,to=nnn)$y; 
        for(i in 1:dim(z)[2]) z[,i]=density(1:nnn,weights=(z[,i] + 1)/sum(z[,i] + 1),bw=smooth,n=nnn,from=1,to=nnn)$y; 
        mx=density(1:nnn,weights=(mx + 1)/sum(mx + 1),bw=smooth,n=nnn,from=1,to=nnn)$y; 
        my=density(1:nnn,weights=(my + 1)/sum(my + 1),bw=smooth,n=nnn,from=1,to=nnn)$y; 
        mz=density(1:nnn,weights=(mz + 1)/sum(mz + 1),bw=smooth,n=nnn,from=1,to=nnn)$y; 
    }

    n = apply(x, 2, sum);
    x = (t(t(x) / n) + 1e-8) / (mx / sum(mx) + 1e-8);
    n = apply(y, 2, sum);
    y = (t(t(y) / n) + 1e-8) / (my / sum(my) + 1e-8);
    n = apply(z, 2, sum);
    z = (t(t(z) / n) + 1e-8) / (mz / sum(mz) + 1e-8);

    if(file.exists(paste(id,".labmap",sep=""))==TRUE) colnames(x)=read.table(paste(id,".labmap",sep=""))[-rm,2];
    colnames(y)=read.table("chrHMM.labmap")[,2];#paste("C",1:dim(y)[2]-1,sep="");
    colnames(z)=read.table("segway.labmap")[,2];#paste("S",1:dim(z)[2]-1,sep="");
    cc = c("red","green","blue")[c(rep(1,dim(x)[2]),rep(2,dim(y)[2]),rep(3,dim(z)[2]))];
    source("heatmap3.R");

    rr = log(cbind(x,y,z));
    rang=range(rr);
#   rang[2]=max(rang[2],abs(rang[1]));
#   rang[1]=-max(abs(rang[1]),rang[2]);
    
    if(type == 0)
    {   x=log(x);
        r=range(x);
        m=cutree(hclust(dist(t(x)),method=method),k=K);
        print(m);
        par(mfrow=c(2,round(K/2+0.1)));
        par(mar=c(2,2,2,2));
        labs=colnames(x);
        for(i in 1:K)
        {   t=which(m==i);
            plot(x[,t[1]], col=i, type="l", ylim=r);    
            abline(0,0,lty=2);
            lines(rep(dim(x)[1]/4+0.5,2),c(-100,100),lty=2);
            lines(rep(dim(x)[1]-dim(x)[1]/4-0.5,2),c(-100,100),lty=2);
            if(length(t)>1)
            {   for(j in 2:length(t)) lines(x[,t[j]],col=i);
            }
            for(j in t)
            {   tx=which(x[,j]==max(x[,j]))[1];
                ty=x[tx,j];
                text(tx,ty,labels=labs[j]);
            }
        }   
    }
    else if(type == 1)
    {   colors=0:100/100*(rang[2]-rang[1])+rang[1];
        my_palette=colorRampPalette(c("white","blue","black","red","yellow"))(n=100);
        heatmap.3(rr,breaks=colors, col = my_palette, Rowv=NA,ColSideCol=cc);
    }
    else
    {   
        tt=c(rep(2,dim(x)[2]),rep(3,dim(y)[2]),rep(4,dim(z)[2]));
        rr=plotPCA(rr, apply(rr,1,mean), 0, tt, pc1=pc1,pc2=pc2);
    }

    return(rr);
}

plotPCA<-function(rr, mm, pcn = 0, col = NULL, pc1=1, pc2=2)
{   rr0 = t(rr); 
rr0=rr0-apply(rr0,1,mean);
    mm0=apply(rr0,2,mean);
mm=c(mm);
mm=mm0;
mm=rep(mean(rr0),dim(rr0)[2]);
#plot(mm0,mm);abline(0,1);return(cbind(mm0,mm));

    rr0 = t(t(rr0) - mm);
    r=svd(rr0);
orr = r;
r=prcomp(rr0);
pp = summary(r)[[6]];
if(pcn==0) pcn = min(which(pp[3,]>=0.9));

pcn=min(pcn,5);

print(range(r$rotation));
print(summary(r));
r$d=r$sdev;
r$v=r$rotation;

r=orr;

    z=rr0%*%r$v;
    lmat=cbind(c(1:pcn),pcn+1,pcn+1,pcn+1);
        layout(lmat);
    par(mar=c(2,2,2,2));
    for(i in 1:pcn) 
    {   f=1;    
        plot(mm,type="l",lty=3,ylim=range(cbind(mm+r$v[,i]*f,mm-r$v[,i]*f)),main=paste("PC",i," (",round(r$d[i]*10)/10,",",round(pp[2,i]*100)/100,")",sep=""));
        lines(mm+r$v[,i]*f,lty=1);
        #lines(mm-r$v[,i]*f,lty=3);

        lines(rep(dim(rr)[1]/4+0.5,2),c(-1000,1000),lty=3);
        lines(rep(dim(rr)[1]*3/4+0.5,2),c(-1000,1000),lty=3);
        abline(0,0,lty=3);
    }
    tt=col;
    if(length(tt)==0) tt = rep(1,dim(z)[1]);
    #plot(z[,pc1],z[,pc2],col=tt,pch=tt);
    plot(0,0,xlim=range(z[,c(pc1,pc2)]),ylim=range(z[,c(pc1,pc2)]),col="white",xlab=paste("PC",pc1,sep=""),ylab=paste("PC",pc2,sep=""));
    abline(0,0, lty=3);lines(c(0,0),c(-100000,100000), lty=3);
    
    t=table(tt);
    t=as.matrix(t);
    n=as.numeric(rownames(t));
    o=order(t);
    t=t[o];
    n=n[o];
    for(i in length(t):1)
    {   f=which(tt==n[i]);
        text(z[f,pc1],z[f,pc2],rownames(rr0)[f],col=tt[f]);
    }

    print(r$d);

    return(z);
}

plotPairTssTdd<-function(showsuf="mineComplete2",rm=NULL, smooth=0, distrib=FALSE, pc1=1, pc2=2, fout=NULL, a=100)
{   
    if(length(fout)>0)
    {   pdf(file=paste(fout,".pdf",sep=""), height=8, width=8, onefile=TRUE, family='Helvetica', paper='a4r', pointsize=12);
    }
    N=tt=M=NULL;
    for(i in 1:length(showsuf))
    {   n=getPairTssTddMatrix(showsuf[i], NULL, smooth, distrib, a);
        nm=n[,dim(n)[2]]; n=n[,-dim(n)[2]];
        N=cbind(N,n);
        tt=c(tt,rep(i+1,dim(n)[2]));
        M=cbind(M,nm);
    }

    if(length(showsuf)==1)
    {   cl=colnames(N);
        t=regexpr(":",cl)[1:dim(N)[2]];
        ta=substring(cl,1,t-1);
        tb=substring(cl,t+1);
        tt[which(ta==tb)]=1;
    }
    else 
    {   M=apply(M,1,mean);
    }
    n=plotPCA(N,M,col=tt,pc1=pc1,pc2=pc2);
    
    if(length(fout)>0) dev.off();

    rt=NULL;
    rt$N=N;
    rt$z=n;
    rt$h=M;
    return(rt);
}

getPairTssTddMatrix<-function(suffix,rm,smooth,distrib, a=100)
{   x=read.table(paste("tddpair.gencode7.", suffix, sep=""));
    x=as.matrix(x);
    sx=x[,dim(x)[2]];
    x=x[,-dim(x)[2]];
    nnn=dim(x)[1];
    mmm=apply(x,1,sum);
    x=x/mmm;

    y=read.table(paste("tdd.gencode7.", suffix, sep=""));
    y=as.matrix(y);
    y=y[,-dim(y)[2]];
    sy=apply(y,2,sum);
    mmmy=apply(y,1,sum);
    y=y/mmmy;
    
    oa=a;a=rep(a,dim(x)[1])/mmm;
    if(smooth>0) 
    { for(i in 1:dim(x)[2]) x[,i]=sum(x[,i])*density(1:nnn,weights=(x[,i] + a)/sum(x[,i] + a),bw=smooth,n=nnn,from=1,to=nnn)$y;     
      if(FALSE)
      {
        k=1;
        for(i in 1:(dim(y)[2])) 
            for(j in 1:(dim(y)[2]))
            {   if(j > i) next;
                a=rep(oa,dim(x)[1])*(as.integer(i!=j)+1)*sy[i]*sy[j]/sum(sy)^2*dim(x)[2]/mmm;
                x[,k]=sum(x[,k])*density(1:nnn,weights=(x[,k] + a)/sum(x[,k] + a),bw=smooth,n=nnn,from=1,to=nnn)$y; 
                k=k+1;
            }
      }
    }
    if(distrib==FALSE) { x=x/apply(x,1,sum); }
    else { x=x/apply(x,1,sum); if(FALSE){x = t(t(x)/apply(x,2,sum));} 
        sx = sx / sum(sx); }
        
    a=oa;a=rep(a,dim(y)[1])/mmmy;
    if(smooth>0) 
    #{ for(i in 1:dim(y)[2]) y[,i]=sum(y[,i])*density(1:nnn,weights=(y[,i] + a)/sum(y[,i] + a),bw=smooth,n=nnn,from=1,to=nnn)$y; }
     { for(i in 1:dim(y)[2]) y[,i]=sum(y[,i])*density(1:nnn,weights=(y[,i] + a*sy[i]/sum(sy)*dim(y)[2])/sum(y[,i] + a*sy[i]/sum(sy)*dim(y)[2]),bw=smooth,n=nnn,from=1,to=nnn)$y; } 
    y=y/apply(y,1,sum);

    labs=1:dim(y)[2]-1;
    if(file.exists(paste(suffix,".labmap",sep=""))==TRUE) labs=read.table(paste(suffix,".labmap",sep=""))[,2];
    if(suffix=="chrHMMf") labs=read.table("chrHMM.labmap")[,2];
    if(suffix=="segwayf") labs=read.table("segway.labmap")[,2];
    n=lab=NULL;
    states = 1:dim(y)[2];
    if(length(rm)==0)
    {   rm=which(sy < 100);
    }
    if(length(rm)>0) states=states[-rm];
    tt=NULL;
    for(i in states)
    {   for(j in states)
        {   if(j > i) next;
            t=(i-1)*i/2+j;
            if(i==j) tt=c(tt,t);
        #   print(c(i,j,t));
            if(distrib==TRUE) { n=cbind(n,log(x[,t]+1e-10));}
            else 
            { #n=cbind(n,(log(x[,t]+1e-10)-log(y[,i]*y[,j]*(as.integer(i!=j)+1)+1e-10))); 
              n=cbind(n, ((x[,t])-(y[,i]*y[,j]))/sqrt(y[,i]*(1-y[,i])*y[,j]*(1-y[,j])));
            }
            lab=c(lab,paste(labs[i],labs[j],sep=":"));
        }
    }
    
    if(distrib==TRUE)
    {   n = cbind(n,sx);
        lab=c(lab, "mean");
    }   
    else
    {   #print(tt);
        n = cbind(n, log(apply(x[,tt],1,sum)+1e-10) - log(apply(y,1,function(x){sum(x^2)})+1e-10));
        lab=c(lab, "mean");
    }

    colnames(n)=lab;
    return((n));
}

plotHomoEnrich<-function(x)
{   cl=colnames(x);
    t=regexpr(":",cl)[1:dim(x)[2]];
    ta=substring(cl,1,t-1);
    tb=substring(cl,t+1);
    tt=which(ta==tb);   

#   x=t(t(x)-apply(x,2,mean));
    
    source("heatmap3.R");
    rang=range(x[,tt]);
    colors=0:100/100*(rang[2]-rang[1])+rang[1];
    my_palette=colorRampPalette(c("white","blue","black","red","yellow"))(n=100);
    heatmap.3(t(x[,tt]),breaks=colors, col = my_palette, Colv=NA);
}

plotCountLength<-function(gfile)
{   g=read.table(gfile);
    print("gfile read");
    celln = dim(g)[2]-4;
    K=max(g[,3+1:celln]);
    k=K+1;
    n=hist(as.matrix(g[,3+1:celln]),breaks=(-1):K,plot=F)[[2]];
    
    len = NULL;
    l = dim(g)[1];
    t=g[2:l,3]-g[2:l-1,3];
    step = min(t[t>0]);
    for(i in 1:celln)
    {   print(i);
        t=which(g[2:l,2]!=g[2:l-1,2] | g[2:l,3]>g[2:l-1,3]+step | g[2:l,3+i]!=g[2:l-1,3+i]);
        s=c(0,t)+1;
        t=c(t,l);
        len=rbind(len, cbind(g[s,3+i],g[t,3]-g[s,3]+step));
    }
    rt=NULL;
    rt$count=n;
    rt$length=len;
    return(rt);
}

plotStateMethylation<-function()
{
    x=read.table("Methyl/mineCompleten4.state_methylation");
    m=0;for(i in 1:25) m[i]=median(x[x[,1]==i-1,2])
    o=order(m)
    map=rep(0,25)
    map[o]=0:24
    r=boxplot(x[,2]~map[x[,1]+1],names=read.table("mineCompleten4.labmap")[o,2],las=2, outline=F,ylim=range(x[,2]));
    gg=oo=NULL;
    for(i in unique(r$group))
    {   t=which(r$group == i);
        u=unique(r$out[t]);
        oo=c(oo,u);
        gg=c(gg,rep(i,length(u)));
    }
    points(gg,oo);
}   

prepareEFstate<-function(g, cut=0.001)
{   cellnames=c("Gm12878","H1hesc","Helas3","Hepg2","Huvec","K562");
    l=dim(g)[1];
    t=which(g[2:l,2]!=g[2:l-1,2] | g[2:l,3]!=g[2:l-1,3]+200);
    st=c(1,t+1);
    ed=c(t,l);
    chr=as.numeric(substring(g[,2],4));
    chr[which(g[,2]=="chrX")]=23;
    chr[which(g[,2]=="chrY")]=24;
    inv=cbind(chr[st],g[st,3],g[ed,3]+200);
    for(i in 1:6)
    {   print(i);
        x=read.table(paste0("enhancerfinder.",cellnames[i]));
        t=t(array(unlist(strsplit(as.character(x[,1]),":")),dim=c(2,dim(x)[1])));
        t2=t(array(unlist(strsplit(t[,2],"-")),dim=c(2,dim(x)[1])));
        chr = as.numeric(substring(t[,1],4));
        chr[which(t[,1]=="chrX")] = 23;
        chr[which(t[,1]=="chrY")] = 24;
        x = cbind(chr,as.numeric(t2[,1]),as.numeric(t2[,2]),as.numeric(x[,2]));
        x = x[x[,4]>=cut,];
        nx = NULL;
        for(j in 1:24)
        {   t=which(x[,1]==j);
            if(length(t)>0)
            {   nx=rbind(nx,x[t[order(x[t,2])],]);
            }
        }
        x=nx;
    
        newinv = inv;
        
        k=1;l=dim(inv)[1];ll=dim(x)[1];
        stin=rep(-1,ll);
        for(j in 1:ll)
        {   while(k<=l & (inv[k,1]<x[j,1] | (inv[k,1]==x[j,1] & inv[k,3]<=x[j,2]))) { k=k+1; if(k>l) break; }   
            if(k>l) break;
            if(inv[k,1]==x[j,1] & inv[k,2] < x[j,2]) stin[j]=k;
            if(k>=l) break;
        }
        t=which(stin>0);
        if(length(t)>0)
        {   
            k=sort(c(1:dim(inv)[1],stin[t]));
            newinv=inv[k,];
            lll=length(t);
            print(lll);
            for(j in 1:lll)
            {   jj=stin[t[j]]
                m=jj+j-1;
                newinv[m,3]=x[t[j],2];
                newinv[m+1,2]=x[t[j],2];
            }
        }

        myinv=newinv;
        k=1;l=dim(myinv)[1];ll=dim(x)[1];
        edin=rep(-1,ll);
        for(j in 1:ll)
        {   while(k<=l & (myinv[k,1]<x[j,1] | (myinv[k,1]==x[j,1] & myinv[k,3]<=x[j,3]))) { k=k+1; if(k>l) break;}  
            if(k>l) break;
            if(myinv[k,1]==x[j,1] & myinv[k,2] < x[j,3]) edin[j]=k;
            if(k>=l) break;
        }
        t=which(edin>0);
        if(length(t)>0)
        {   
            k=sort(c(1:dim(myinv)[1],edin[t]));
            newinv=myinv[k,];
            lll=length(t);
            print(lll);
            for(j in 1:lll)
            {   jj=edin[t[j]]
                m=jj+j-1;
                newinv[m,3]=x[t[j],3];
                newinv[m+1,2]=x[t[j],3];
            }
        }
        k=1;l=dim(newinv)[1];ll=dim(x)[1];

        states=rep(1,l);
        print(l);
        for(j in 1:l)
        {   while(k<=ll & (x[k,1]<newinv[j,1] | (x[k,1]==newinv[j,1] & x[k,3]<newinv[j,3]))) { k=k+1;}  
            if(x[k,1]==newinv[j,1] & x[k,2] <= newinv[j,2]) states[j]=0;
            if(k>=ll) break;
        }
        print("done");

        options(scipen=10);
        write.table(cbind(newinv, states), paste0(cellnames[i], ".enhfinder"), quote=F,row.names=F,col.names=F);
    }
}

compareSignal<-function(mstates=t(19), cstates=t(24), sstates=t(10:12), estates=t(0), showcol = 6, log10=F, datasource=paste0(rep(paste0("P300/",c("Gm12878","H1hesc","Hepg2"),".p300.rep"),each=2),c(1,2)), pref="enhancer_p300_list", remove0 = -1, orient = -1, sign = 1)
{
    if(is.null(dim(mstates))==T) { mstates = t(mstates); }
    if(is.null(dim(cstates))==T) { cstates = t(cstates); }
    if(is.null(dim(sstates))==T) { sstates = t(sstates); }
    if(is.null(dim(estates))==T) { estates = t(sstates); }
    X = NULL;
    Lm = Lc = Ls = Le = NULL;
    for(i in 1:length(datasource))
    {   
        if(FALSE)
        {
        x = read.table(datasource[i]);
        chr=NULL;
        if(length(grep(":",x[1,1]))!=0)
        {   t=t(array(unlist(strsplit(as.character(x[,1]),":")),dim=c(2,dim(x)[1])));
            t2=t(array(unlist(strsplit(t[,2],"-")),dim=c(2,dim(x)[1])));
            chr = as.numeric(substring(t[,1],4));
            chr[which(t[,1]=="chrX")] = 23;
            chr[which(t[,1]=="chrY")] = 24;
            x = cbind(1:dim(x)[1],chr,as.numeric(t2[,1]),as.numeric(t2[,2]),as.numeric(x[,2]));
            nx = NULL;
            for(j in 1:24)
            {   t=which(x[,2]==j);
                if(length(t)>0)
                {   nx=rbind(nx,x[t[order(x[t,3])],]);
                }
            }
            x=nx;
        }
        else
        {   chr = as.integer(substr(x[,2],4,100));
            chr[which(x[,2]=="chrX")] = 23;
            chr[which(x[,2]=="chrY")] = 24;
        }
        nx = NULL;
        for(k in 1:24)
        {   t = which(chr == k);
            if(length(t) == 0) next;
            nx = rbind(nx, x[t[order(x[t,3])],]);
        }
        x = nx;
        X = rbind(X, x);
print(i);
print(dim(x)[1]);
        }

        y = read.table(paste(pref,"m",i-1, sep="."));
print(dim(y)[1]);
        sel = array(0, dim=c(dim(y)[1],dim(mstates)[1]));
        for(k in 1:dim(mstates)[1]) 
        {   t = unique(mstates[k,]);
            t = t[t>=0];
            for(l in t)
            {   sel[which(y[,3+l] > 0),k] = 1;
            }
        }
        Lm = rbind(Lm, sel);
X=rbind(X,as.matrix(y[,2]));

        y = read.table(paste(pref,"c",i-1, sep="."));
print(dim(y)[1]);
        sel = array(0, dim=c(dim(y)[1],dim(cstates)[1]));
        for(k in 1:dim(cstates)[1]) 
        {   t = unique(cstates[k,]);
            t = t[t>=0];
            for(l in t)
            {   sel[which(y[,3+l] > 0),k] = 1;
            }
        }
        Lc = rbind(Lc, sel);

        y = read.table(paste(pref,"s",i-1, sep="."));
print(dim(y)[1]);
        sel = array(0, dim=c(dim(y)[1],dim(sstates)[1]));
        for(k in 1:dim(sstates)[1]) 
        {   t = unique(sstates[k,]);
            t = t[t>=0];
            for(l in t)
            {   sel[which(y[,3+l] > 0),k] = 1;
            }
        }
        Ls = rbind(Ls, sel);

        if(max(estates) >= 0)
        {   y = read.table(paste(pref,"e",i-1, sep="."));
print(dim(y)[1]);
            sel = array(0, dim=c(dim(y)[1],dim(estates)[1]));
            for(k in 1:dim(estates)[1]) 
            {   t = unique(estates[k,]);
                t = t[t>=0];
                for(l in t)
                {   sel[which(y[,3+l] > 0),k] = 1;
                }
            }
            Le = rbind(Le, sel);
        }
    }
    
    if(FALSE)
    {
    id = showcol;
    print(id);print(dim(X));
    }
id=1;

    if(remove0>=0)
    {   t=which(X[,id]<remove0);
        X=X[-t,];
        Lm=Lm[-t];Lc=Lc[-t];Ls=Ls[-t];
        if(max(estates) >= 0) Le = Le[-t];
    }
    if(sign < 0) X[,id]=1/X[,id];
    if(log10==T) X[,id]=log10(X[,id]+min(X[X[,id]>0,id]/1));
    
    let = c("-","o");
    if(max(estates) >= 0)
    {   L = apply(cbind(Lm, Lc, Ls, Le),1,function(x){paste(rev(let[x+1]),collapse=" ")});
    }
    else
    {   L = apply(cbind(Lm, Lc, Ls),1,function(x){paste(rev(let[x+1]),collapse=" ")});
    }
    ll = unique(L);
    xx = yy = mm = NULL;
    for(i in 1:length(ll))
    {   t = which(L == ll[i]);
        xx = c(xx, rep(i, length(t)));
        yy = c(yy, X[t, id]);
        if(log10==T) { mm[i] = mean(X[t,id]); }
        else { mm[i] = median(X[t,id]); }
    }
    
    o=order(mm);
    if(orient==-1) { o = rev(o); }
    map=rep(0,length(ll))
    map[o]=1:length(ll);
    r=boxplot(yy~map[xx],names=ll[o],las=2, outline=F,ylim=range(yy));
    gg=oo=NULL;
    for(i in unique(r$group))
    {   t=which(r$group == i);
        u=unique(r$out[t]);
        oo=c(oo,u);
        gg=c(gg,rep(i,length(u)));
    }
    points(gg,oo);
    
    if(max(estates) >= 0) { return(cbind(X, Lm, Lc, Ls, Le)); }
    else { return(cbind(X,Lm, Lc, Ls)); }
}

openPDFdevice<-function(fname, height=8, width=8)
{
    pdf(file=fname, height=height, width=width, onefile=TRUE, family='Helvetica', paper='a4r', pointsize=12);
}

plotPieRadius<-function(x, r, cols, labels = names(x), edges = 200, startpos = 0, std = T, orient=-1, ma=-1)
{
    if (is.null(labels)) { labels <- as.character(seq_along(x)); }
    if(std) x=x/sum(x);
    cx=c(0,cumsum(x));
    nx = length(x);
    if(ma<0) { ma=max(r); }
    plot(-1000,-1000,xlim=c(-1,1)*ma, ylim=c(-1,1)*ma);
    for (i in 1:nx) 
    {   
        n <- max(2, floor(edges * x[i]))
            t2p <- orient * 2 * pi * seq(cx[i], cx[i + 1], length = n) + startpos
            xc <- c(cos(t2p) * r[i], 0)
            yc <- c(sin(t2p) * r[i], 0)
            polygon(xc, yc, col = cols[i])
            t2p <- orient * 2 * pi * mean(cx[i + 0:1]) + startpos
            xc <- cos(t2p) * r[i]
            yc <- sin(t2p) * r[i]
        lines(c(1,1.05)*xc,c(1,1.05)*yc);
        text(1.1 * xc, 1.1 * yc, labels[i], xpd = TRUE, adj = ifelse(xc < 0, 1, 0));
    }
}

plotCTCFcompare<-function()
{   pref=c("ctcf_manualset","ctcf_nonencode");
    names=c("mineCompleten4","chrHMMf","segwayf");
    par(mfrow=c(2,3));
    par(mar=c(2,2,2,2));
    
    for(i in 1:length(pref))
    {   ma=NULL;
        for(j in 1:length(names))
        {   x=read.table(paste0(pref[i],".collapse.",names[j]));
            p=c(x[1,5],diff(x[,5]));
            p=p/sum(p);
            r=x[,8];
            p=p[-length(p)];
            r=r[-length(r)];
            if(j==1) { ma=max(r); }
            cols=array(255,dim=c(length(p),3));
            cols[,-j] = rep((1:length(p))/length(p)*200,2);
            cols=rgb(cols,maxColorValue=255);
            plotPieRadius(p,r,cols,std=F,ma=ma);
        }
    }
}

geneReg_each<-function(x)
{
y=(c(t(log2(x[2:13]+0.01))))
z=t(array(c(t(x[-(1:13)])),dim=c(gN,6)));z=z[rep(1:dim(z)[1],each=2),]
#z=t(t(z)-apply(z,2,mean));
k<<-k+1;
print(k);
#r=svd(z);
#z=z%*%r$v[,1:pcn];
#f=lm(y~z)$fitted;
Z=cbind(1,z);
f=Z%*%solve(t(Z)%*%Z+diag(0.01,dim(Z)[2]))%*%t(Z)%*%y;
return(var(f)/var(y));
}

compareGeneReg_each<-function(scut=2, pn=1, fout=NULL, n=1000000, suf="", expfile = "a6.txt")
{
    if(length(fout)>0)
    {   openPDFdevice(paste(fout,".pdf",sep=""), height=4, width=16);
    }
    names=c("IDEAS","ChromHMM","Segway");
    pcn<<-pn;
    par(mfrow=c(1,3));

    g=read.table(paste0("GeneExp/",expfile)); #a6 is all gene exp, c6 is tss2k exp
    g=as.matrix(g);
    mm=read.table("GeneExp/gencode7.filter500.map");
    mm=as.matrix(mm)[,1];
    g=g[mm,];
    s=apply(log2(g+0.01),1,sd);
    s=rep(s,each=100);
    s=s[1:n];
    for(i in 0:2)
    {   if(i>0) { x=read.table(paste0("gencode7exp.",i),nrows=n); }
        else { x=read.table(paste0("gencode7exp",suf,".",i),nrows=n); } 
        x[,2:13]=g[rep(1:(n/100),each=100),];
        #s=apply(log2(x[,2:13]+1),1,sd);
        gN<<-(dim(x)[2]-13)/6;
for(k in 0:99)
{ t=which(s>scut & x[,1]==k);
nx=t(array(t(x[t,-(1:13)]),dim=c(gN, length(t)*6)));
for(j in 1:dim(nx)[2]) nx[,i]=nx[,i]-mean(nx[,i]);
r=svd(nx);
z=nx%*%r$v[,1:pcn];
z=t(array(t(z),dim=c(6*pcn,length(t))));
x[t,13+1:(pcn*6)] = z;
}
x=x[,1:(13+pcn*6)]
gN<<-pcn;
        k<<-0;
        r2=apply(x[s>scut,],1,geneReg_each);
        m=rep(0,100);
        for(j in unique(x[s>scut,1]))
        {   m[j+1]=mean(r2[which(x[s>scut,1]==j)]);
        }
        boxplot(r2~x[s>scut,1], main=names[i+1]);
        lines(1:length(m)-1,m,col=i+2);
    }
    if(length(fout)>0) { dev.off(); }
}

compareGeneReg_all<-function(scut=2, fout=NULL, n=1000000, suf="", expfile="a6.txt")
{
    if(length(fout)>0)
    {   openPDFdevice(paste(fout,".pdf",sep=""), height=8, width=16);
    }
    pref=c("mineCompleten4","chrHMM","segway");

    g=read.table(paste0("GeneExp/",expfile));
    g=as.matrix(g);
    mm=read.table("GeneExp/gencode7.filter500.map");
    mm=as.matrix(mm)[,1];
    g=g[mm,];
    s=apply(log2(g+0.01),1,sd);
    s=rep(s,each=100);
    s=s[1:n];

#library(e1071);
#library("randomForest");
    for(i in 0:2)
    {   if(i>0) { x=read.table(paste0("gencode7exp.",i),nrows=n); }
        else { x=read.table(paste0("gencode7exp",suf,".",i),nrows=n); } 
        #x=read.table(paste0("gencode7exp.",i),nrows=n);
        x[,2:13]=g[rep(1:(n/100),each=100),];
#nx=x[,-(1:13)];
#for(j in 2:99)
#{  t=(1:(n/100)-1)*100+j;
#   for(k in c(-1,1))
#   {   
#       x[t,-(1:13)] = x[t,-(1:13)]+nx[t+k,];
#   }
#}
        #s=apply(log2(x[,2:13]+1),1,sd);
        gN <<- (dim(x)[2]-13)/6;
        scut<<-scut;
        #r2=apply(x[s>scut,],1,geneReg_each);
        #rr=t(array(r2,dim=c(100,length(r2)/100)));
        #m=cutree(hclust(dist(rr),method="ward"),k=10);
        a=b=coef=NULL;
        for(j in 0:99)
        {   t=which(s>scut & x[,1]==j);
            y=c(t(log2(x[t,2:13]+1)));
            z=t(array(c(t(x[t,-(1:13)])),dim=c(gN,length(t)*6)))[rep(1:(length(t)*6),each=2),];

        #   l=dim(z)[1];t1=sample(l,size=round(0.9*l));t2=(1:l)[-t1];
        #   beta=lm(y[t1]~z[t1,])$coefficients;
        #   if(length(which(is.na(beta)==T))>0) { beta[which(is.na(beta)==T)] = 0; }
        #   f=cbind(1,z[t2,])%*%beta;
        #   a[j+1]=var(f)/var(y[t2]);

        #   f=rep(0,length(y));
        #   for(k in 1:max(m))
        #   {   tt=which(m==k);
        #       ttt=rep((tt-1)*12,each=12)+rep(1:12,length(tt));
        #       f[ttt]=lm(y[ttt]~z[ttt,])$fitted;   
        #   }
        #   a[j+1]=var(f)/var(y);
            
            ta=rep(0,6);
            tcoef = NULL;
            for(k in 1:6)
            {   tt=rep((1:length(t)-1) * 12, each=2) + rep((k-1)*2+1:2,length(t));
                r=(lm(y[tt]~z[tt,]));
                #r=0;
                #for(kkk in 1:10)
                #{  ttt=sample(1:length(tt),size=1000);
            #       ff=svm(y[tt[ttt]]~z[tt[ttt],])$fitted;
            #       r=r+var(ff)/var(y[tt[ttt]]);
            #   }
            #   ta[k]=r/10;
            #   print(c(j,k,r/10));
                ta[k]=summary(r)$r.squared;
                tcoef=c(tcoef, as.numeric(r$coef));
#ty=y[tt];tz=z[tt,];
#data=data.frame(ty=ty,tz=tz);
#ttt=sample(1:length(tt),size=length(tt)*9/10);
#r=randomForest(ty~.,data=data,ntree=50,subset=ttt);
#f=predict(r,data);
#print(c(j,k,ta[k],1-var(ty[ttt]-f[ttt])/var(ty[ttt]),1-var(ty[-ttt]-f[-ttt])/var(ty[-ttt])));
#ta[k+6]=1-var(ty[ttt]-f[ttt])/var(ty[ttt]);
#ta[k+12]=1-var(ty[-ttt]-f[-ttt])/var(ty[-ttt]);
            }
            a=rbind(a,ta);
            coef=rbind(coef,tcoef);

            #r=summary(lm(y~z));
            #a[j+1]=r$r.squared;
#           b=rbind(b,r$coefficients[,1]);
        }
        write.table(cbind(a,b),paste0(pref[i+1],".geneReg_all.sum"),quote=F,row.names=F,col.names=F);
        write.table(coef,paste0(pref[i+1],".geneReg_all.coef"), quote=F,row.names=F,col.names=F);
        if(i==0)
        {   plot(a[,1],ylim=c(0,1),type="o",col=i+2);
        }
        else
        {   lines(a[,1],type="o",col=i+2);
        }
    }
    if(length(fout)>0) { dev.off(); }
}

plotBox<-function(x, outlier = NULL, cols = "white", ylim=NULL, border="black", fence=TRUE)
{
    ncols=cols;if(is.na(border)==FALSE) { ncols[1]="black"; }
    x=as.matrix(x);
    n = length(x) / 5;
    if(n == 1) { x = array(c(x), dim=c(1,5)); }
    ooo = NULL;
    outlier = as.character(outlier);
    if(length(outlier)>0)
    {   for(i in 1:n)
        {   ooo = c(ooo, as.numeric(unlist(strsplit(outlier[i],","))));
        }
    }
    if(length(ylim)==0)
    {   ylim=c(min(min(x[,(2-as.integer(fence)):(4+as.integer(fence))]),min(ooo)), max(max(x[,(2-as.integer(fence)):(4+as.integer(fence))]),max(ooo)));
    }
    for(i in 1:n)
    {   if(i == 1)
        {   plot(-100,-100,xlim=c(1,n),ylim=ylim);
        }
        ooo = NULL;
        if(length(outlier)>0)
        {   ooo = as.numeric(unlist(strsplit(outlier[i],",")));
        }
        if(fence) lines(rep(i,2),c(x[i,1],x[i,5]),col=ncols[(i-1)%%length(cols)+1]);
        dd=0.5*as.integer((i%%3)!=1);
        rect(i-0.5-dd,x[i,2],i+0.5-dd+2*as.integer(i%%3==1),x[i,4],col=cols[(i-1)%%length(cols)+1],border=border);
        if((i%%3)!=1)lines(c(i-0.5,i+0.5)-dd,rep(x[i,3],2),lwd=2);
        if(length(ooo)>0) points(rep(i,length(ooo)),ooo, col=ncols[(i-1)%%length(ncols)+1]);
    }
}

stateGWAS<-function(threshold=0.05, usep=F, cellnames=paste("v",1:6,".bed.gwasv",sep=""))
{   
    rt=NULL;
    I=X=P=F=NULL;
    for(i in 1:length(cellnames))
    {   x=read.table((cellnames[i]));
        X=rbind(X,x);
        I=c(I,rep(i,dim(x)[1]));
        p=apply(as.matrix(x[,c(2,6:dim(x)[2])]),1,function(x){length(which(x[-1]>=x[1]))/(length(x)-1);});
        m=apply(as.matrix(x[,-(1:5)]),1,mean);
        f=(x[,2] + 0.5) /(m+0.5);
        P=c(P,p);
        F=c(F,f);
        #fdr=p.adjust(p,"fdr");
        #t=which(fdr<=threshold);
        #if(length(t)>0)
        #{  row=cbind(i,x[t,2:3], p[t], fdr[t], f[t], x[t,5]);
    #       rt=rbind(rt, row);  
    #   }
    }
    fdr=p.adjust(P,"fdr");
    t=which(fdr<=threshold);
    if(usep==T) { t=which(P<=threshold)}
    rt=cbind(I[t],X[t,2:3],P[t],fdr[t], F[t], X[t,5]);
    dnames=sort(unique(as.character(rt[,7])));
    l=length(dnames);
    Pp = Fdr = Fold = array(0,dim=c(l,length(cellnames)));
    for(i in 1:length(cellnames))
    {   #x=read.table((cellnames[i]));
        #p=apply(as.matrix(x[,c(2,6:dim(x)[2])]),1,function(x){length(which(x[-1]>=x[1]))/(length(x)-1);});
        #m=apply(as.matrix(x[,-(1:5)]),1,mean);
        #f=(x[,2] + 0.5) /(m+0.5);
        #fdr=p.adjust(p,"fdr");
        for(j in dnames)
        {   t=which(as.character(X[,5])==j);
            tt=which(dnames==j);
            for(k in t)
            {   Pp[tt,I[k]] = X[k,2]/X[k,3];
                Fdr[tt,I[k]] = fdr[k];
                Fold[tt,I[k]] = F[k];
            }
        }    
    }
    rownames(Pp)=rownames(Fdr)=rownames(Fold)=dnames;
    colnames(Pp)=colnames(Fdr)=colnames(Fold)=cellnames;
    RT=NULL;
    RT$Pp=Pp;
    RT$Fdr=Fdr;
    RT$Fold=Fold;
    RT$rt=rt;
    return(RT);
}

heatmapwrap<-function(rr, simrange=F, celltext=-1, method="complete", cols=c("white","blue","black","red","yellow"), rowv=T,colv=NA, rang=NULL, colsep=NULL,rowsep=NULL, cx=NULL,rowcol=NULL,colcol=NULL)
{
    if(length(rang)==0) { rang = range(rr); }
    if(simrange==T) { rang=c(-1,1)*max(abs(rang));}
    colors=0:100/100*(rang[2]-rang[1])+rang[1];
        my_palette=colorRampPalette(cols)(n=100);
    myhclust<-function(x){return(hclust(x,method=method))};

    if(length(rowsep)>0)
    {   rcl = NULL;
        for(i in 1:length(rowsep)) 
        {   if(i==1) { rcl = c(rcl, rep(i+1,rowsep[i])); }
            else { rcl = c(rcl, rep(i+1,rowsep[i] - rowsep[i-1])); }
        }
        rcl = c(rcl, rep(1, dim(rr)[1]-rowsep[length(rowsep)]));
        rcl = c("grey",rainbow(length(rowsep)))[rcl];
    }
    if(length(rowcol)>0 | length(colcol)>0)
    {   if(length(rowcol)>0 & length(colcol)>0)
        {   o=heatmap.3(rr, breaks=colors, col=my_palette,hclustfun=myhclust,Colv=colv,Rowv=rowv,colsep=colsep,rowsep=rowsep, RowSideColors=rowcol, ColSideColors=colcol); 
        } else if(length(rowcol)>0)
        {   o=heatmap.3(rr, breaks=colors, col=my_palette,hclustfun=myhclust,Colv=colv,Rowv=rowv,colsep=colsep,rowsep=rowsep, RowSideColors=rowcol);
        } else
        {   o=heatmap.3(rr, breaks=colors, col=my_palette,hclustfun=myhclust,Colv=colv,Rowv=rowv,colsep=colsep,rowsep=rowsep, RowSideColors=rowcol);
        }
    } else if(celltext >= 0)
        {       if(length(cx)==0) cx = round(rr * (10^celltext))/(10^celltext);
                if(length(rowsep)==0)
        {   o=heatmap.3(rr, cellnote = cx, notecol = array(c("white","black")[2+(sign(rr)-1)/2],dim=dim(rr)), breaks=colors, col=my_palette,hclustfun=myhclust,Colv=colv,Rowv=rowv,colsep=colsep,rowsep=rowsep);
        }
        else
        {   
            o=heatmap.3(rr, cellnote = cx, notecol = array(c("white","black")[2+(sign(rr)-1)/2],dim=dim(rr)), breaks=colors, col=my_palette,hclustfun=myhclust,Colv=colv,Rowv=rowv,colsep=colsep,rowsep=rowsep, RowSideColors=rcl);
        }
        }
        else
        {       if(length(rowsep)==0)
        {   o=heatmap.3(rr, breaks=colors, col=my_palette,hclustfun=myhclust,Colv=colv,Rowv=rowv,colsep=colsep,rowsep=rowsep);
        }
        else
        {   o=heatmap.3(rr, breaks=colors, col=my_palette,hclustfun=myhclust,Colv=colv,Rowv=rowv,colsep=colsep,rowsep=rowsep, RowSideColors=rcl);
        }
    }
    
    return(o);
}

plotStatePair<-function(log=F,fold=F,method="complete")
{
    names=c("mineCompleten4","chrHMMf","segwayf");
    rang=c(0,0);
    for(i in 1:3)
    {   x=read.table(paste(names[i],".pair",sep=""));
        x=as.matrix(x);
        p=apply(x,1,sum)/5;
        e=t(t(p))%*%t(p)/sum(p);
        e1=t(t(1-p/sum(p)))%*%t(1-p/sum(p));
        rr=(x/5-e)/sqrt(e)/sqrt(e1)/sqrt(sum(p));
        if(fold)
        {   rr= (x/5+1)/(e+1);
            rr=log2(abs(rr))*sign(rr);
        }
        colnames(rr)=rownames(rr)=read.table(paste(names[i],".labmap",sep=""))[,2];
        l=dim(rr)[1];
        if(i==1) {l=24; }
        rr=rr[1:l,1:l];
        if(log & !fold)
        {   rr=log2(abs(rr*100)+1)*sign(rr);
        }
        if(rang[1]>min(rr)) {rang[1]=min(rr);}
        if(rang[2]<max(rr)) {rang[2]=max(rr);}
rang=c(-1,1)*max(abs(rang));
        cols=c("blue","dark blue","black","yellow","yellow");
        colors=0:100/100*(rang[2]-rang[1])+rang[1];
            my_palette=rainbow(100);#colorRampPalette(cols)(n=100);
            myhclust<-function(x){return(hclust(x,method=method))};
        heatmap.3(rr, breaks=colors, col=my_palette,hclustfun=myhclust,symm=T);
    }
}

plotStateCell<-function()
{   names=c("mineCompleten4","chrHMMf","segwayf");
        for(i in 1:3)
    {   x=read.table(paste(names[i],".statehist",sep=""));
        ns=x[,1];
        x=as.matrix(x[,-1]);
x=x+100;
        x=x/apply(x,1,sum);
        rownames(x)=ns;
        if(i==1) { x=x[-25,]; }

x=log2(x*6);
rang=range(x);
rang=c(-1,1)*2.;
cols=c("blue","black","yellow");
colors=0:100/100*(rang[2]-rang[1])+rang[1];
my_palette=colorRampPalette(cols)(n=100);
        heatmap.3(x,Colv=NA,breaks=colors,col=my_palette)#
    }
}

getTssTestoInterval<-function(inv, n=25, glist=NULL,weight=F)
{   g=read.table("GeneExp/GenCodeV7hg19_filter500.bed");
#t=which(g[,3]-g[,2]>1000);
#   g=g[t,];
    G=NULL;
    inv=as.matrix(inv);
    if(weight==T)
    {   G=read.table("G.txt");G=log(as.matrix(G)+0.01);
#G=G[t,];
    }
    if(length(glist)>0) {g=g[glist,];}
    tss=g[,2]; tes=g[,3];
    t=which(g[,6]=="-");
    tss[t] = g[t,3];
    tes[t] = g[t,2];
    
    tsz = step = 100000 ^ (1 / n);
    wsz = rep(200, n);
    for(i in 1:n)
    {   wsz[i] = max(wsz[i], round(tsz));
        tsz = tsz * step;
    }
    sidesz = sum(wsz);
    osz=wsz;
    wsz = cumsum(wsz);

    ntss=ntes=r2tss=r2tes= counttss = counttes = rep(0, 4 * n);

    chrname = paste("chr",1:22,sep="");
    chrname = c(chrname, "chrX", "chrY");
    l=dim(inv)[1];
    W = W2 = N = N2 = 0;
    for(i in 1:l)
    {   t = which(as.character(g[,1]) == chrname[inv[i,1]] & tss > inv[i,2] - sidesz & tss < inv[i,3]+sidesz);
        tw = rep(1,length(t));
        tr = rep(0,length(t));
        if(weight==T & length(t)>0)
        {   if(sd(inv[i,6:11])==0) { tw = rep(1e-5, length(t)); }
            else
            {
                for(j in 1:length(t))
                {#  print(c(t[j],inv[i,6:11]));
                    r=summary(lm(G[t[j],]~as.factor(rep(inv[i,6:11],each=2))));
                    f=r$f;
                    tr[j]=r$r.square;
                    tw[j]=-log10(1-pf(f[1],f[2],f[3])+1e-16);
if(tss[t[j]]<inv[i,2])
{   k=inv[i,2]-tss[t[j]];
    tt=min(which(k<wsz));
    ii=n-tt+1;if(tss[t[j]]>tes[t[j]]){ii=4*n-ii+1;}
}
else if(tss[t[j]]>inv[i,3])
{   k=tss[t[j]]-inv[i,3];
    tt = min(which(k < wsz));
    ii=3*n+tt;if(tss[t[j]]>tes[t[j]]){ii=4*n-ii+1;}
}
else
{   tsz = (inv[i,3]-inv[i,2] + 1) / 2 / n;
    k = as.integer((tss[t[j]] - inv[i,2]) / tsz) + 1;
    ii=n+k;if(tss[t[j]]>tes[t[j]]){ii=4*n-ii+1;}
}
print(c(ii,var(G[t[j],]),tr[j]));
                }
            }
        }
        a = which(tss[t] < inv[i,2]);
        if(length(a)>0)
        {   for(j in a)
            {   jj=t[j];
                k = inv[i,2] - tss[jj];
                tt = min(which(k < wsz));
                ii=n-tt+1;if(tss[jj]>tes[jj]){ii=4*n-ii+1;}
                w=tw[j];
                counttss[ii] = counttss[ii] + w / osz[tt];
                r2tss[ii]=r2tss[ii]+tr[j];
                ntss[ii]=ntss[ii]+1;
                W = W + w; N = N + 1;
            }
        }
        b = which(tss[t] > inv[i,3]);
        if(length(b)>0)
        {   for(j in b)
            {   jj=t[j];
                k = tss[jj] - inv[i,3];
                tt = min(which(k < wsz));
                ii=3*n+tt;if(tss[jj]>tes[jj]){ii=4*n-ii+1;}
                w=tw[j];
                counttss[ii] = counttss[ii] + w / osz[tt];
                r2tss[ii]=r2tss[ii]+tr[j];
                ntss[ii]=ntss[ii]+1;
                W = W + w; N = N + 1;
            }
        }
        m = which(tss[t] >= inv[i,2] & tss[t] <= inv[i,3]);
        tsz = (inv[i,3]-inv[i,2] + 1) / 2 / n;
        if(length(m)>0)
        {   for(j in m)
            {   jj=t[j];
                k = as.integer((tss[jj] - inv[i,2]) / tsz) + 1;
                ii=n+k;if(tss[jj]>tes[jj]){ii=4*n-ii+1;}
                w=tw[j];
                counttss[ii] = counttss[ii] + w / tsz;
                r2tss[ii]=r2tss[ii]+tr[j];
                ntss[ii]=ntss[ii]+1;
                W = W + w; N = N + 1;
            }
        }

        t = which(as.character(g[,1]) == chrname[inv[i,1]] & tes > inv[i,2] - sidesz & tes < inv[i,3]+sidesz);
        tw = rep(1,length(t));
        tr = rep(0,length(t));
        if(weight==T & length(t)>0)
        {   if(sd(inv[i,6:11])==0) { tw = rep(1e-5, length(t)); }
            else
            {   for(j in 1:length(t))
                {#  print(c(t[j],inv[i,6:11]));
                    r=summary(lm(G[t[j],]~as.factor(rep(inv[i,6:11],each=2))));
                    f=r$f;
                    tr[j] = r$r.square;
                    tw[j]=-log10(1-pf(f[1],f[2],f[3])+1e-16);
                }
            }
        }
        a = which(tes[t] < inv[i,2]);
        if(length(a)>0)
        {   for(j in a)
            {   jj=t[j];
                k = inv[i,2] - tes[jj];
                tt = min(which(k < wsz));
                ii=n-tt+1;if(tss[jj]>tes[jj]){ii=4*n-ii+1;}
                w=tw[j];
                counttes[ii] = counttes[ii] + w / osz[tt];
                r2tes[ii]=r2tes[ii]+tr[j];
                ntes[ii]=ntes[ii]+1;
                W2 = W2 + w; N2 = N2 + 1;
            }
        }
        b = which(tes[t] > inv[i,3]);
        if(length(b)>0)
        {   for(j in b)
            {   jj=t[j];
                k = tes[jj] - inv[i,3];
                tt = min(which(k < wsz));
                ii=3*n+tt;if(tss[jj]>tes[jj]){ii=4*n-ii+1;}
                w=tw[j];
                counttes[ii] = counttes[ii] + w / osz[tt];
                r2tes[ii]=r2tes[ii]+tr[j];
                ntes[ii]=ntes[ii]+1;
                W2 = W2 + w; N2 = N2 + 1;
            }
        }
        m = which(tes[t] >= inv[i,2] & tes[t] <= inv[i,3]);
        if(length(m)>0)
        {   for(j in m)
            {   jj=t[j];
                k = as.integer((tes[jj] - inv[i,2]) / tsz) + 1;
                ii=n+k;if(tss[jj]>tes[jj]){ii=4*n-ii+1;}
                w=tw[j];
                counttes[ii] = counttes[ii] + w / tsz;
                r2tes[ii]=r2tes[ii]+tr[j];
                ntes[ii]=ntes[ii]+1;
                W2 = W2 + w; N2 = N2 + 1;
            }
        }

    }
    W=W/(N+1e-100);W2=W2/(N2+1e-100);
print(c(W,W2,N,N2));
    counttss = counttss / dim(inv)[1]/W;
    counttes = counttes / dim(inv)[1]/W2;

    r2tss = r2tss / (ntss+1e-10);
    r2tes = r2tes / (ntes+1e-10);
    
if(FALSE)
{
    lo1=loess(counttss~I(1:100),span=0.2,degree=1);
    lo2=loess(counttes~I(1:100),span=0.2,degree=1);
    plot(predict(lo1));
    lines(predict(lo2),type="h");
    lines(rep(n,2),c(0,1));
    lines(rep(3*n,2),c(0,1));
}
    pos=c(rev(wsz),rep(1,2*n),wsz);
    return(cbind(pos,counttss, counttes, r2tss, r2tes));
}

plotTssTestoInterval<-function(lcut=5000000,zcut=10)
{
    x=read.table("mineCompleten4.longregionTSSTES");
    x=as.matrix(x[0:118*2+1,-1])*1000000;
    x1=read.table("mineCompleten4.longregionTSSTES_wlogP");
    x1=as.matrix(x1[0:118*2+1,-1])*1000000;
    p=read.table("encodeComplete.All.n4.popsum");
    p=as.matrix(p);
    l=pat=NULL;
    lt=c("0","1","2","3","4","5");
    for(i in 0:118)
    {   t=which(p[,12]==i & p[,3]-p[,2]>1000);
        l=c(l,sum(p[t,3]-p[t,2]));
        a=p[which(p[,12]==i)[1],6:11];
        b=table(a);
        o=order(as.numeric(b),decreasing=T);
        na=a;
        for(j in 1:length(o))
        {   na[which(as.numeric(names(b)[o[j]])==a)]=j;
        }
        pat=c(pat,paste(lt[na],collapse=""));
    }
    z=apply(x,1,function(x){length(which(x==0))});

    t=which(l>=lcut & z<=zcut);
    
    y=y1=NULL;
    for(i in t)
    {   y=rbind(y,predict(loess(x[i,]~I(1:100),span=0.2,degree=1)));
        y1=rbind(y1,predict(loess(x1[i,]~I(1:100),span=0.2,degree=1)));
    }
    #y=y-apply(y,1,mean);
    #y1=y1-apply(y1,1,mean);
#   y=(y)/apply(y,1,sum);
#   y1=(y1)/apply(y1,1,sum);
    rownames(y)=pat[t];

#return(cbind(y,y1));
    heatmapwrap(log2((y1+1e-6)/(y+1e-6)),simrange=T, method="ward",colsep=100);
    return(y);
}

plotIntervalState<-function(state, plot=T)
{   x=read.table(paste("mineCompleten4.longregion_state",state,sep=""));
    x=as.matrix(x)
    x=t(t(x)/apply(x,2,sum));
    x=x/x[,26];
    nx=NULL;
    for(i in 1:24)
    {   nx=rbind(nx,predict(loess(x[,i]~I(1:100),span=0.2,degree=1)));
    }
    rownames(nx)=read.table("mineCompleten4.labmap")[-25,2];
    colnames(nx)=read.table("mineCompleten4.longregionTSSTES")[239,-1];

    if(plot) heatmapwrap(nx,method="ward",rang=c(0,2));

    return(nx);
}

pairEnrich_clusterheat<-function(rt, labs=NULL, K=8, method="ward", fout=NULL, pseudo=100)
{
    if(length(fout)>0)
    {   pdf(file=paste(fout,".pdf",sep=""), height=8, width=8, onefile=TRUE, family='Helvetica', paper='a4r', pointsize=12);
    }
    xl=labs;
#remove last column of rt, which is the sum of rows
    rt=as.matrix(rt[,-dim(rt)[2]]);
    if(length(labs)==0) xl=1:(sqrt(dim(rt)[2]))-1;

    sz=length(xl);
    totaln = apply(rt,2,sum)/5;
    rt=t(t(rt + pseudo)/apply(rt+pseudo,2,sum));
    rt=log2(rt/0.01)
    r=apply(rt,2,function(x){max(x)-min(x)});
    u=which(r>1);
    mm=cc=NULL;
    h=hclust(dist(t(rt[,u])),method=method);
    m=rep(1,dim(rt)[2]);
    m[u]=1+cutree(h,K-1);
    
    me=se=NULL;
    for(i in 1:K)
    {   if(length(which(m==i))>1)
        {   me=cbind(me,apply(rt[,m==i],1,mean));
            se=cbind(se,apply(rt[,m==i],1,sd)); 
        }
        else
        {   me=cbind(me,rt[,m==i]);
            se=cbind(se,0);
        }
    }
    me<-me;
    se<-se;

    dd=array(dim=c(sz,sz));
    ccc=array(dim=c(sz,sz));
    colnames(dd)=rownames(dd)=xl;
    k=1;
    for(i in 1:sz) 
    {   ff = 1:sz;
        for(j in ff)
        {   dd[i,j]= m[k];
            ccc[i,j]= round(log10(totaln[k]+1));; 
            k=k+1;
        }
    }
    mydist<-function(x)
    { d=array(0,dim=rep(dim(x)[1],2));
      for(i in 1:dim(x)[1])
        for(j in 1:i)
        {   d[i,j]=0;
            for(k in 1:dim(x)[2])
            {   d[i,j]=d[i,j]+sum((me[,x[i,k]]-me[,x[j,k]])^2);
            }
            d[i,j]=sqrt(d[i,j]);

            d[i,j]=(dim(x)[1]-length(which(x[i,]==x[j,])))
        }
        return(as.dist(d))
    }   
    
    source("heatmap3.R");
    #lmat=rbind(c(5,4,4,4),c(5,1,1,1),c(3,2,2,2),c(9,8,10,10),c(7,6,10,10));
    lmat=rbind(c(5,5,4,9,9,8),c(3,1,2,7,6,6),c(3,1,2,7,6,6),c(3,1,2,10,10,10),c(3,1,2,10,10,10),c(3,1,2,10,10,10));
    layout((lmat), heights = c(5,5,5,1,5,5), widths = c(5,1,10,3,5,5), respect = FALSE);
#5,1,10?
    rang=range(rt);
    rang=c(-1,1)*max(abs(rang));
    colors=0:100/100*(rang[2]-rang[1])+rang[1];
    my_palette=colorRampPalette(c("white","blue","black","red","yellow"))(n=100);
    myhclust<-function(x){return(hclust(x,method=method))};
    #heatmap.3(rt,Rowv=NA,breaks=colors,col=my_palette,hclustfun=myhclust,ColSideColors=rainbow(K)[m],drawlayout=F);

    heatmap.3(t(rt[,u]),Colv=NA,dendrogram="row",breaks=colors,col=my_palette,hclustfun=myhclust,RowSideColors=c(rainbow(K-1))[m[u]-1],drawlayout=F);
    
    if(length(ccc)>0)
    {   heatmap.3(dd,trace="none",distfun=mydist,breaks=0:K,col=c("black",rainbow(K-1)),drawlayout=F, cellnote=ccc, notecol="black", notecex=1,key=F);
    }
    else
    {   heatmap.3(dd,trace="none",distfun=mydist,breaks=0:K,col=c("black",rainbow(K-1)),drawlayout=F,key=F);
    }
    plot(-100,-100,xlim=c(1,dim(me)[1]), ylim=range(me));
    for(i in 2:K)
    {   lines(me[,i],col=rainbow(K-1)[i-1],lwd=1.5);
    }
    abline(0,0,lty=2)
    lines(rep(dim(rt)[1]/4+0.5,2),c(-10,10),lty=2);
    lines(rep(dim(rt)[1]-dim(rt)[1]/4-0.5,2),c(-10,10),lty=2);


    if(length(fout)>0)
    {   dev.off();
    }

    return(m);
}

stateGWAS_stateEnrich<-function(type = 0, back = 1)
{       Fold = Fdr = H = Hl = NULL;
        if(type == 0)
        {       Fold = read.table("mineCompleten4.varGwasp0.01.fold");
                Fdr = read.table("mineCompleten4.varGwasp0.01.fdr");
        }
        if(type == 1)
        {       Fold = read.table("chrHMM.varGwasp0.01.fold");
                Fdr = read.table("chrHMM.varGwasp0.01.fdr");
        }
        if(type == 2)
        {       Fold = read.table("segway.varGwasp0.01.fold");
                Fdr = read.table("segway.varGwasp0.01.fdr");
        }
    f = apply(Fold, 1, max);
    p = apply(Fdr, 1, min);
        #t = which(apply(Fold,1,max)>=6 & apply(Fdr,1,min)<=0.01);
    t = which(f >= 2 & p <= 0.05 & -log10(p+1e-8)+f>=10);
        Fold = Fold[t,];
        Fdr = Fdr[t,];
        disease = as.character(read.table("disease.txt")[,1]);

        a = c(2,1,5,6,15,3);
        t = NULL;
        for(i in 1:6)
        {       if(type == 0)
        {   x = read.table(paste("v",i,".bed.hist",sep=""));
                    H = cbind(H, x[,2]);
                    x = read.table(paste("long",a[i],".bed.hist",sep=""));
                    Hl = cbind(Hl, x[,2]);
        }
        if(type == 1)
        {   x = read.table(paste("c",i,".bed.hist",sep=""));
                    H = cbind(H, x[,2]);
        }
        if(type == 2)
        {   x = read.table(paste("s",i,".bed.hist",sep=""));
                    H = cbind(H, x[,2]);
        }
        }
if(back==0)
{
if(type == 0)
{   H = Hl = as.matrix(read.table("mineCompleten4.statehist")[,-1]);
}
else if(type == 1)
{   H = as.matrix(read.table("chrHMMf.statehist")[,-1]);
}
else
{   H = as.matrix(read.table("segwayf.statehist")[,-1]);
}
}
        gN = 25;
        gNc = 25;
        gNs = 55;
        rn = rownames(Fold);
        SE0 = SE1 = SE2 = FID = NULL;
    for(i in 1:dim(Fold)[1])
        {       did = which(rn[i] == disease)[1] - 1;
                fid = which(Fold[i,] == max(Fold[i,]))[1];
        FID = c(FID, fid);
                if(type == 0)
        {   fid = as.integer((fid + 1) / 2);
fid=min(fid,6);
            if(fid>6)
            {
                            x = y = read.table(paste("long00",".bed.dhits",sep=""));
            }
            else
            {   x = read.table(paste("v",fid,".bed.dhits",sep=""));
                            y = read.table(paste("long",a[fid],".bed.dhits",sep=""));
            }
                }
                if(type == 1)
                {       x = read.table(paste("c",fid,".bed.dhits",sep=""));
                }
                if(type == 2)
                {       x = read.table(paste("s",fid,".bed.dhits",sep=""));
                }
                u = which(x[,2] == 0 & x[,1] == did);
                h = hist(x[u,2+fid],breaks=0:gN-1,plot=F)[[2]];
                SE0 = rbind(SE0, h);#(h + 0.5) / sum(h + 0.5) / (H[1:gN,fid] + 0.5) * sum(H[1:gN,fid] + 0.5));
                if(type == 0)
                {       u = which(y[,2] == 0 & y[,1] == did);
                        hl = hist(y[u,2+fid],breaks=0:gN-1,plot=F)[[2]];
                        SE0 = rbind(SE0, hl);#(hl + 0.5) / sum(hl + 0.5) / (Hl[1:gN,fid] + 0.5) * sum(Hl[1:gN,fid] + 0.5));
                }
                u = which(x[,2] == 1 & x[,1] == did);
                h = hist(x[u,2+fid],breaks=0:gNc-1,plot=F)[[2]];
                SE1 = rbind(SE1, h);#(h + 0.5) / sum(h + 0.5) / (H[gN+1:gNc,fid] + 0.5) * sum(H[gN+1:gNc,fid] + 0.5));
                if(type == 0)
                {       u = which(y[,2] == 1 & y[,1] == did);
                        hl = hist(y[u,2+fid],breaks=0:gNc-1,plot=F)[[2]];
                        SE1 = rbind(SE1, hl);#(hl + 0.5) / sum(hl + 0.5) / (Hl[gN+1:gNc,fid] + 0.5) * sum(Hl[gN+1:gNc,fid] + 0.5));
                }
                u = which(x[,2] == 2 & x[,1] == did);
                h = hist(x[u,2+fid],breaks=0:gNs-1,plot=F)[[2]];
                SE2 = rbind(SE2, h);#(h + 0.5) / sum(h + 0.5) / (H[gN+gNc+1:gNs,fid] + 0.5) * sum(H[gN+gNc+1:gNs,fid] + 0.5));
                if(type == 0)
                {       u = which(y[,2] == 2 & y[,1] == did);
                        hl = hist(y[u,2+fid],breaks=0:gNs-1,plot=F)[[2]];
                        SE2 = rbind(SE2, hl);#(hl + 0.5) / sum(hl + 0.5) / (Hl[gN+gNc+1:gNs,fid] + 0.5) * sum(Hl[gN+gNc+1:gNs,fid] + 0.5));
                }
        }
        colnames(SE0) = read.table("mineCompleten4.labmap")[,2];
        colnames(SE1) = read.table("chrHMM.labmap")[,2];
        colnames(SE2) = read.table("segway.labmap")[,2];
        rownames(SE0) = rep(rn, each=1+as.integer(type == 0));
        rownames(SE1) = rep(rn, each=1+as.integer(type == 0));
        rownames(SE2) = rep(rn, each=1+as.integer(type == 0));

        rt = NULL;
        rt$SE0 = SE0;
        rt$SE1 = SE1;
        rt$SE2 = SE2;
    rt$H = H;
    rt$Hl = Hl;
    rt$FID = FID;
        return(rt);
}

plotStateGWAS_stateEnrich<-function(type = 0, back=1, rowv=F, colv=T)#SE, FID, H, rowv=F, colv=T)
{
    rt=stateGWAS_stateEnrich(type, back=back);
    if(type == 0)
    {   n = dim(rt$SE0)[1]/2;
        SE = cbind(rt$SE0[1:n*2-1,],rt$SE0[1:n*2,]);    
        FID = rt$FID;
        H = rbind(rt$H[1:25,],rt$Hl[1:25,]);
        o = read.table("mineCompleten4.varGwasp0.01.order");
        o = as.numeric(as.matrix(o));
        SE = SE[o,rep(1:24,each=2)+rep(c(0,25),24)];
        FID = FID[o];
        H = H[rep(1:24,each=2)+rep(c(0,25),24),];
if(FALSE)
{
x=read.table("mineCompleten4.gwas.statehist.cellspecific");
xn = as.character(x[,1]);
x=as.matrix(x[,-1]);
n=rownames(SE);
nx=NULL;
for(i in 1:length(n))
{   
    t=which(xn == n[i])[1];
    f=min(6,round(FID[i]/2 + 0.1));
    nx=rbind(nx, x[t,(f-1)*25+1:24]);
}
colnames(nx) = read.table("mineCompleten4.labmap")[1:24,2];
rownames(nx) = n;
SE=nx;
H=as.matrix(read.table("mineCompleten4.statehist.cellspecific")[1:24,-1]);
SE=SE[,rep(1:24,each=2)];
H=H[rep(1:24,each=2),];
}   
#H=as.matrix(read.table(paste(ccc,".statehist",sep=""))[1:24,-1]);
#H=H[rep(1:24,each=2),]; 
    }
    else if(type == 1)
    {   SE = rt$SE1;
        FID = rt$FID;
        H = rt$H[25+1:25,];
        o = read.table("chrHMM.varGwasp0.01.order");
        o = as.numeric(as.matrix(o));
        SE = SE[o,];
        FID = FID[o];
    }
    else 
    {   SE = rt$SE2;
        FID = rt$FID;
        H = rt$H[50+1:55,];
        o = read.table("segway.varGwasp0.01.order");
        o = as.numeric(as.matrix(o));
        SE = SE[o,];
        FID = FID[o];
    }

#below ignores cell type specificity, but still use cell types that are most relevant to the phenotype
cnn = c("mineCompleten4","chrHMM","segway");
ttt=type+1;ccc=cnn[ttt];ppn=25;
tt=read.table(paste(ccc,".gwas.statehist",sep=""));
print(dim(tt));
oo=read.table(paste(ccc,".varGwasp0.01.order",sep=""));nn=rownames(oo);oo=as.numeric(as.matrix(oo));
for(i in 1:length(nn))
{   t=which(as.character(tt[,1]) == nn[i])[1];
    if(ttt==1)
    {   f=min(6,round(FID[oo[i]]/2 + 0.1));
        SE[oo[i],1:24*2]=as.numeric(tt[t,1+(f-1)*25+1:24]);
        SE[oo[i],1:24*2-1]=SE[oo[i],1:24*2];
    }
    else
    {   f=FID[oo[i]];
        SE[oo[i],]=as.numeric(tt[t,1+(f-1)*ppn+1:ppn]);
    }
}
H=as.matrix(read.table(paste(ccc,".statehist",sep=""))[,-1]);
if(ttt==1) { H=H[rep(1:24,each=2),]; }

    SE=as.matrix(SE);
    H=as.matrix(H);
    repn = dim(SE)[1] / length(FID);
    fid = rep(FID,each=repn);
    if(max(fid)>=12) { fid = as.integer((fid+1)/2); }
    fid[fid>6]=1;
        
    x=t(apply(cbind(fid,SE),1,function(x){(x[-1]+0.5)/(sum(x[-1])*(H[,x[1]])/sum(H[,x[1]])+0.5);}));

#p=t(apply(cbind(fid,SE),1,function(x){1-pbinom(x[-1],sum(x[-1]),H[,x[1]]/sum(H[,x[1]]))}));
#p=c(p);p=p.adjust(p,method="fdr");
#p=array(p,dim=dim(x));
#x=x-log10(p+1e-8);

tx=x;
if(type==0)
{   tx=rbind(x[,0:23*2+1],x[,0:23*2+2]);
}
oo=hclust(dist(t(log(tx))),method="ward")$order;
if(type==0)
{   oo = rep((oo - 1) * 2,each=2)+rep(1:2,length(oo));
    x=x[,oo];
}
else
{   x=x[,oo];
}
colv=F;

#   heatmapwrap((log(x)),simrange=T,colv=colv,rowv=rowv,cols=c("blue","black","yellow"),rang=c(-2,2),method="ward");
    heatmapwrap(((x)),simrange=F,rang=c(0,20),cols=c("black","blue","green","red","yellow"),colv=colv,rowv=rowv);
    return(x);
}   

plotDiseaseEnrich<-function(type = 0, showcell=F)
{   cellnames=NULL;
    suff="gwas";
    if(type == 0)
    {   cellnames=paste(c("v1","long2","v2","long1","v3","long5","v4","long6","v5","long15","v6","long3","long00"),".bed.",suff,sep="");
    }
    else if(type == 1)
    {   cellnames=paste("c",1:6,".bed.",suff,sep="");
    }
    else if(type == 2)
    {   cellnames=paste("s",1:6,".bed.",suff,sep="");
    }
    rt=stateGWAS(0.01,usep=T,cellnames=cellnames);
    f=apply(rt$Fold,1,max);
    p=apply(rt$Fdr,1,min);
    t=which(f>=2 & p<=0.05 & -log10(p+1e-8)+f>=10);

    OO = oo = NULL;
    if(length(t)>0) 
    {   x=as.matrix(rt$Fold[t,]);#-log10(rt$Fdr[t,]+1e-8));
        nn = round(dim(x)[2]/(1+as.integer(type==0)) + 0.2);
        mm = apply(x, 1, function(x){which(x==max(x))[1]});
        nx = nsz = ncx = NULL;
        for(i in 1:nn)
        {   if(type == 0) { tt = which(mm == (i-1)*2+1 | mm == (i-1)*2+2); }
            else { tt = which(mm == i); }
            if(length(tt>1)) 
            {   o = hclust(dist(x[tt,]),method="ward")$order;
            }
            else { o = 1:length(o); }
            o = tt[o];
            nx = rbind(nx, x[o,]);
            OO = c(OO, o);
            nsz[i] = dim(nx)[1];
            if(showcell == T)
            {   ncx = rbind(ncx, rt$Fold[o,]);
            }
        }
        if(showcell==T)
        {   #heatmapwrap(nx,cols=c("black","blue","red","yellow"),rang=c(0,20),method="ward",celltext=0,cx=ncx);
            oo = heatmapwrap(nx,cols=c("black","blue","green","red","yellow"),rang=c(0,20),rowv=F,celltext=0,cx=ncx, rowsep=nsz[-length(nsz)]);
        }
        else
        {   #heatmapwrap(nx,cols=c("black","blue","red","yellow"),rang=c(0,20),method="ward");
            oo = heatmapwrap(nx,cols=c("black","blue","green","red","yellow"),rang=c(0,20),rowv=F, rowsep=nsz[-length(nsz)]);
        }
    }
    rt$order=OO;
    rt$sel=t;
    return(rt);
}
