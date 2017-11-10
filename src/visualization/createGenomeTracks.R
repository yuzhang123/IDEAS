#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("-a", "--adult_mort"), action="store", dest="adult_mort", type="integer", help="Adjustment rate for adult mortality"),
    make_option(c("-b", "--adult_accum"), action="store", dest="adult_accum", type="integer", help="Adjustment of DD accumulation (old nymph->adult)"),
    make_option(c("-c", "--egg_mort"), action="store", dest="egg_mort", type="integer", help="Adjustment rate for egg mortality"),
    make_option(c("-d", "--latitude"), action="store", dest="latitude", type="double", help="Latitude of selected location"),
    make_option(c("-e", "--location"), action="store", dest="location", help="Selected location"),
    make_option(c("-f", "--min_clutch_size"), action="store", dest="min_clutch_size", type="integer", help="Adjustment of minimum clutch size"),
    make_option(c("-i", "--max_clutch_size"), action="store", dest="max_clutch_size", type="integer", help="Adjustment of maximum clutch size"),
    make_option(c("-j", "--nymph_mort"), action="store", dest="nymph_mort", type="integer", help="Adjustment rate for nymph mortality"),
    make_option(c("-k", "--old_nymph_accum"), action="store", dest="old_nymph_accum", type="integer", help="Adjustment of DD accumulation (young nymph->old nymph)"),
    make_option(c("-o", "--output"), action="store", dest="output", help="Output dataset"),
    make_option(c("-p", "--oviposition"), action="store", dest="oviposition", type="integer", help="Adjustment for oviposition rate"),
    make_option(c("-q", "--photoperiod"), action="store", dest="photoperiod", type="double", help="Critical photoperiod for diapause induction/termination"),
    make_option(c("-s", "--replications"), action="store", dest="replications", type="integer", help="Number of replications"),
    make_option(c("-t", "--se_plot"), action="store", dest="se_plot", help="Plot SE"),
    make_option(c("-u", "--year"), action="store", dest="year", type="integer", help="Starting year"),
    make_option(c("-v", "--temperature_dataset"), action="store", dest="temperature_dataset", help="Temperature data for selected location"),
    make_option(c("-y", "--young_nymph_accum"), action="store", dest="young_nymph_accum", type="integer", help="Adjustment of DD accumulation (egg->young nymph)")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser, positional_arguments=TRUE)
opt <- args$options

#markcol=rbind(c(1,0,0),c(1,1,0),c(0,1,1),c(0,0,1));
stateColor<-function(statemean, markcolor=NULL)
{   
    if(length(markcolor)==0)
    {   markcolor=rep("",dim(statemean)[2]);
        markcolor[order(apply(statemean,2,sd),decreasing=T)]=hsv((1:dim(statemean)[2]-1)/dim(statemean)[2],1,1)
        markcolor=t(col2rgb(markcolor));
    }

    rg=apply(statemean,1,range);
    mm=NULL;
    for(i in 1:dim(statemean)[1])
    {   mm=rbind(mm,(statemean[i,]-rg[1,i])/(rg[2,i]-rg[1,i]+1e-10));
    }
    mm = mm^6; 
    mm = mm / (apply(mm, 1, sum)+1e-10);
    mycol=mm%*%markcolor;
    s=apply(statemean,1,max);
    s=(s-min(s))/(max(s)-min(s)+1e-10);
#s=s^1.5;
    
    h=t(apply(mycol,1,function(x){rgb2hsv(x[1],x[2],x[3])}));
    h[,2]=h[,2]*s;
    h=apply(h,1,function(x){hsv(x[1],x[2],x[3])});
    rt=cbind(apply(t(col2rgb(h)),1,function(x){paste(x,collapse=",")}),h);
    
    return(rt);
}

createTrack<-function(statefiles, genomefile, outpref, statecolor, header, statename=NULL)
{   message("Reading state file: ", appendLF=FALSE);
    library("data.table");
    genomesz = read.table(genomefile);
    g=NULL;
    for(i in statefiles)
    {   message(paste(i," ",sep=""),appendLF=F);
        #tg=as.matrix(read.table(i, comment="!", header=(length(header)==0)));
        tg=as.matrix(fread(i));
        t=NULL;
        for(j in 1:dim(genomesz)[1])
        {   t=c(t,which(tg[,2]==as.character(genomesz[j,1]) & as.numeric(tg[,4])>as.numeric(genomesz[j,2])));
        }   
        if(length(t)>0) { tg = tg[-t,]; }
        t=which(is.na(match(tg[,2], genomesz[,1]))==T);
        if(length(t)>0) { tg = tg[-t,]; }   
print(c(dim(g),dim(tg)));
        g=rbind(g,tg);
    }
    message("Done");
    uchr = sort(unique(as.character(g[,2])));
    g1=NULL;
    for(i in uchr)
    {   t=which(g[,2]==i);
        g1=rbind(g1,g[t[order(as.integer(g[t,3]))],]);
    }
    g=NULL;

    chr=as.character(g1[,2]);
    posst=as.numeric(g1[,3]);
    posed=as.numeric(g1[,4]);
    state=as.matrix(g1[,5:(dim(g1)[2]-1)]);
    if(length(statename)==0) statename=0:max(state);
    L=dim(g1)[1];
    n=dim(state)[2];
    if(length(header) > 0) colnames(g1) = header;
    cells=as.character(colnames(g1)[5:(dim(g1)[2]-1)]);
    g1=NULL;
    message("Generating tracks");
    options(scipen=999);

    tt = which(chr[2:L]!=chr[2:L-1]);
    tt = c(tt,which(posst[2:L]!=posed[2:L-1]));
    tt = sort(unique(tt));

    for(i in 1:n)
    {   tstate = state[,i];
        #print(c(i,L,length(tstate),length(chr),length(posst),length(posed)));

        t=c(tt,which(tstate[2:L]!=tstate[2:L-1]));
        t=sort(unique(t));
        t0=c(0,t)+1;
        t=c(t,L);
        np=cbind(chr[t],posst[t0],posed[t],tstate[t]);

        #print("make track");
        x = cbind(np[,1:3],statename[as.integer(np[,4])+1],1000,".",np[,2:3],statecolor[as.numeric(np[,4])+1]);
        write.table(as.matrix(x),paste(outpref,i,"bed1",sep="."),quote=F,row.names=F,col.names=F);
print(x[1,]);
        #x = apply(x,1,function(x){paste(x,collapse="\t")});
        #write.table(x,paste(outpref,i,"bed",sep="."),quote=F,row.names=F,col.names=F);
        #system(paste("sort-bed ", outpref, ".", i, ".bed > ", outpref, ".", i, ".bed1", sep=""));
        system(paste("bedToBigBed ", outpref, ".", i, ".bed1 ", genomefile, " ", outpref, ".",i,".bb",sep=""));
        #system(paste("rm ", paste(outpref, i,"bed",sep=".")));
        system(paste("rm ", paste(outpref, i,"bed1",sep=".")));
    }
    return(cells);
}

#cellinfo: shortid match with those in states, long id, cell description, text color
run<-function(statefiles, hubid, genomeid, genomefile, statecolor, targetURL="http://bx.psu.edu/~yuzhang/tmp/", trackfolder=NULL, hubname = NULL, cellinfo = NULL, header=NULL,statename=NULL)
{
#   if(file.exists("bedToBigBed") == FALSE) { message("Cannot find bedToBigBed in the folder where you run this function."); return(0); }
#   if(file.exists("sort-bed") == FALSE) { message("Cannot find sort-bed in the folder where you run this function."); return(0); }
    if(length(hubname) == 0) hubname = hubid;
    if(length(trackfolder) == 0) trackfolder = paste("tracks_", hubid, "/", sep="");
    if(substring(trackfolder, nchar(trackfolder)) != "/") trackfolder = paste(trackfolder, "/", sep="");
    dir.create(trackfolder,showWarnings=FALSE);

    cells=createTrack(statefiles, genomefile, paste(trackfolder, hubid, sep=""), statecolor, header=header,statename=statename);
    if(length(cellinfo) == 0)
    {   cellinfo = cbind(cells, cells, cells, "#000000");
        cellinfo = array(cellinfo, dim=c(length(cells),4));
    }
    cellinfo = as.matrix(cellinfo);

    trackDb=NULL;
    for(i in 1:length(cells))
    {   ii=which(cells[i] == cellinfo[,1]);
        if(length(ii)==0) next;
        ii=ii[1];
        trackDb=c(trackDb, paste("track bigBed", i, sep=""));
        trackDb=c(trackDb, paste("priority",ii));
        trackDb=c(trackDb, "type bigBed 9 .");
        trackDb=c(trackDb, "itemRgb on");
        trackDb=c(trackDb, "maxItems 100000");
        trackDb=c(trackDb, paste("bigDataUrl ", targetURL, hubid,".",i,".bb",sep=""));
        trackDb=c(trackDb, paste("shortLabel", cellinfo[ii,2]));
        trackDb=c(trackDb, paste("longLabel", paste(hubname, cellinfo[ii,3])));
        trackDb=c(trackDb, paste("color",paste(c(col2rgb(cellinfo[ii,4])),collapse=",")));
        trackDb=c(trackDb, "visibility dense");
        trackDb=c(trackDb, ""); 
    }

    write.table(trackDb, paste(trackfolder, "trackDb_", hubid, ".txt",sep=""),quote=F,row.names=F,col.names=F);

    write.table(c(paste("genome", genomeid), paste("trackDb trackDb_", hubid, ".txt", sep="")), paste(trackfolder, "genomes_", hubid, ".txt",sep=""), quote=F,row.names=F,col.names=F);

    write.table(c(paste("hub", hubid), paste("shortLabel", hubid), paste("longLabel", hubname), paste("genomesFile genomes_", hubid, ".txt", sep=""), "email yzz2 at psu.edu"), paste(trackfolder, "hub_", hubid, ".txt",sep=""), quote=F,row.names=F,col.names=F);

    return(1);
}

createHeatmap<-function(parafile, statecolor = NULL, markcolor = NULL, cols=c("white","dark blue"), fout=NULL)
{   x=read.table(parafile,comment="!",header=T);
    k=dim(x)[2];
    l=dim(x)[1];
    p=(sqrt(9+8*(k-1))-3)/2;
    m=as.matrix(x[,1+1:p]/x[,1]);
    marks=colnames(m);
    rownames(m)=paste(1:l-1," (",round(x[,1]/sum(x[,1])*10000)/100,"%)",sep="");

    if(length(fout)!=0)
    {   pdf(fout);  }   
    par(mar=c(6,1,1,6));
    rg=range(m);
    colors=0:100/100*(rg[2]-rg[1])+rg[1];
        my_palette=colorRampPalette(cols)(n=100);
    defpalette=palette(my_palette);

    plot(NA,NA,xlim=c(0,p+0.7),ylim=c(0,l),xaxt="n",yaxt="n",xlab=NA,ylab=NA,frame.plot=F);
    axis(1,at=1:p-0.5,labels=colnames(m),las=2);
    axis(4,at=1:l-0.5,labels=rownames(m),las=2);
    rect(rep(1:p-1,l),rep(1:l-1,each=p),rep(1:p,l),rep(1:l,each=p),col=round((t(m)-rg[1])/(rg[2]-rg[1])*100));#,border=NA);

    if(length(statecolor)==0)
    {   if(length(markcolor)==0)
        {   markcolor=t(col2rgb(terrain.colors(ceiling(p))[1:p]));  
            for(i in 1:length(marks))
            {   if(regexpr("h3k4me3",tolower(marks[i]))>0)
                {   markcolor[i,]=c(255,0,0);   }
                if(regexpr("h3k4me2",tolower(marks[i]))>0)
                {   markcolor[i,]=c(250,100,0); }
                if(regexpr("h3k4me1",tolower(marks[i]))>0)
                {   markcolor[i,]=c(250,250,0); }
                if(regexpr("h3k36me3",tolower(marks[i]))>0)
                {   markcolor[i,]=c(0,150,0);   }
                if(regexpr("h2a",tolower(marks[i]))>0)
                {   markcolor[i,]=c(0,150,150); }
                if(regexpr("dnase",tolower(marks[i]))>0)
                {   markcolor[i,]=c(0,200,200); }
                if(regexpr("h3k9ac",tolower(marks[i]))>0)
                {   markcolor[i,]=c(250,0,200); }
                if(regexpr("h3k9me3",tolower(marks[i]))>0)
                {   markcolor[i,]=c(100,100,100);   }
                if(regexpr("h3k27ac",tolower(marks[i]))>0)
                {   markcolor[i,]=c(250,150,0); }
                if(regexpr("h3k27me3",tolower(marks[i]))>0)
                {   markcolor[i,]=c(0,0,200);   }
                if(regexpr("h3k79me2",tolower(marks[i]))>0)
                {   markcolor[i,]=c(200,0,200); }
                if(regexpr("h4k20me1",tolower(marks[i]))>0)
                {   markcolor[i,]=c(50,200,50); }
                if(regexpr("ctcf",tolower(marks[i]))>0)
                {   markcolor[i,]=c(200,0,250); }
            }
        }
        statecolor=stateColor(m,markcolor)[,2];
    }
    rect(rep(p+0.2,l),1:l-0.8,rep(p+0.8,l),1:l-0.2,col=statecolor);

    palette(defpalette);
    if(length(fout)!=0)
    {   dev.off();  }
    return(statecolor);
}
