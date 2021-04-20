
apt_quality_control <- function(exec, cellist, analysispath, xmlconfig, 
    outpath, verbose = FALSE) {
    run <- paste0(exec, 
        " --analysis-files-path ", analysispath,  
        " --xml-file ", xmlconfig, 
        " --cel-files ", cellist, 
        " --out-dir ", outpath, 
        " --out-file ", file.path(outpath, "apt-geno-qc.txt ", 
          fsep=.Platform$file.sep ), 
        " --log-file ", file.path(outpath, "apt-geno-qc.log", 
          fsep=.Platform$file.sep ))

    if (verbose) { 
      message("[ apt_quality_control ]: Running apt-geno-qc with:\n", run, "\n") 
    }
    system(run)

    return(file.path(outpath, "apt-geno-qc.txt", fsep=.Platform$file.sep))
}


apt_genotype_ax <- function(exec, cellist, analysispath, xmlconfig, 
    outpath, verbose = FALSE) {
    run <- paste0(exec, 
        " --log-file ", file.path(outpath, "apt-probeset-genotype.log", 
          fsep=.Platform$file.sep), 
        " --xml-file ", file.path(analysispath, xmlconfig, 
          fsep=.Platform$file.sep), 
        " --out-dir ", outpath, 
        " --cel-files ", cellist, 
        " --analysis-files-path ", analysispath, 
        " --summaries --write-models")
    
    if (verbose) { 
      message("[ apt_genotype ]: Running apt-probeset-genotype:\n", run, "\n")
    }
    
    system(run)
}

# metrix.fun
# ------------------------------------------------------------------------------
# From Affymetrix. 
# Function used to perform a metrix file to filter the SNPs in base different
# metrix: Call rate, FLD, HetSO and/or HomRO.

metrics.fun <- function( posteriors.file, calls.file, metrics.out.file, GTC ) {

    #Size of chunks of calls file to read and process at a time.  Consider reducing if calls file is very large or system has low RAM.  Do not set lower than 10,000.  Default set at 100,000 for Linux users. Windows users may need to reduce to 10,000.
    size<-100000

    processed<-"probeset_id"
    probesetid.f<-NULL
    CR.f<-NULL
    FLD.f<-NULL
    HomFLD.f<-NULL
    HetSO.f<-NULL
    HomRO.f<-NULL
    nMinorAllele.f<-NULL
    Nclus.f<-NULL
    nAA.f<-NULL
    nAB.f<-NULL
    nBB.f<-NULL
    nNC.f<-NULL

    #read input
    posteriors<-read.table(posteriors.file,header=T,stringsAsFactors=F)
    calls.header<-as.character(read.table(calls.file,stringsAsFactors=F,sep="\t",nrows=1)[1,])
    if (GTC==1) calls.header[1]<-"probeset_id"

    for (i in 0:99) {
        obj<-try(calls<-read.table(calls.file,stringsAsFactors=F,sep="\t",nrows=size,skip=i*size),silent=T)
        if (is(obj, "try-error")) break
        if (GTC==1 & i==0) calls[1,1]<-"probeset_id"
        names(calls)<-calls.header
        calls.1<-calls[!(calls[,1] %in% processed),]
        rm(calls)
        processed<-append(processed,calls.1[,1])
        if (nrow(calls.1)==0) break

        #clean up input
        posteriors.1<-posteriors[posteriors$id %in% calls.1$probeset_id,]
        calls.1<-calls.1[calls.1$probeset_id %in% posteriors$id,]
        posteriors.1<-posteriors.1[order(posteriors.1$id,posteriors.1$BB),]
        if (GTC==1) {
            calls.1<-calls.1[,c(1,grep("Call Codes",names(calls.1)))]
            if (is.null(dim(calls.1))) {
                print ("No A/B output")
                break
            }
            calls.1<-replace(calls.1, calls.1=="AA", 0)
            calls.1<-replace(calls.1, calls.1=="AB", 1)
            calls.1<-replace(calls.1, calls.1=="BB", 2)
            calls.1<-replace(calls.1, calls.1=="NoCall", -1)
        }
        calls.1<-calls.1[order(calls.1$probeset_id,-as.numeric(calls.1[,2])),]
        posteriors.1<-posteriors.1[!duplicated(posteriors.1$id),]
        calls.1<-calls.1[!duplicated(calls.1$probeset_id),]
        if (nrow(calls.1)==0) break

        #parse posterior values
        BB<-unlist(strsplit(posteriors.1$BB,","))
        BB.meanX<-BB[seq(1,length(BB),7)]
        BB.varX<-BB[seq(2,length(BB),7)]
        BB.nObsMean<-BB[seq(3,length(BB),7)]
        BB.nObsVar<-BB[seq(4,length(BB),7)]
        BB.meanY<-BB[seq(5,length(BB),7)]
        BB.varY<-BB[seq(6,length(BB),7)]
        BB.covarXY<-BB[seq(7,length(BB),7)]
        AB<-unlist(strsplit(posteriors.1$AB,","))
        AB.meanX<-AB[seq(1,length(AB),7)]
        AB.varX<-AB[seq(2,length(AB),7)]
        AB.nObsMean<-AB[seq(3,length(AB),7)]
        AB.nObsVar<-AB[seq(4,length(AB),7)]
        AB.meanY<-AB[seq(5,length(AB),7)]
        AB.varY<-AB[seq(6,length(AB),7)]
        AB.covarXY<-AB[seq(7,length(AB),7)]
        AA<-unlist(strsplit(posteriors.1$AA,","))
        AA.meanX<-AA[seq(1,length(AA),7)]
        AA.varX<-AA[seq(2,length(AA),7)]
        AA.nObsMean<-AA[seq(3,length(AA),7)]
        AA.nObsVar<-AA[seq(4,length(AA),7)]
        AA.meanY<-AA[seq(5,length(AA),7)]
        AA.varY<-AA[seq(6,length(AA),7)]
        AA.covarXY<-AA[seq(7,length(AA),7)]
        posteriors.2<-data.frame(id=posteriors.1$id,BB.meanX,BB.varX,BB.nObsMean,BB.nObsVar,BB.meanY,BB.varY,BB.covarXY,AB.meanX,AB.varX,AB.nObsMean,AB.nObsVar,AB.meanY,AB.varY,AB.covarXY,AA.meanX,AA.varX,AA.nObsMean,AA.nObsVar,AA.meanY,AA.varY,AA.covarXY,stringsAsFactors=F)

        #calculate potential fld values
        fld.1<-(as.numeric(posteriors.2$AA.meanX)-as.numeric(posteriors.2$AB.meanX))/sqrt(as.numeric(posteriors.2$AB.varX))
        fld.2<-(as.numeric(posteriors.2$AB.meanX)-as.numeric(posteriors.2$BB.meanX))/sqrt(as.numeric(posteriors.2$AB.varX))
        posteriors.3<-data.frame(id=posteriors.2$id,"fld.AA.AB"=fld.1,"fld.AB.BB"=fld.2,stringsAsFactors=F)

        #calculate number of calls
        calls.c0<-apply(calls.1[,2:ncol(calls.1)],1,function(x){length(x[x==0])})
        calls.c1<-apply(calls.1[,2:ncol(calls.1)],1,function(x){length(x[x==1])})
        calls.c2<-apply(calls.1[,2:ncol(calls.1)],1,function(x){length(x[x==2])})
        calls.c3<-apply(calls.1[,2:ncol(calls.1)],1,function(x){length(x[x==-1])})
        calls.2<-data.frame(probeset_id=calls.1$probeset_id,"AA"=calls.c0,"AB"=calls.c1,"BB"=calls.c2,"NN"=calls.c3,stringsAsFactors=F)
        rm(calls.1)
        gc()

        #determine which fld value to assign based on configuration of populated clusters
        findfld<-function(x){
            if(x[1]>0 & x[2]>0 & x[3]>0) return(4)
            if(x[1]>0 & x[2]>0 & x[3]==0) return(1)
            if(x[1]>0 & x[2]==0 & x[3]>0) return(0)
            if(x[1]==0 & x[2]>0 & x[3]>0) return(2)
            if(x[1]>0 & x[2]==0 & x[3]==0) return(0)
            if(x[1]==0 & x[2]>0 & x[3]==0) return(0)
            if(x[1]==0 & x[2]==0 & x[3]>0) return(0)
            if(x[1]==0 & x[2]==0 & x[3]==0) return(0)
        }
        whichfld<-apply(calls.2[,2:5],1,findfld)
        calls.2<-cbind(calls.2,whichfld,stringsAsFactors=F)
        posteriors.3<-cbind(posteriors.3,fld=-1,stringsAsFactors=F)
        if(nrow(calls.2[calls.2$whichfld==4,])>0)posteriors.3[posteriors.3$id %in% calls.2[calls.2$whichfld==4,1],]$fld<-apply(posteriors.3[posteriors.3$id %in% calls.2[calls.2$whichfld==4,1],2:3],1,min)
        if(nrow(calls.2[calls.2$whichfld==1,])>0)posteriors.3[posteriors.3$id %in% calls.2[calls.2$whichfld==1,1],]$fld<-posteriors.3[posteriors.3$id %in% calls.2[calls.2$whichfld==1,1],2]
        if(nrow(calls.2[calls.2$whichfld==2,])>0)posteriors.3[posteriors.3$id %in% calls.2[calls.2$whichfld==2,1],]$fld<-posteriors.3[posteriors.3$id %in% calls.2[calls.2$whichfld==2,1],3]

        #calculate homRO
        homRO<-data.frame(id=posteriors.3$id,"homROa"=rep(100,nrow(posteriors.3)),"homROb"=rep(-100,nrow(posteriors.3)),"homRO"=rep(0,nrow(posteriors.3)),stringsAsFactors=F)
        homRO[homRO$id %in% calls.2[calls.2$BB>0,1],"homROb"]<-posteriors.2[posteriors.2$id %in% calls.2[calls.2$BB>0,1],"BB.meanX"]
        homRO[homRO$id %in% calls.2[calls.2$AA>0,1],"homROa"]<-posteriors.2[posteriors.2$id %in% calls.2[calls.2$AA>0,1],"AA.meanX"]
        homRO$homRO<-apply(homRO[,2:3],1,function(x){min(abs(as.numeric(x[1])),abs(as.numeric(x[2])))})
        homRO[as.numeric(homRO$homROa)<0,"homRO"]<-homRO[as.numeric(homRO$homROa)<0,"homROa"]
        homRO[as.numeric(homRO$homROb)>0,"homRO"]<--1*as.numeric(homRO[as.numeric(homRO$homROb)>0,"homROb"])
        homRO[as.numeric(homRO$homRO)==100,"homRO"]<--10

        #calculate number of clusters
        findnc<-function(x){
            if(x[1]>0 & x[2]>0 & x[3]>0) return(3)
            if(x[1]>0 & x[2]>0 & x[3]==0) return(2)
            if(x[1]>0 & x[2]==0 & x[3]>0) return(2)
            if(x[1]==0 & x[2]>0 & x[3]>0) return(2)
            if(x[1]>0 & x[2]==0 & x[3]==0) return(1)
            if(x[1]==0 & x[2]>0 & x[3]==0) return(1)
            if(x[1]==0 & x[2]==0 & x[3]>0) return(1)
            if(x[1]==0 & x[2]==0 & x[3]==0) return(0)
        }
        nc<-apply(calls.2[,2:5],1,findnc)
        calls.2<-cbind(calls.2,nc,stringsAsFactors=F)

        #calculate hetSO
        hetSO<-data.frame(id=posteriors.3$id,"ye"=rep(0,nrow(posteriors.3)),"hetSO"=rep(0,nrow(posteriors.3)),stringsAsFactors=F)
        ye<-apply(posteriors.2,1,function(x){as.numeric(x[6])+(as.numeric(x[20])-as.numeric(x[6]))*((as.numeric(x[9])-as.numeric(x[2]))/(as.numeric(x[16])-as.numeric(x[2])))})
        hetSO$ye<-ye
        ye.1<-as.numeric(posteriors.2$AB.meanY)-ye
        hetSO$hetSO<-ye.1
        hetSO[hetSO$id %in% calls.2[calls.2$AB==0,1],3]<--10

        #calculate CR
        cr<-(calls.2$AA+calls.2$AB+calls.2$BB)/(calls.2$AA+calls.2$AB+calls.2$BB+calls.2$NN)*100

        #calculate HomFLD
        HomFLD<-data.frame(id=posteriors.3$id,"HomFLD"=rep(-1,nrow(posteriors.3)),stringsAsFactors=F)
        HomFLD$HomFLD<-(as.numeric(posteriors.2$AA.meanX)-as.numeric(posteriors.2$BB.meanX))/sqrt(as.numeric(posteriors.2$AB.varX))
        HomFLD[HomFLD$id %in% calls.2[calls.2$AA==0 | calls.2$BB==0,1],2]<--1

        #calculate nMinorAllele
        countMinorAllele<-function(x){
            if(which.max(x)==1) return(x[2]+x[3])
            if(which.max(x)==2) return(x[1]+x[3])
            if(which.max(x)==3) return(x[1]+x[2])
        }
        nMinorAllele<-apply(calls.2[,2:4],1,countMinorAllele)

        probesetid.f<-append(probesetid.f,posteriors.3$id)
        CR.f<-append(CR.f,cr)
        FLD.f<-append(FLD.f,posteriors.3$fld)
        HomFLD.f<-append(HomFLD.f,HomFLD$HomFLD)
        HetSO.f<-append(HetSO.f,hetSO$hetSO)
        HomRO.f<-append(HomRO.f,homRO$homRO)
        nMinorAllele.f<-append(nMinorAllele.f,nMinorAllele)
        Nclus.f<-append(Nclus.f,calls.2$nc)
        nAA.f<-append(nAA.f,calls.2$AA)
        nAB.f<-append(nAB.f,calls.2$AB)
        nBB.f<-append(nBB.f,calls.2$BB)
        nNC.f<-append(nNC.f,calls.2$NN)

    } #for-loop end

    #write out metric table
    write.table(data.frame("probeset_id"=probesetid.f,CR=CR.f,"FLD"=FLD.f,HomFLD=HomFLD.f,"HetSO"=HetSO.f,"HomRO"=HomRO.f,"nMinorAllele"=nMinorAllele.f,"Nclus"=Nclus.f,"n_AA"=nAA.f,"n_AB"=nAB.f,"n_BB"=nBB.f,"n_NC"=nNC.f),metrics.out.file,quote=F,sep="\t",row.names=F)
}