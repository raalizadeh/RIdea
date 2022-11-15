################################################################
#### This app is written by Dr. Reza Aalizadeh (C) 2022 ########
#### For any questions, drop email in raalizadeh@chem.uoa.gr ###
################################################################
jscode <- "shinyjs.closeWindow = function() { window.close(); }"
limitUse1<-5000 # batch mode limit

modelpH <-mx.model.load("www/models/pH_strict_30/model_final_pH", 662)
modelpredRI <-mx.model.load("www/models/model_reg_RI/model_DPL_RI_reg12", 1500)
modelpredRIhpH <-mx.model.load("www/models/model_reg_RI_negPH/model_DPL_RIneg", 1500)
modelChrom<-mx.model.load("www/models/LC/Chrom3d", 250)

FPmatneg<-readRDS("www/FPmatneg.rds")
FPmat<-readRDS("www/FPmat.rds")
dataX<-readRDS("www/dataX.rds")
databaseRI<-readRDS("www/databaseRI.rds")

shinyServer(function(input, output, session) {
  
  RTINDEA_cal<-function(Suspect,DataCal){
    DataCal<-DataCal[which(!is.na(DataCal$tR)),]
    IDn<-which(DataCal$tR<Suspect)
    IDn<-IDn[length(IDn)]
    
    if(rapportools::is.empty(IDn)){
      nDRT<-0
      nDRTDEA<-0}else{
        nDRT<-DataCal$tR[IDn]
        nDRTDEA<-DataCal$N[IDn]
      }
    IDnplus1<-which(DataCal$tR>=Suspect)[1]
    if(rapportools::is.empty(IDnplus1)){
      nPlus1DRT<-NA
      nPlus1DRTDEA<-NA
    }else{
      nPlus1DRT<-DataCal$tR[IDnplus1]
      nPlus1DRTDEA<-DataCal$N[IDnplus1]
    }

    
    RTI=(nDRTDEA+((Suspect-nDRT)/(nPlus1DRT-nDRT)))*100
    if(rapportools::is.empty(RTI)){RTI=NA}else{RTI<-round(RTI,2)}
    
    return((RTI))} 
  
  docanoSingleRI<-function(Compopunds){
    opsMW2<-function(SMILES){
      m <- parse.smiles(SMILES,kekulise=TRUE)[[1]]
      do.aromaticity(m) 
      canSMI<-get.smiles(m, flavor = smiles.flavors(c('Canonical')), smigen = NULL)
      return(as.character(canSMI))
    }
    allmv<-data.frame("SMILES"=as.array("NA"),"err"=as.array("NA"),"nm"=1:length(Compopunds))
    #allmv$MW<-0
    for (i in 1:length(Compopunds)) {
      result = tryCatch({
        result=list(result=opsMW2(as.character(Compopunds[i])),result2="Processed")
        
      }, error = function(e) {
        result=list(result=as.character(Compopunds[i]),result2="Failed")
        # allmv$err[i]<-"Failed"
      })
      allmv$SMILES[i]<-result$result
      allmv$err[i]<-result$result2
    }
    
    #  return(list(SMILES_given=Compopunds,SMILES_after=allmv$SMILES,statuss=allmv$err))
    return(allmv$SMILES[1])
  }
  # Get SMILES
  toSmiles <- function(mol){
    Smiles  <- J("org.openscience.cdk.smiles.SmilesGenerator")
    #Canonical
    smigen  <- new(Smiles)$unique()
    smigen$create(mol)
    return(smigen$create(mol))
  }
  # Get InChI
  toInChi <- function(mol){
    result = tryCatch({
      factory <- .jnew("org.openscience.cdk.inchi.InChIGeneratorFactory")$getInstance()
      result     <- factory$getInChIGenerator(mol)
      
    }, error = function(e) {
      result<-list(getInchiKey=NULL,stdinchi=NULL)
    })
    return(list(c(inchikey=result$getInchiKey(),stdinchi=result$getInchi())))
  }
  # Get InChIKey
  toInChikey <- function(mol){
    result = tryCatch({
      factory <- .jnew("org.openscience.cdk.inchi.InChIGeneratorFactory")$getInstance()
      result <- factory$getInChIGenerator(mol)
    }, warning = function(w) {
      result<-list(getInchiKey="",stdinchi="")
    },  error = function(e) {
      result<-list(getInchiKey="",stdinchi="")
    })
    return(result$getInchiKey())
  }
  
  SmilesToMW<-function(SMILES){
    out<-  tryCatch({
      mols <- parse.smiles(SMILES)
      do.typing(mols[[1]])
      convert.implicit.to.explicit(mols[[1]])
      do.aromaticity(mols[[1]])

      Des<-eval.desc(mols[[1]], which.desc=get.desc.names(type = "constitutional"), verbose = FALSE)
      
      out<-Des$MW
    },error=function(e){out<-NA})
    return(out)
  }
  
  FPFundes<-function(SMILES){
    pathRowrite<-paste(paste("www","HashFP","temp.mol",sep = "/"))
    FGfun2<-tryCatch({
      mols <- parse.smiles(SMILES)
      do.typing(mols[[1]])
      convert.implicit.to.explicit(mols[[1]])
      do.aromaticity(mols[[1]])
      
      mol2D<-generate.2d.coordinates(mols[[1]])
      write.molecules(mol2D,pathRowrite)
        comands<-paste("www/checkmol -c ",pathRowrite,sep="")
        gnps<-system(comands,intern = TRUE)
      unlink("www/temp.mol",recursive = TRUE)
      ######################################
      ### FP
      ######################################
      
      
      FGfun2<-data.frame("000000T2"=0, "000000T1"=0, "C2O10000"=0,"C2O1H000"=0,
                         "C2O1C000"=0,"C2S10000"=0,"C2S1H000"=0,"C2S1C000"=0,
                         "C2N10000"=0,"C2N1N000"=0,"C2NNC4ON"=0,"C2NNC4SN"=0,
                         "C2N1OH00"=0,"C2N1OC00"=0,"C3OC0000"=0,"C3OCC000"=0,
                         "C2O2H200"=0,"C2O2HC00"=0,"C2O2CC00"=0,"C2NOHC10"=0,
                         "C2N2CC10"=0,"C2NSHC10"=0,"C2S2CC00"=0,"C2CNH000"=0,
                         "C2COH000"=0,"C2COC000"=0,"O1H00000"=0,"O1H0C000"=0,
                         "O1H1C000"=0,"O1H2C000"=0,"O1H3C000"=0,"O1H0CO1H"=0,
                         "O1H0CN1C"=0,"O1H1A000"=0,"O1H2A000"=0,"C2COH200"=0,
                         "O1C00000"=0,"O1C0CC00"=0,"O1C0CA00"=0,"O1C0AA00"=0,
                         "S1C00000"=0,"S1S1C000"=0,"O1O1C000"=0,"O1O1H000"=0,
                         "N1N10000"=0,"N1O1H000"=0,"N1C00000"=0,"N1C10000"=0,
                         "N1C1C000"=0,"N1C1A000"=0,"N1C20000"=0,"N1C2CC00"=0,
                         "N1C2AC00"=0,"N1C2AA00"=0,"N1C30000"=0,"N1C3CC00"=0,
                         "N1C3AC00"=0,"N1C3AA00"=0,"N1C400T2"=0,"N0O10000"=0,
                         "XX000000"=0,"XX00C000"=0,"XF00C000"=0,"XC00C000"=0,
                         "XB00C000"=0,"XI00C000"=0,"XX00A000"=0,"XF00A000"=0,
                         "XC00A000"=0,"XB00A000"=0,"XI00A000"=0,"000000MX"=0,
                         "000000ML"=0,"000000MM"=0,"C3O20000"=0,"C3O2H000"=0,
                         "C3O200T1"=0,"C3O2C000"=0,"C3O2CZ00"=0,"C3ONC000"=0,
                         "C3ONC100"=0,"C3ONC200"=0,"C3ONC300"=0,"C3ONCZ00"=0,
                         "C3ONN100"=0,"C3ONN200"=0,"C3ONOH00"=0,"C3N2H000"=0,
                         "C3NNN100"=0,"C3N00000"=0,"C3OXX000"=0,"C3OXF000"=0,
                         "C3OXC000"=0,"C3OXB000"=0,"C3OXI000"=0,"C2OC3N00"=0,
                         "C3NOC000"=0,"C3NXX000"=0,"C3SO0000"=0,"C3SOH000"=0,
                         "C3SOC000"=0,"C3SOCZ00"=0,"C3SNH000"=0,"C3SNCZ00"=0,
                         "C3NSC000"=0,"C3ONAZ00"=0,"C3SNAZ00"=0,"C3NNAZ00"=0,
                         "C3O30000"=0,"C3O3C000"=0,"C3O3NC00"=0,"C3O2C3O2"=0,
                         "C3ONC000"=0,"C3ONCH10"=0,"C3ONCC10"=0,"C4000000"=0,
                         "C4O30000"=0,"C4O3C100"=0,"C4O3C200"=0,"C4O3CX00"=0,
                         "C4SO0000"=0,"C4SOC100"=0,"C4SOC200"=0,"C4SOX_00"=0,
                         "C4O2N000"=0,"C4O2NH00"=0,"C4O2NC00"=0,"C4O2NX00"=0,
                         "C4SN0000"=0,"C4SNOH00"=0,"C4SNOC00"=0,"C4SNXX00"=0,
                         "C4O1N200"=0,"C4N2O100"=0,"C4S1N200"=0,"C4N2S100"=0,
                         "C4N30000"=0,"C4ON2N00"=0,"C4SN2N00"=0,"N4N20000"=0,
                         "N2N10000"=0,"N3N100T2"=0,"N3C10000"=0,"C4NO1000"=0,
                         "C4NO2000"=0,"C4NS1000"=0,"C4NS2000"=0,"C4N20000"=0,
                         "N2O10000"=0,"N4O20000"=0,"N3O20000"=0,"N4O30000"=0,
                         "S6O00000"=0,"S6O4H000"=0,"S6O4HC00"=0,"S6O4CC00"=0,
                         "S6O3NC00"=0,"S6O3N100"=0,"S6O2N200"=0,"S6O3XX00"=0,
                         "S5O00000"=0,"S5O3H000"=0,"S5O3C000"=0,"S5O2N000"=0,
                         "S5O2XX00"=0,"S4O20000"=0,"S2O10000"=0,"S3O00000"=0,
                         "S3O2H000"=0,"S3O2C000"=0,"S3O1XX00"=0,"S3O1N000"=0,
                         "S1O00000"=0,"S1O1H000"=0,"S1O1C000"=0,"S1O0XX00"=0,
                         "S1O0N100"=0,"S1H10000"=0,"S1H1C000"=0,"S1H1A000"=0,
                         "P5O0H000"=0,"P5O4H200"=0,"P5O4HC00"=0,"P5O3HX00"=0,
                         "P5O3HN00"=0,"P5O0S000"=0,"P5O3SH00"=0,"P5O3SC00"=0,
                         "P5O2SX00"=0,"P5O2SN00"=0,"P4O30000"=0,"P4O3H000"=0,
                         "P4O3C000"=0,"P3000000"=0,"P2O00000"=0,"B2O20000"=0,
                         "B2O2H000"=0,"B2O2C000"=0,"000C2C00"=0,"000C3C00"=0,
                         "0000A000"=0,"0000CZ00"=0,"C3O2HN1C"=0,"C3O2HO1H"=0,
                         check.names = FALSE)
      
      foundS<-strsplit(gnps,";")
      ID<-which(colnames(FGfun2)%in%foundS[[1]])
      FGfun2[ID]<-1
      FGfun2
    },error=function(e){
      FGfun2<-data.frame("000000T2"=NA,"000000T1"=NA, "C2O10000"=NA,"C2O1H000"=NA,
                         "C2O1C000"=NA,"C2S10000"=NA,"C2S1H000"=NA,"C2S1C000"=NA,
                         "C2N10000"=NA,"C2N1N000"=NA,"C2NNC4ON"=NA,"C2NNC4SN"=NA,
                         "C2N1OH00"=NA,"C2N1OC00"=NA,"C3OC0000"=NA,"C3OCC000"=NA,
                         "C2O2H200"=NA,"C2O2HC00"=NA,"C2O2CC00"=NA,"C2NOHC10"=NA,
                         "C2N2CC10"=NA,"C2NSHC10"=NA,"C2S2CC00"=NA,"C2CNH000"=NA,
                         "C2COH000"=NA,"C2COC000"=NA,"O1H00000"=NA,"O1H0C000"=NA,
                         "O1H1C000"=NA,"O1H2C000"=NA,"O1H3C000"=NA,"O1H0CO1H"=NA,
                         "O1H0CN1C"=NA,"O1H1A000"=NA,"O1H2A000"=NA,"C2COH200"=NA,
                         "O1C00000"=NA,"O1C0CC00"=NA,"O1C0CA00"=NA,"O1C0AA00"=NA,
                         "S1C00000"=NA,"S1S1C000"=NA,"O1O1C000"=NA,"O1O1H000"=NA,
                         "N1N10000"=NA,"N1O1H000"=NA,"N1C00000"=NA,"N1C10000"=NA,
                         "N1C1C000"=NA,"N1C1A000"=NA,"N1C20000"=NA,"N1C2CC00"=NA,
                         "N1C2AC00"=NA,"N1C2AA00"=NA,"N1C30000"=NA,"N1C3CC00"=NA,
                         "N1C3AC00"=NA,"N1C3AA00"=NA,"N1C400T2"=NA,"N0O10000"=NA,
                         "XX000000"=NA,"XX00C000"=NA,"XF00C000"=NA,"XC00C000"=NA,
                         "XB00C000"=NA,"XI00C000"=NA,"XX00A000"=NA,"XF00A000"=NA,
                         "XC00A000"=NA,"XB00A000"=NA,"XI00A000"=NA,"000000MX"=NA,
                         "000000ML"=NA,"000000MM"=NA,"C3O20000"=NA,"C3O2H000"=NA,
                         "C3O200T1"=NA,"C3O2C000"=NA,"C3O2CZ00"=NA,"C3ONC000"=NA,
                         "C3ONC100"=NA,"C3ONC200"=NA,"C3ONC300"=NA,"C3ONCZ00"=NA,
                         "C3ONN100"=NA,"C3ONN200"=NA,"C3ONOH00"=NA,"C3N2H000"=NA,
                         "C3NNN100"=NA,"C3N00000"=NA,"C3OXX000"=NA,"C3OXF000"=NA,
                         "C3OXC000"=NA,"C3OXB000"=NA,"C3OXI000"=NA,"C2OC3N00"=NA,
                         "C3NOC000"=NA,"C3NXX000"=NA,"C3SO0000"=NA,"C3SOH000"=NA,
                         "C3SOC000"=NA,"C3SOCZ00"=NA,"C3SNH000"=NA,"C3SNCZ00"=NA,
                         "C3NSC000"=NA,"C3ONAZ00"=NA,"C3SNAZ00"=NA,"C3NNAZ00"=NA,
                         "C3O30000"=NA,"C3O3C000"=NA,"C3O3NC00"=NA,"C3O2C3O2"=NA,
                         "C3ONC000"=NA,"C3ONCH10"=NA,"C3ONCC10"=NA,"C4000000"=NA,
                         "C4O30000"=NA,"C4O3C100"=NA,"C4O3C200"=NA,"C4O3CX00"=NA,
                         "C4SO0000"=NA,"C4SOC100"=NA,"C4SOC200"=NA,"C4SOX_00"=NA,
                         "C4O2N000"=NA,"C4O2NH00"=NA,"C4O2NC00"=NA,"C4O2NX00"=NA,
                         "C4SN0000"=NA,"C4SNOH00"=NA,"C4SNOC00"=NA,"C4SNXX00"=NA,
                         "C4O1N200"=NA,"C4N2O100"=NA,"C4S1N200"=NA,"C4N2S100"=NA,
                         "C4N30000"=NA,"C4ON2N00"=NA,"C4SN2N00"=NA,"N4N20000"=NA,
                         "N2N10000"=NA,"N3N100T2"=NA,"N3C10000"=NA,"C4NO1000"=NA,
                         "C4NO2000"=NA,"C4NS1000"=NA,"C4NS2000"=NA,"C4N20000"=NA,
                         "N2O10000"=NA,"N4O20000"=NA,"N3O20000"=NA,"N4O30000"=NA,
                         "S6O00000"=NA,"S6O4H000"=NA,"S6O4HC00"=NA,"S6O4CC00"=NA,
                         "S6O3NC00"=NA,"S6O3N100"=NA,"S6O2N200"=NA,"S6O3XX00"=NA,
                         "S5O00000"=NA,"S5O3H000"=NA,"S5O3C000"=NA,"S5O2N000"=NA,
                         "S5O2XX00"=NA,"S4O20000"=NA,"S2O10000"=NA,"S3O00000"=NA,
                         "S3O2H000"=NA,"S3O2C000"=NA,"S3O1XX00"=NA,"S3O1N000"=NA,
                         "S1O00000"=NA,"S1O1H000"=NA,"S1O1C000"=NA,"S1O0XX00"=NA,
                         "S1O0N100"=NA,"S1H10000"=NA,"S1H1C000"=NA,"S1H1A000"=NA,
                         "P5O0H000"=NA,"P5O4H200"=NA,"P5O4HC00"=NA,"P5O3HX00"=NA,
                         "P5O3HN00"=NA,"P5O0S000"=NA,"P5O3SH00"=NA,"P5O3SC00"=NA,
                         "P5O2SX00"=NA,"P5O2SN00"=NA,"P4O30000"=NA,"P4O3H000"=NA,
                         "P4O3C000"=NA,"P3000000"=NA,"P2O00000"=NA,"B2O20000"=NA,
                         "B2O2H000"=NA,"B2O2C000"=NA,"000C2C00"=NA,"000C3C00"=NA,
                         "0000A000"=NA,"0000CZ00"=NA,"C3O2HN1C"=NA,"C3O2HO1H"=NA,
                         check.names = FALSE)
      
    })
    
    return(FGfun2)
  }
  
  FPQdes<-function(SMILES){
    mols<-tryCatch({
      mols <- parse.smiles(SMILES)
      do.typing(mols[[1]])
      convert.implicit.to.explicit(mols[[1]])
      do.aromaticity(mols[[1]])
      mols
    },error=function(e){mols=NULL})
    ######################################
    ### FP
    ######################################
    titID2<- tryCatch({
      titID<- paste("FP",seq(1,1024,1), sep="")
      fp_<- get.fingerprint(mols[[1]], type='standard',fp.mode = "bit",depth=6, size=1024)
      titID2<-as.data.frame(matrix(0,1,1024))
      titID2[fp_@bits]<-rep(1,length(fp_@bits))
      colnames(titID2)<-titID
      titID2
    },error=function(e){
      titID<- paste("FP",seq(1,1024,1), sep="")
      titID2<-as.data.frame(matrix(NA,1,1024))
      colnames(titID2)<-titID
      titID2
    }) 
    FP<-titID2
    ######################################
    ### ExtendFP
    ######################################
    titID2<- tryCatch({
      titID<- paste("ExtFP",seq(1,1024,1), sep="")
      fp_extend<- get.fingerprint(mols[[1]], type='extended',fp.mode = "bit",depth=6, size=1024)
      titID2<-as.data.frame(matrix(0,1,1024))
      titID2[fp_extend@bits]<-rep(1,length(fp_extend@bits))
      colnames(titID2)<-titID
      titID2
    },error=function(e){
      titID<- paste("ExtFP",seq(1,1024,1), sep="")
      titID2<-as.data.frame(matrix(NA,1,1024))
      colnames(titID2)<-titID
      titID2
    }) 
    ExtendFP<-titID2
    ######################################
    ### graphFP
    ######################################
    titID2<- tryCatch({
      titID<- paste("GraphFP",seq(1,1024,1), sep="")
      fp_GraphFP<- get.fingerprint(mols[[1]], type='graph',fp.mode = "bit",depth=6, size=1024)
      titID2<-as.data.frame(matrix(0,1,1024))
      titID2[fp_GraphFP@bits]<-rep(1,length(fp_GraphFP@bits))
      colnames(titID2)<-titID
      titID2
    },error=function(e){
      titID<- paste("GraphFP",seq(1,1024,1), sep="")
      titID2<-as.data.frame(matrix(NA,1,1024))
      colnames(titID2)<-titID
      titID2
    }) 
    graphFP<-titID2
    ######################################
    ### PubChem
    ######################################  
    titID2<- tryCatch({
      titID<- paste("PubchemFP",seq(0,880,1), sep="")
      fp_pubChem<- get.fingerprint(mols[[1]], type='pubchem',fp.mode = "bit")
      titID2<-as.data.frame(matrix(0,1,881))
      titID2[fp_pubChem@bits]<-rep(1,length(fp_pubChem@bits))
      colnames(titID2)<-titID
      titID2
    },error=function(e){
      titID<- paste("PubchemFP",seq(0,880,1), sep="")
      titID2<-as.data.frame(matrix(NA,1,881))
      colnames(titID2)<-titID
      titID2
    })
    PubChem<-titID2
    ######################################
    ### EStateFP
    ######################################  
    titID2<- tryCatch({
      titID<- paste("EStateFP",seq(1,79,1), sep="")
      fp_EStateFP<- get.fingerprint(mols[[1]], type='estate',fp.mode = "bit")
      titID2<-as.data.frame(matrix(0,1,79))
      titID2[fp_EStateFP@bits]<-rep(1,length(fp_EStateFP@bits))
      colnames(titID2)<-titID
      titID2
    },error=function(e){
      titID<- paste("EStateFP",seq(1,79,1), sep="")
      titID2<-as.data.frame(matrix(NA,1,79))
      colnames(titID2)<-titID
      titID2
    })
    EStateFP<-titID2
    ######################################
    ### maccs
    ######################################  
    titID2<- tryCatch({
      titID<- paste("MACCSFP",seq(1,166,1), sep="")
      fp_MACCSFP<- get.fingerprint(mols[[1]], type='maccs',fp.mode = "bit")
      titID2<-as.data.frame(matrix(0,1,166))
      titID2[fp_MACCSFP@bits]<-rep(1,length(fp_MACCSFP@bits))
      colnames(titID2)<-titID
      titID2
    },error=function(e){
      titID<- paste("MACCSFP",seq(1,166,1), sep="")
      titID2<-as.data.frame(matrix(NA,1,166))
      colnames(titID2)<-titID
      titID2
    })
    maccs<-titID2
    ######################################
    ### Klekota and Roth
    ######################################  
    titID2<- tryCatch({
      titID<- paste("KRFP",seq(1,4860,1), sep="")
      fp_KRFP<- get.fingerprint(mols[[1]], type='kr',fp.mode = "bit")
      titID2<-as.data.frame(matrix(0,1,4860))
      titID2[fp_KRFP@bits]<-rep(1,length(fp_KRFP@bits))
      colnames(titID2)<-titID
      titID2
    },error=function(e){
      titID<- paste("KRFP",seq(1,4860,1), sep="")
      titID2<-as.data.frame(matrix(NA,1,4860))
      colnames(titID2)<-titID
      titID2
    })
    KR<-titID2
    ######################################
    ### shortestpath
    ######################################  
    titID2<- tryCatch({
      titID<- paste("SPFP",seq(1,1024,1), sep="")
      fp_shortestpath<- get.fingerprint(mols[[1]], type='shortestpath',fp.mode = "bit")
      titID2<-as.data.frame(matrix(0,1,1024))
      titID2[fp_shortestpath@bits]<-rep(1,length(fp_shortestpath@bits))
      colnames(titID2)<-titID
      titID2
    },error=function(e){
      titID<- paste("SPFP",seq(1,1024,1), sep="")
      titID2<-as.data.frame(matrix(NA,1,1024))
      colnames(titID2)<-titID
      titID2
    })
    SHT<-titID2
    ######################################
    ### Combine data
    ######################################
    out<-cbind.data.frame(FP,ExtendFP,graphFP,PubChem,EStateFP,maccs,KR,SHT)
    return(out)
  }
  
  DoAdduct_formula<-function(mol){
    allmv<- tryCatch({
      allmv<-data.frame("Formula"=as.array(0),"Monoisotopic Mass"=as.array(0),
                        "[M]+"=as.array(0),"[M+H]+"=as.array(0),"[M+NH4]+"=as.array(0),
                        "[M+Na]+"=as.array(0),"[M+K]+"=as.array(0),"[M+CH3OH+H]+"=as.array(0),"[M]2+"=as.array(0),
                        "[M+H]2+"=as.array(0),"[2M+H]+"=as.array(0),"[3M+H]+"=as.array(0),"[M+2H]+"=as.array(0),"[M+3H]+"=as.array(0),
                        "[M+H+1]+"=as.array(0),"[M+H+2]+"=as.array(0),"[M+NH4]2+"=as.array(0),
                        "[M]-"=as.array(0),"[M-H]-"=as.array(0),"[M-H2O-H]-"=as.array(0),"[M+Cl]-"=as.array(0),
                        "[M+HCOOH-H]-"=as.array(0),"[M+CH3COOH-H]-"=as.array(0),"[2M-H]-"=as.array(0),
                        "[3M-H]-"=as.array(0),"[M]2-"=as.array(0),"[M-2H]-"=as.array(0),
                        "[M-3H]-"=as.array(0),"[M-H-1]-"=as.array(0),"[M-H-2]-"=as.array(0),check.names = FALSE)
      
      #  SMILES<-as.character(mol)
      #  molecule <- parse.smiles(SMILES)[[1]]
      #  convert.implicit.to.explicit(molecule)
      formula <- get.formula(mol,charge=+1)
      formula2 <- get.formula(mol,charge=+2)
      formula3 <- get.formula(mol,charge=+1)
      formula4 <- get.formula(mol,charge=0)
      allmv$Formula<-as.array(formula@string)
      mz<-formula@mass+1.007825032
      allmv$`[M+H]+`<-mz
      allmv$`[M+2H]+`<-mz+1.007825032
      allmv$`[M+3H]+`<-mz+1.007825032+1.007825032
      allmv$`[M+H+1]+`<-mz+1
      allmv$`[M+H+2]+`<-mz+2
      mzForm<-formula@string<-paste("[",formula@string,"]+H",sep = "")
      allmv$`Monoisotopic Mass`<-formula4@mass
      
      
      SMILES<-as.character("[H][N]([H])([H])")
      molecule <- parse.smiles(SMILES)[[1]]
      convert.implicit.to.explicit(molecule)
      formula <- get.mol2formula(molecule,charge=0)
      mzNH<-formula@mass # M+NH4
      allmv$`[M+NH4]+`<-mzNH+mz
      
      
      SMILES<-as.character("[Na]")
      molecule <- parse.smiles(SMILES)[[1]]
      convert.implicit.to.explicit(molecule)
      formula <- get.mol2formula(molecule,charge=0)
      mzNa<-formula@mass # M+Na
      allmv$`[M+Na]+`<-mzNa+formula3@mass
      
      SMILES<-as.character("[K]")
      molecule <- parse.smiles(SMILES)[[1]]
      convert.implicit.to.explicit(molecule)
      formula <- get.mol2formula(molecule,charge=0)
      mzK<-formula@mass # M+K
      allmv$`[M+K]+`<-mzK+formula3@mass
      
      SMILES<-as.character("[H]C([H])([H])O") 
      molecule <- parse.smiles(SMILES)[[1]]
      convert.implicit.to.explicit(molecule)
      formula <- get.mol2formula(molecule,charge=0)
      mzMeOH<-formula@mass
      MeOHaduct<-mz+mzMeOH # M+CH3OH+H
      allmv$`[M+CH3OH+H]+`<-MeOHaduct
      
      mzM<-formula3@mass # [M+]
      allmv$`[M]+`<-mzM
      
      mz2M<-formula2@mass/2 # [M2+]
      allmv$`[M]2+`<-mz2M
      
      
      mz2MH<-(formula2@mass/2)+1.007825032 # [M+H2+]
      allmv$`[M+H]2+`<-mz2MH
      allmv$`[M+NH4]2+`<-(mzNH+mz2MH)
      mz22mh<-(2*formula3@mass)+1.007825032+(formula4@mass-formula3@mass)# [2M+H]+
      allmv$`[2M+H]+`<-mz22mh
      mz32mh<-(3*formula3@mass)+1.007825032+2*(formula4@mass-formula3@mass)# [3M+H]+
      allmv$`[3M+H]+`<-mz32mh
      ################## Do similar task for neg
      charg<-formula4@mass-formula3@mass
      #  SMILES<-as.character((mol))
      #  molecule <- parse.smiles(SMILES)[[1]]
      # convert.implicit.to.explicit(molecule)
      formula <- get.formula(mol,charge=-1)
      formula2 <- get.formula(mol,charge=-2)
      formula3 <- get.formula(mol,charge=-1)
      formula4 <- get.formula(mol,charge=0)
      mz<-formula@mass-1.007825032
      allmv$`[M-2H]-`<-mz-1.007825032
      allmv$`[M-3H]-`<-mz-1.007825032-1.007825032
      allmv$`[M-H-1]-`<-mz-1
      allmv$`[M-H-2]-`<-mz-2
      mz_minus<-formula3@mass # [M-]
      allmv$`[M]-`<-mz_minus
      
      mz_H_minus<-mz # [M-H]-
      allmv$`[M-H]-`<-mz
      
      SMILES<-as.character("O[H][H]") 
      molecule <- parse.smiles(SMILES)[[1]]
      convert.implicit.to.explicit(molecule)
      formula <- get.mol2formula(molecule,charge=0)
      mzHOH<-formula@mass
      mz_adduct_HOH<-mz_minus-mzHOH # M-H2O-H
      allmv$`[M-H2O-H]-`<-mz_adduct_HOH
      
      SMILES<-as.character("[Cl]") 
      molecule <- parse.smiles(SMILES)[[1]]
      convert.implicit.to.explicit(molecule)
      formula <- get.mol2formula(molecule,charge=0)
      mzcl<-formula@mass
      mz_adduct_cl<-mz_minus+mzcl # M+cl
      allmv$`[M+Cl]-`<-mz_adduct_cl
      
      SMILES<-as.character("[H]C(O)=O") 
      molecule <- parse.smiles(SMILES)[[1]]
      convert.implicit.to.explicit(molecule)
      formula <- get.mol2formula(molecule,charge=0)
      mzformic<-formula@mass
      mz_adduct_formic<-mz_H_minus+mzformic # M+HCOOH-H
      allmv$`[M+HCOOH-H]-`<-mz_adduct_formic
      
      SMILES<-as.character("CC(O)=O") 
      molecule <- parse.smiles(SMILES)[[1]]
      convert.implicit.to.explicit(molecule)
      formula <- get.mol2formula(molecule,charge=0)
      mzAcetic<-formula@mass
      mz_adduct_Acetic<-mz_H_minus+mzAcetic # M+CH3COOH-H
      allmv$`[M+CH3COOH-H]-`<-mz_adduct_Acetic
      
      mz_2MH_minus<-(2*mz_minus)-1.007825032-charg # 2M-H
      allmv$`[2M-H]-`<- mz_2MH_minus
      
      mz_3MH_minus<-(3*mz_minus)-1.007825032-2*charg # 3M-H
      allmv$`[3M-H]-`<-mz_3MH_minus
      
      mz_M2_minus<-(formula2@mass)/2 # [M2-]
      allmv$`[M]2-`<-mz_M2_minus
      
      allmv
    },error=function(e){
      allmv<-data.frame("Formula"="Error","Monoisotopic Mass"="Error","[M]+"="Error",
                        "[M+H]+"="Error","[M+NH4]+"="Error",
                        "[M+Na]+"="Error","[M+K]+"="Error","[M+CH3OH+H]+"="Error","[M]2+"="Error",
                        "[M+H]2+"="Error","[2M+H]+"="Error","[3M+H]+"="Error","[M+2H]+"="Error",
                        "[M+3H]+"="Error","[M+H+1]+"="Error","[M+H+2]+"="Error","[M+NH4]2+"="Error",
                        "[M]-"="Error","[M-H]-"="Error","[M-H2O-H]-"="Error","[M+Cl]-"="Error",
                        "[M+HCOOH-H]-"="Error","[M+CH3COOH-H]-"="Error","[2M-H]-"="Error",
                        "[3M-H]-"="Error","[M]2-"="Error","[M-2H]-"="Error","[M-3H]-"="Error",
                        "[M-H-1]-"="Error","[M-H-2]-"="Error",check.names = FALSE)
      
    })
    
    return(allmv)
  }
  
  sourceCpp(file="www/Ex_6.cpp")
  CompareFP<-function(comp1,comp2){
    res<-Ex_6(comp1, comp2)
    commonXSus=length(res[[2]])
    all<-length(res[[1]])
    Score<-commonXSus/all
    return (Score)
  }
  CompareFPSmiles<-function(comp1,comp2){
    
    out1<-FPQdes(as.character(comp1))
    out11<-FPFundes(as.character(comp1))
    out01<-cbind(out1,out11)
    out2<-FPQdes(as.character(comp2))
    out22<-FPFundes(as.character(comp2))
    out02<-cbind(out2,out22)
    
    OneX<-which(out01==1)
    OneSus<-which(out02==1)
    commonXSus<-c(intersect(OneX, OneSus),intersect(OneSus,OneX))
    commonXSus<-commonXSus[!duplicated(commonXSus)]
    UncommonXSus<-c(setdiff(OneX,OneSus), setdiff(OneSus,OneX))
    UncommonXSus<-UncommonXSus[!duplicated(UncommonXSus)]
    
    Score<-length(commonXSus)/(length(commonXSus)+length(UncommonXSus)) 
    
    
    
    return (Score)
  }
  MedTP5 <- function(indf) {
    #do.call(rbind, with(indf, by(Value, Name, FUN=function(x) {
    Vals <- head(sort(indf, decreasing=TRUE), 10)
    if (length(Vals) < 3) {
      out<-c(Mean = NA, Median = NA)
    } else {
      out<-c(Mean = mean(Vals), Median = median(Vals))
    }
    #}
    
    return(out)
  }
  extract_feature <-function(dir_path, width, height, AD0 = TRUE,AD1 = TRUE,
                             add_label = TRUE, cols) {
    img_size <- width*height
    ## List images in path
    images_names <- list.files(dir_path)
    if (add_label) {
     
      if (AD0 == "TRUE"){
        label <- 0
        
      }else if (AD1 == "TRUE"){label <- 1}
      
    }
    
    
    images_na<- data.frame("images_na"='png',stringsAsFactors=TRUE )
    images_na[1:length(images_names),1]<-'png'
    
    print(paste("Start processing", length(images_names), "RI"))
    ## This function will resize an image, turn it into greyscale
    feature_list <- pblapply(images_names, function(imgname) {
      ## Read image
      img2 <- readImage(file.path(dir_path, imgname))
      img<-img2
      img_resized <- resize(img, w = width, h = height)  
      grayimg <- channel(img_resized, cols)
      # display(grayimg)
      ## Get the image as a matrix
      img_matrix <- grayimg@.Data
      ## Coerce to a vector
      img_vector <- as.vector(t(img_matrix))
      return(img_vector)
    })
    ## bind the list of vector into matrix
    feature_matrix <- do.call(rbind, feature_list)
    feature_matrix <- as.data.frame(feature_matrix)
    ## Set names
    names(feature_matrix) <- paste0("pixel", c(1:img_size))
    if (add_label) {
      ## Add label
      feature_matrix <- cbind(label = label, feature_matrix)
    }
    return(list(aa=feature_matrix,
                bb=images_names))
  }

  #### UI code --------------------------------------------------------------
  output$ui <- renderUI({
      ##### UI code for login page

      #### Your app's UI code goes here!

      fluidPage(tags$head(tags$style(type = "text/css", '.well{width: 365px}')),
        
        
            mainPanel(#tags$head(tags$style(type = "text/css", '.well{width: 365px}')),
              tabsetPanel(
                tabPanel("About",
                         h4("Citation:"),
                         h6("R. Aalizadeh, V. Nikolopoulou and N. S. Thomaidis,"),
                         h6("Development of Liquid Chromatographic Retention Index based on Cocamide Diethanolamine Homologous Series (C(n)-DEA), Analytical Chemistry, 2022. DOI:10.1021/acs.analchem.2c02893"),
                         downloadButton("ExportCit", "Export Citation",style="color: black; font-size: 14px; border-color: #2e6da4;width:150px"),
                         h6(""),
                         tags$head(
                         ),tags$img(src="TOC.png"),h6(""),
                         downloadButton("manRI", "Download Manual",style="color: black; font-size: 14px; border-color: #2e6da4;width:340px"),
                         
                         h6("___________________________________________________________________________________________________________________________________________________________"),
                         h4("Are you interested to purchase RI mix?"),
                         h6("Contact information:"),
                         h6("Dr. Reza Aalizadeh (raalizadeh@chem.uoa.gr)"),
                         h6("Dr. Varvara Nikolopoulou (vnikol@chem.uoa.gr)"),
                         h6("Prof. Nikolaos Thomaidis (ntho@chem.uoa.gr)"),
                         
                         h6(""),
                         tags$head(
                         ),tags$img(src="all_product1.png"),
                                      ),
                tabPanel("Prediction of RI",                 
                         h6(""),
                         
                         h6("___________________________________________________________________________________________________________________________________________________________"),
                         
                         radioButtons("method_pRI", "pH of mobile phase:",
                                      c("pH >= 6"="neg",
                                        "pH <= 4"="pos")),
                         h6(""),
                         
                         downloadButton("downloadRIfile", "Download CDEA RI calibrants Template",style="color: black; font-size: 14px; border-color: #2e6da4;width:340px"),
                         h6(""),
                         fileInput("file1RI", "Upload the tR data of CDEA calibrant...",multiple = FALSE,accept =".csv", width="490px"),
                         textAreaInput("caption", "Insert SMILES and tR here...", "", width = "800px",height="320px"),
                         textOutput("limitUseRI"),
                         textOutput("CDF1"),
                         textOutput("CDF2"),
                         textOutput("CDF3"),
                         textOutput("CDF4"),
                        
                         tags$head(tags$style("#limitUseRI{color: red;
                                                                   font-size: 14px;
                                                                   font-style: italic;
    }"
                         )
                         ),
                         
                         h6("___________________________________________________________________________________________________________________________________________________________"),
                         
                         h6(""),
                         actionButton("butt1CRI","Submit the calculation",
                                      style="color: black; font-size: 14px; border-color: #2e6da4;width:320px"),
                         h6(""),
                         DT::dataTableOutput("contents1CRI"),
                         downloadButton("downsusCRI", "Download the list",
                                        style="color: black; font-size: 14px;background-color: #aecde8; 
                                 border-color: #2e6da4;width:320px"),
                         
                         #h6("___________________________________________________________________________________________________________________________________________________________"),
                         h6("")
                ),
                
                tabPanel("RI bank",
                         textInput("MZSuspect","Enter the m/z value (Da) here:",""),
                         tags$head(tags$style("#MZSuspect{color: black;
                                                font-size: 14px;
                                                font-style: italic;
  }"
                         )
                         ),
                         textInput("MZacc","Enter mass accuracy (Da):",0.002),
                         tags$head(tags$style("#MZacc{color: black;
                              font-size: 14px;
                              font-style: italic;
  }"
                         )
                         ),
                         selectInput("IonTypeCRI", "Ion Type:",
                                     c("[M]+"="Mp","[M+H]+"="MpH","[M+NH4]+"="MpNH",
                                       "[M+Na]+"="MpNa","[M+K]+"="MpK","[M+CH3OH+H]+"="MpAct","[M]2+"="Mp2",
                                       "[M+H]2+"="MpH2p","[M]-"="Mm","[M-H]-"="MmH","[M-H2O-H]-"="MmH2Omh","[M+Cl]-"="MpClm",
                                       "[M+HCOOH-H]-"="MpFAmH","[M+CH3COOH-H]-"="MpActmH","[2M-H]-"="2MmH",
                                       "[M]2-"="M2m","[M-2H]-"="Mm2H")),
                         
                         actionButton("butt_libCRI","Search Database",
                                      style="color: black; font-size: 14px; border-color: #2e6da4;width:300px"),
                         h6(""),
                         textOutput("warn1"),
                         tags$head(tags$style("#warn1{color: red;
                              font-size: 14px;
                              font-style: bold;
  }"
                         )),
                         DT::dataTableOutput("contents_libCRI"),
                         downloadButton("downsus_libCRI", "Download the list",
                                        style="color: black; font-size: 14px; border-color: #2e6da4;width:300px")
                )
              )
            ), 
      h6("___________________________________________________________________________________________________________________________________________________________"),
      
tags$head(
),tags$img(src="logo.jpg"),
h6("The shiny app is desinged by Dr. Reza Aalizadeh"),
h6("Laboratory of Analytical Chemistry") ,
h6("National and Kapodistrian University of Athens"),
h6(""),
useShinyjs(),
extendShinyjs(text = jscode, functions = c("closeWindow")),
actionButton("close", "Close App",style="color: red; font-size: 14px; border-color: #2e6da4;width:320px"),
h6(""),h6(""),h6("")) 
  })
  ######################################################################
  ### Codes for CDEA RI
  ######################################################################
  output$limitUseRI <- renderText({""})
  
  
  valuesout0CRI<-reactiveValues(MolID0=NULL)
  
  
  
  #################### Read Accross
  valuesoutCRI<-reactiveValues(MolID=NULL,RI=NULL,libF=NULL)
  output$TheSMILES<- renderText(
    paste("Canonical SMILES: ",as.character(docanoSingleRI(as.character(input$SMILES))[[1]]))
  )
  
  observeEvent(input$butt1CRI,{
    
    compound=unlist(strsplit(input$caption,"\n"))
    print(compound)
    n<-length(compound)
    Output<-data.frame("ID"=1:n,"SMILES"=NA,"tR"=NA,"Exp RI"=NA,"Pred RI"=NA,"Probability Level 1"=NA,
                       "Probability Level 2"=NA,"pH chemical Space"=NA,"RI reliability"=NA,"AD model"=NA, check.names = FALSE)
    #################################################################
    #### Conditions 
    ################################################################# 
    if(n>limitUse1){
      output$limitUseRI <- renderText({paste0("Input is more than ", limitUse1, " enteries...Denied")})
    }else{
      output$limitUseRI <- renderText({paste0("Input is less than ", limitUse1, " enteries...Acceptable")})
      
      #################################################################
      #### Calculate experimental RI values 
      #################################################################  
      # if(input$method_S=="ExpRI"){
      inFile <- input$file1RI   
      if (is.null(inFile))
        return(NULL)
      calibrant<- read.csv(inFile$datapath, header = TRUE)
      RI_exp=NULL
      for (pp in 1:length(compound)) {
        tRT<-as.numeric(strsplit(compound[pp]," ")[[1]][2])
        Output$SMILES[pp]<-as.character(strsplit(compound[pp]," ")[[1]][1])
        Output$tR[pp]<-tRT
        RI_exp[pp]<-RTINDEA_cal(tRT,calibrant)
        Output$`Exp RI`[pp]<-RI_exp[pp]
      }
      valuesoutCRI$RI<-RI_exp
      # }
      
      
      #################################################################
      #### Calculate reliability values 
      ################################################################# 

      
      do.call("file.remove", list(list.files(paste("www","HashFP",sep = "/"), full.names = TRUE)))
      
      simo<-NULL
      MW<-NULL
      progress <- shiny::Progress$new()
      progress$set(message = "Running job...", value = 0)
      nameID=NULL
      comSus=NULL
      for (pp in 1:length(compound)) {
        out1<-FPQdes(as.character(compound[pp]))
        out<-FPFundes(as.character(compound[pp]))
        out<-cbind(out1,out)
        fingerprintsMP<-out
        comSus[[pp]]<-which(fingerprintsMP==1) 
        nameID[pp]<-paste("M",pp,sep = "")
        
        XG2<-lapply(split(as.list(fingerprintsMP), cut(1:ncol(fingerprintsMP), 101, labels = FALSE)), as.data.frame)
        nameExport<-paste(getwd(),"/www/","HashFP","/S",pp,".png",sep="")  
        png(nameExport,
            width     = 5,
            height    = 5,
            units     = "in",
            res       = 100,
            pointsize = 2)
        plot(c(0, 102), c(0, 101), type= "n", xlab = "", ylab = "",axes = TRUE,yaxt="n",xaxt="n")
        for (i in 1:101) {
          data<-XG2[[i]]
          for(j in 1:length(data)) {
            rect(j-0.5, i-1, j, i,col = if(data[1,j]==1){"black"}else{"white"},border = NA) # transparent
          }
        }
        dev.off()
        graphics.off()
        #}
        unlink(paste("www","HashFP","temp.mol",sep = "/"), recursive=TRUE)
        
      }
      width <- 40
      height <- 40
      
      progress$inc(0.30,detail = "Chemical Fingerprints are calculated ...")
      images_names <- list.files(paste("www/","HashFP/",sep=""),pattern = "png")
      
      suspec <- extract_feature(dir_path =paste("www/","HashFP",sep=""), width = 40, height = 40, 
                                AD0 = TRUE,AD1 = FALSE,add_label = TRUE, cols="grey")
      
      suspect<- suspec$aa
      suspect_name<- suspec$bb
      suspect_data <- data.matrix(suspect)
      suspect_x <- t(suspect[,-1])
      suspect_y <- suspect[,1]
      suspect_array<- suspect_x
      dim(suspect_array) <- c(40, 40, 1, ncol(suspect_x))
      
      suspect_probs <- predict(modelpH, suspect_array)
      predicted_labels_suspect <- max.col(t(suspect_probs)) - 1
      table(suspect_data[, 1], predicted_labels_suspect)
      sum(diag(table(suspect_data[, 1], predicted_labels_suspect)))/length(suspect_data[, 1])
      
      io<- 1
      labso<- as.array(0)
      labso[1:length(predicted_labels_suspect)]<- as.array(0)
      for (io in 1:length(predicted_labels_suspect)){
        if (predicted_labels_suspect[io]==0){label <- labso[io]<- 'Level 1 (RI may not be affected by pH)'}
        if (predicted_labels_suspect[io]==1){label <- labso[io]<- 'Level 2 (RI can be affected by pH)'}
        if (predicted_labels_suspect[io]==2){label <- labso[io]<- 'Outside of the pH chemical space domain'}
      }
      
      
      suspect_name2<-gsub("S","\\1",suspect_name)
      suspect_name2<-as.numeric(gsub(".png","\\1",suspect_name2))
      suspect_propability<- data.frame("index"=suspect_name2,
                                       "Name"=suspect_name,"Probability_Level1"=suspect_probs[1,],
                                       "Probability_Level2"=suspect_probs[2,],"Probability_OutAD"=suspect_probs[3,],
                                       "Predicted_label_cat"=labso,check.names = FALSE)
      newdata <- suspect_propability[order(suspect_propability$index),] 
      newdata$Probability_Level1=round(newdata$Probability_Level1,digit=3)
      newdata$Probability_Level2=round(newdata$Probability_Level2,digit=3)
      newdata$Probability_OutAD=round(newdata$Probability_OutAD,digit=3)
      
      Output$`Probability Level 1`<-newdata$Probability_Level1
      Output$`Probability Level 2`<-newdata$Probability_Level2
      Output$`pH chemical Space`<-newdata$Probability_OutAD
      Output$`RI reliability`<-newdata$Predicted_label_cat
      
      #################################################################
      #### Calculate Chemical Space
      #################################################################
      FPmatSpect2<-NULL
      ADSim<-NULL
      ADSim2<-NULL
      jj=1
      if(input$method_pRI=="pos"){
        for (jj in 1:length(compound)) {
          FPmatSpect<-NULL
          for (pp in 1:2290) {
            FPmatSpect[pp]<-CompareFP(FPmat[[pp]],comSus[[jj]])
          }
          FPmatSpect2[[jj]]<-FPmatSpect
          
          ADSim2[jj]<- MedTP5(FPmatSpect2[[jj]])[1]
          if(MedTP5(FPmatSpect2[[jj]])[1]>=0.3){
            ADSim[jj]<-"Within the Chemical Space Domain"
          }else{
            ADSim[jj]<-"Outside the Chemical Space Domain"
          }
        }
      }else if(input$method_pRI=="neg"){
        
        for (jj in 1:length(compound)) {
          FPmatSpect<-NULL
          for (pp in 1:728) {
            FPmatSpect[pp]<-CompareFP(FPmatneg[[pp]],comSus[[jj]])
          }
          FPmatSpect2[[jj]]<-FPmatSpect
          
          ADSim2[jj]<- MedTP5(FPmatSpect2[[jj]])[1]
          if(MedTP5(FPmatSpect2[[jj]])[1]>=0.3){
            ADSim[jj]<-"Within the Chemical Space Domain"
          }else{
            ADSim[jj]<-"Outside the Chemical Space Domain"
          }
        }
      }
      Output$`AD model`<-ADSim
      #################################################################
      #### Check if LC amendable
      #################################################################  
      suspect_probsChr <- predict(modelChrom, suspect_array)
      suspect_name2<-gsub("S","\\1",suspect_name)
      suspect_name2<-as.numeric(gsub(".png","\\1",suspect_name2))
      suspect_prred<- data.frame("index"=suspect_name2,
                                 "Name"=suspect_name,
                                 "ProbGC"=round(suspect_probsChr[1,],3),
                                 check.names = FALSE)
      newdata0 <- suspect_prred[order(suspect_prred$index),] 
      rownames(newdata0)<-1:nrow(suspect_prred)
      IOP<-which(newdata0$ProbGC>=0.7)
      if (!(rapportools::is.empty(IOP))){
        Output$`AD model`[IOP]<-"Analyte may not be appropriate for LC"
      }
      
      #################################################################
      #### Calculate predicted RI values
      #################################################################
      if(input$method_pRI=="pos"){
        MXX=2291.176 # maximum RI observed C(n=23)DEA
        suspect_predRIs <- as.numeric (predict(modelpredRI, suspect_array))*MXX
        suspect_name2<-gsub("S","\\1",suspect_name)
        suspect_name2<-as.numeric(gsub(".png","\\1",suspect_name2))
        suspect_prred<- data.frame("index"=suspect_name2,
                                   "Name"=suspect_name,
                                   "Pred RI"=round(suspect_predRIs,digit=3),
                                   check.names = FALSE)
        newdata0 <- suspect_prred[order(suspect_prred$index),] 
        Output$`Pred RI`<-newdata0$`Pred RI`
        

      }else if(input$method_pRI=="neg"){
        MXX=1774.07 # maximum RI observed C(n=23)DEA
        suspect_predRIs <- as.numeric (predict(modelpredRIhpH, suspect_array))*MXX
        suspect_name2<-gsub("S","\\1",suspect_name)
        suspect_name2<-as.numeric(gsub(".png","\\1",suspect_name2))
        suspect_prred<- data.frame("index"=suspect_name2,
                                   "Name"=suspect_name,
                                   "Pred RI"=round(suspect_predRIs,digit=3),
                                   check.names = FALSE)
        newdata0 <- suspect_prred[order(suspect_prred$index),] 
        Output$`Pred RI`<-newdata0$`Pred RI`
      }
      
      ########### deep learner output
      
      ################################  
      progress$close()
      
      
    }
    #################################################################
    #### prepare output
    #################################################################
    valuesoutCRI$MolID<-Output
    
  })
  
  
  #### show the read across table
  output$contents1CRI <- DT::renderDataTable({
    DT::datatable(valuesoutCRI$MolID,
                  filter = list(position = 'top', clear = FALSE),
                  options = list(lengthMenu = c(10, 30, 50), pageLength = 10, searchCols = list(NULL))
                  , rownames = FALSE,
                  escape = FALSE, selection = 'none') })
  #### donwload the results of read across   
  output$downsusCRI <- downloadHandler(
    filename = function() {
      paste('Candidates', '.csv', sep='')
    },
    content = function(file) {
      ## write the info here
      write.csv(valuesoutCRI$MolID, file, sep = ",",row.names = FALSE,na = "")
    }      
  )  
  
  
  #################################################################
  #### Search Database
  #################################################################
  observeEvent(input$butt_libCRI,{
    
    code<-c("[M]+"="Mp","[M+H]+"="MpH","[M+NH4]+"="MpNH",
            "[M+Na]+"="MpNa","[M+K]+"="MpK","[M+CH3OH+H]+"="MpAct","[M]2+"="Mp2",
            "[M+H]2+"="MpH2p","[M]-"="Mm","[M-H]-"="MmH","[M-H2O-H]-"="MmH2Omh","[M+Cl]-"="MpClm",
            "[M+HCOOH-H]-"="MpFAmH","[M+CH3COOH-H]-"="MpActmH","[2M-H]-"="2MmH",
            "[M]2-"="M2m","[M-2H]-"="Mm2H")
    MZSusV<-as.numeric(input$MZSuspect)
    code1<-input$IonTypeCRI
    MZacc33<-as.numeric(input$MZacc)
    
    databaseRI3<-databaseRI[,which(colnames(databaseRI)==(names(code)[which(code==code1)]))]
    # print(databaseRI3)
    if(MZacc33 > 1 || MZacc33 < -1){
      output$warn1<-renderText({"Broad mass accuracy value. Try to narrow down the mass accuracy value."})
      valuesoutCRI$libF<-NULL
    }else {
    valuesoutCRI$libF<-databaseRI[which(abs(databaseRI3-MZSusV)<=as.numeric(input$MZacc)),c(1,2,3,5,9,7,15,16,17,18,19)]
    }
    
    
  })
  
  #################################################################
  #### Search Database
  #################################################################
  output$contents_libCRI <- DT::renderDataTable({
    DT::datatable(valuesoutCRI$libF,
                  filter = list(position = 'top', clear = FALSE),
                  options = list(lengthMenu = c(10, 30, 50), pageLength = 10, searchCols = list(NULL))
                  , rownames = FALSE,
                  escape = FALSE, selection = 'none') })
  #### donwload the results of read across   
  output$downsus_libCRI <- downloadHandler(
    filename = function() {
      paste('LibFound', '.csv', sep='')
    },
    content = function(file) {
      ## write the info here
      write.csv(valuesoutCRI$libF, file, sep = ",",row.names = FALSE,na = "")
    }      
  )  
  
  #################################################################
  #### few examples/templates
  #################################################################
  
  output$CDF1<- renderText({
    "Example: SMILES tR"
  })
  output$CDF2<- renderText({
    "c1ccc(c(c1)CC(=O)OCC(=O)O)Nc2c(cccc2Cl)Cl 9.2"
  })
  output$CDF3<- renderText({
    "COC1=C(OC)C=C(C(C2=C(N)N=C(N)N=C2)=O)C=C1OC 5.49"
  })
  output$CDF4<- renderText({
    "CCCCCCCCn1c(=O)c(c(s1)Cl)Cl 12.49"
  })
  
  
  SuslistCRI <- reactive({ 
    # Downloadable csv of +ESI ----
    data.frame("N"=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23),
               "Name"=c("C0DEA","C1DEA","C2DEA","C3DEA","C4DEA","C5DEA","C6DEA","C7DEA",
                        "C8DEA","C9DEA","C10DEA","C11DEA","C12DEA","C13DEA","C14DEA","C15DEA",
                        "C16DEA","C17DEA","C18DEA","C19DEA","C20DEA","C21DEA","C22DEA","C23DEA"),
               "SMILES"=c("C(=O)N(CCO)CCO","CC(=O)N(CCO)CCO", "CCC(=O)N(CCO)CCO", "CCCC(=O)N(CCO)CCO", "CCCCC(=O)N(CCO)CCO", "CCCCCC(=O)N(CCO)CCO",  
                          "CCCCCCC(=O)N(CCO)CCO","CCCCCCCC(=O)N(CCO)CCO", "CCCCCCCCC(=O)N(CCO)CCO", "CCCCCCCCCC(=O)N(CCO)CCO", "CCCCCCCCCCC(=O)N(CCO)CCO",  
                          "CCCCCCCCCCCC(=O)N(CCO)CCO","CCCCCCCCCCCCC(=O)N(CCO)CCO", "CCCCCCCCCCCCCC(=O)N(CCO)CCO", "CCCCCCCCCCCCCCC(=O)N(CCO)CCO",  
                          "CCCCCCCCCCCCCCCC(=O)N(CCO)CCO","CCCCCCCCCCCCCCCCC(=O)N(CCO)CCO", "CCCCCCCCCCCCCCCCCC(=O)N(CCO)CCO",  
                          "CCCCCCCCCCCCCCCCCCC(=O)N(CCO)CCO","CCCCCCCCCCCCCCCCCCCC(=O)N(CCO)CCO", "CCCCCCCCCCCCCCCCCCCCC(=O)N(CCO)CCO", 
                          "CCCCCCCCCCCCCCCCCCCCCC(=O)N(CCO)CCO", "CCCCCCCCCCCCCCCCCCCCCCC(=O)N(CCO)CCO", "CCCCCCCCCCCCCCCCCCCCCCCC(=O)N(CCO)CCO)"),
               "tR"=c(1.28, 1.35, 2.51, 3.39, 4.38, 5.46, 6.68, 8.01, 9.08, 10.12, 10.94, 11.71, 12.36, 12.95, 13.39, 13.86, 14.16, 14.51, 14.74, 14.99, 15.34, 15.51, 15.79, 16.13),
               check.names = FALSE,stringsAsFactors =TRUE) 
    
  })
  output$downloadRIfile <- downloadHandler(
    filename = function() {
      paste('Template_calibrants', '.csv', sep='')
    },
    content = function(file) {
      ## write the info here
      write.csv(SuslistCRI(), file, sep = ",",row.names = FALSE,na = "")
    }      
  ) 
  #################################################################
  #### About page
  #################################################################
  output$ExportCit <- downloadHandler(
    filename = "acs.analchem.2c02893.ris",
    content = function(file) {
      file.copy("www/acs.analchem.2c02893.ris", file)
      
    }    
  )  
  output$manRI <- downloadHandler(
    filename = "Manual_RI_CDEA.pdf",
    content = function(file) {
      file.copy("www/manual.pdf", file)
      
    }    
  )  
  
  
  
  observeEvent(input$close, {
    js$closeWindow()
    stopApp()
  })
  
  
  
  })
################################################################
#### This app is written by Dr. Reza Aalizadeh (C) 2022 ########
#### For any questions, drop email in raalizadeh@chem.uoa.gr ###
################################################################

