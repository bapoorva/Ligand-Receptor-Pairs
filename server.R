library(shiny)
library(shinyBS)
library(RColorBrewer)
library(ggplot2)
library(png)
library(d3heatmap)
library(dplyr)
library(tidyr)
library(plotly)
library(shinyjs)
library(htmlwidgets)
library(DT)
library(shinyRGL)
library(rgl)
library(rglwidget)
library(ggrepel)
library(readxl)
library(Biobase)
library(pathview)
library(gage)
library(gageData)
library(org.Mm.eg.db)

#load entrez id's for kegg pathway
data(kegg.sets.mm)
#load indexes for signaling and metabolic pathways
data(sigmet.idx.mm)
#get entrez id's for signaling and metabolic pathwaysmera
kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]

data(go.sets.mm)
data(go.subs.mm)

server <- function(input, output) {
  
  ###################################################
  ###################################################
  ####### LOAD EXCEL AND POPULATE DROP DOWN #########
  ###################################################
  ###################################################
  
  #Read the parameter file
  readexcel = reactive({
    file = read.csv("data/param.csv")
  })
  
  #Get Project list and populate drop-down
  output$projects = renderUI({
    excel=readexcel()
    excel=excel[excel$type=="scrna",]
    prj=excel$projects
    selectInput("projects","Select a project",as.list(sort(as.character(prj))))
  })
  
  output$recprj = renderUI({
    excel=readexcel()
    excel=excel[excel$type==input$receptor,]
    prj=excel$projects
    selectInput("recprj","Select a project",as.list(sort(as.character(prj))))
  })
  
  output$ligprj = renderUI({
    excel=readexcel()
    excel=excel[excel$type==input$ligand,]
    prj=excel$projects
    selectInput("ligprj","Select a project",as.list(sort(as.character(prj))))
  })
  ###################################################
  ###################################################
  ############# READ OR LOAD DATA FILES #############
  ###################################################
  ###################################################
  recfileload <- reactive({
    inFile = paste('data/',as.character(input$recprj),'.RData',sep = '')
    load(inFile)
    loaddata=results
    return(loaddata)
  })
  ligfileload <- reactive({
    inFile = paste('data/',as.character(input$ligprj),'.RData',sep = '')
    load(inFile)
    loaddata=results
    return(loaddata)
  })
  
  recfileread <- reactive({
    inFile = paste('data/',as.character(input$recprj),'.csv',sep = '')
    loaddata=read.csv(inFile)
    return(loaddata)
  })
  
  ligfileread <- reactive({
    inFile = paste('data/',as.character(input$ligprj),'.csv',sep = '')
    loaddata=read.csv(inFile)
    return(loaddata)
  })
  
  ###################################################
  ###################################################
  ### POPULATE DROPDOWN FOR MAINEFFECT AND CLUSTER ##
  ###################################################
  ###################################################
  output$rectype = renderUI({
    if(input$receptor=="scrna"){
      loaddata=recfileread()
      maineffect=unique(loaddata$ident)
    }else if(input$receptor=="rna" | input$receptor=="microarray"){
      results=recfileload()
      pd=pData(results$eset)
      maineffect=pd$maineffect
    }
    selectInput("rectype","Select subtype (maineffect for RNAseq and cluster for scRNA)",as.list(sort(as.character(maineffect))),selected = 1)
  })
  
  output$ligtype = renderUI({
    if(input$ligand=="scrna"){
      loaddata=ligfileread()
      maineffect=unique(loaddata$ident)
    }else if(input$ligand=="rna" | input$ligand=="microarray"){
      results=ligfileload()
      pd=pData(results$eset)
      maineffect=pd$maineffect
    }
    selectInput("ligtype","Select subtype (maineffect for RNAseq and cluster for scRNA)",as.list(sort(as.character(maineffect))),selected = 1)
  })
  
  ###################################################
  ###################################################
  ####### COMPARE EFFECTS AND DISPLAY RESULTS #######
  ###################################################
  ###################################################
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  ligrecpairs = reactive({
    if(input$receptor=="scrna"){
      scrna=recfileread()
      scrna=scrna[scrna$ident==input$rectype,]
      rownames(scrna)=scrna$X
      scrna=scrna %>% dplyr::select(-X:-nGene)
      recumi=input$recumi
      recsamp=(input$recsamp)/100
      keep= colSums(scrna>recumi)>=recsamp*dim(scrna)[1]
      scrna2=scrna[,keep]
      rec_avg=NULL
      rec_genes=colnames(scrna2)
    }else if(input$receptor=="rna" | input$receptor=="microarray"){
      results=recfileload()
      loaddata=as.data.frame(exprs(results$eset))
      pd=pData(results$eset)
      #minexpr=unique(pd$minexpr)
      minexpr=input$exprec
      sel=input$rectype
      samp=as.character(pd$sample_name[pd$maineffect==sel])
      loaddata2=loaddata %>% dplyr::select(samp)
      loaddata2$avg=rowMeans(loaddata2)
      fd=fData(results$eset)
      fd$id=rownames(fd)
      loaddata2$id=rownames(loaddata2)
      loaddata2=inner_join(loaddata2,fd,by="id")
      loaddata2=loaddata2[loaddata2$avg > minexpr,]
      rec_avg=loaddata2 %>% dplyr::select(avg,SYMBOL)
      rec_genes=loaddata2$SYMBOL
    }

    if(input$ligand=="scrna"){
      scrna=ligfileread()
      scrna=scrna[scrna$ident==input$ligtype,]
      rownames(scrna)=scrna$X
      scrna=scrna %>% dplyr::select(-X:-nGene)
      ligumi=input$ligumi
      ligsamp=(input$ligsamp)/100
      keep= colSums(scrna>ligumi)>=ligsamp*dim(scrna)[1]
      scrna2=scrna[,keep]
      lig_avg=NULL
      lig_genes=colnames(scrna2)
    }else if(input$ligand=="rna" | input$ligand=="microarray"){
      results=ligfileload()
      loaddata=as.data.frame(exprs(results$eset))
      pd=pData(results$eset)
      #minexpr=unique(pd$minexpr)
      minexpr=input$explig
      sel=input$ligtype
      samp=as.character(pd$sample_name[pd$maineffect==sel])
      loaddata2=loaddata %>% dplyr::select(samp)
      loaddata2$avg=rowMeans(loaddata2)
      fd=fData(results$eset)
      fd$id=rownames(fd)
      loaddata2$id=rownames(loaddata2)
      loaddata2=inner_join(loaddata2,fd,by="id")
      loaddata2=loaddata2[loaddata2$avg > minexpr,]
      lig_avg=loaddata2 %>% dplyr::select(avg,SYMBOL)
      lig_genes=loaddata2$SYMBOL
    }
    rl=read.csv("data/lig-rec.csv")
    if(is.null(lig_avg)==F){rl=left_join(rl,lig_avg,by=c("Ligand.ApprovedSymbol"="SYMBOL")) %>% dplyr::rename(Ligand_AvgExpr=avg)}
    if(is.null(rec_avg)==F){rl=left_join(rl,rec_avg,by=c("Receptor.ApprovedSymbol"="SYMBOL")) %>% dplyr::rename(Receptor_AvgExpr=avg)}
    rl=rl[rl$Ligand.ApprovedSymbol %in% lig_genes & rl$Receptor.ApprovedSymbol %in% rec_genes,]
    
    
    if(input$liggene ==T){
      validate(
              need(input$liggeneli, "Please Upload genelist")
            )
      lgene=input$liggeneli
      lgenes=read.table(lgene$datapath,stringsAsFactors = F)#get complete gene list as string
      lgenes=as.vector(lgenes$V1)
      lgenes=tolower(lgenes)
      lgenes=firstup(lgenes)
    }else{lgenes=NULL}
    
    if(input$recgene ==T){
      validate(
        need(input$recgeneli, "Please Upload genelist")
      )
      rgene=input$recgeneli
      rgenes=read.table(rgene$datapath,stringsAsFactors = F)#get complete gene list as string
      rgenes=as.vector(rgenes$V1)
      rgenes=tolower(rgenes)
      rgenes=firstup(rgenes)
    }else{rgenes=NULL}
    
    if(is.null(lgenes)==F & is.null(rgenes)==T ){
      rl=rl[rl$Ligand.ApprovedSymbol %in% lgenes,]
    }else if(is.null(lgenes)==T & is.null(rgenes)==F){
      rl=rl[rl$Receptor.ApprovedSymbol %in% rgenes,]
    }else if(is.null(lgenes)==F & is.null(rgenes)==F ){
      rl=rl[(rl$Receptor.ApprovedSymbol %in% rgenes) & (rl$Ligand.ApprovedSymbol %in% lgenes),]
    }else{rl=rl}
    
    validate(
      need(nrow(rl)>=1, "No Ligand-Receptor Pairs for the combination chosen")
    )
    recid=toupper(rl$Receptor.ApprovedSymbol)
    ligid=toupper(rl$Ligand.ApprovedSymbol)
    urlr= paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",recid,sep = "")
    urll= paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",ligid,sep = "")
    rl$Receptor.ApprovedSymbol=paste0("<a href='",urlr,"'target='_blank'>",rl$Receptor.ApprovedSymbol,"</a>")
    rl$Ligand.ApprovedSymbol=paste0("<a href='",urll,"'target='_blank'>",rl$Ligand.ApprovedSymbol,"</a>")
    return(rl)
  })



  #print data TABLE
  output$ligrecpairs = DT::renderDataTable({
    input$receptor
    input$ligand
    input$recprj
    input$ligprj
    input$rectype
    input$ligtype
    input$liggene
    input$recgene
    input$liggeneli
    input$recgeneli
    DT::datatable(ligrecpairs(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=FALSE,escape = F,selection = list(mode = 'single', selected =1))
  })

  output$dwldtab = renderUI({
    downloadButton('dwldtable','Download')
  }) 
  
  output$dwldtable <- downloadHandler(
    filename = function() { paste(input$ligprj,"_",input$ligtype,"_",input$recprj,"_",input$rectype,'.csv', sep='') },
    content = function(file) {
      write.csv(ligrecpairs(), file,row.names=FALSE)
    })
  
  ###################################################
  ###################################################
  ####### UPLOAD GENELISTS AND SELECT CLUSTER #######
  ###################################################
  ###################################################
  output$list1 <- renderUI({
    fileInput('genelist1', 'Upload Receptor Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
  })
  
  output$list2 <- renderUI({
    fileInput('genelist2', 'Upload Ligand Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
  })
  
  output$clust1 <- renderUI({
    selectInput("clust1","Pick cluster1",c(0:11),selected=1)
  })
  
  output$clust2 <- renderUI({
    selectInput("clust2","Pick cluster2",c(0:11),selected=2)
  })
  
  output$list1.1 <- renderUI({
    fileInput('genelist1.1', 'Upload Gene List1',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
  })
  
  output$list2.1 <- renderUI({
    fileInput('genelist2.1', 'Upload Gene List2',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
  })
  
  output$clust1.1 <- renderUI({
    selectInput("clust1.1","Pick cluster1",c(0:11),selected=1)
  })
  
  output$clust2.1 <- renderUI({
    selectInput("clust2.1","Pick cluster2",c(0:11),selected=2)
  })
  ###################################################
  ###################################################
  ###### DISPLAY TSNE PLOT FOR SELECTED PROJECT######
  ###################################################
  ###################################################
  output$tsneprj = renderUI({
    excel=readexcel()
    excel=excel[excel$type=="scrna",]
    prj=excel$projects
    selectInput("tsneprj","Select a project",as.list(sort(as.character(prj))))
  })
  
  output$imp_pdf <- renderUI({
    if(input$seuratclus==F){
    PDFfile=paste(input$tsneprj,"_TSNE.pdf",sep="")
    }else{
      PDFfile=paste(input$tsneprj,"_Clusters.pdf",sep="")
    }
    #PDFfile="test.pdf"
    tags$iframe(src=PDFfile,width="100%",height="1000px")
  })
  ###################################################
  ###################################################
  #### Load lig-receptor list and display results ###
  ###################################################
  ###################################################
  
  datasetInput = reactive({
    prj=paste("data/",input$projects,".csv",sep="")
    rl=read.csv("data/lig-rec.csv")
    my.data=read.csv(prj)
    # clust1.1=input$clust1.1
    # clust2.1=input$clust2.1
    # 
    result=data.frame()
    res=data.frame()
    for(i in 0:(length(unique(my.data$clust))-1)){
      for(j in 0:(length(unique(my.data$clust))-1)){
        if(i!=j){
          test=my.data[my.data$clust==i | my.data$clust==j,]
          R_c1=test[test$clust==i ,(colnames(test) %in% rl$Receptor.ApprovedSymbol)]
          L_c2=test[test$clust==j , (colnames(test) %in% rl$Ligand.ApprovedSymbol)]
          keep1 = colSums(R_c1>1)>=.5*dim(R_c1)[1]
          keep2 = colSums(L_c2>1)>=.5*dim(L_c2)[1]

          R_c1=R_c1[,keep1]
          L_c2=L_c2[,keep2]
          res=rl[(rl$Ligand.ApprovedSymbol %in% colnames(L_c2)) & (rl$Receptor.ApprovedSymbol %in% colnames(R_c1)),]

        }
        else{}
        if(nrow(res)!=0){
          res$Receptor_cluster=i
          res$Lig_cluster=j
          result=rbind(result,res)
        }else{result=result}
      }
    }
    result=result[result$Receptor_cluster!=result$Lig_cluster,]
    
    if(input$clust=="clust" & input$gene=="allgene"){
      clusters=c(input$clust1,input$clust2)
      result=result[(result$Receptor_cluster %in% clusters) & (result$Lig_cluster%in% clusters),]
    }else if(input$clust=="clust" & input$gene=="genelist"){
      clusters.1=c(input$clust1.1,input$clust2.1)
      result=result[(result$Receptor_cluster %in% clusters.1) & (result$Lig_cluster%in% clusters.1),]
      }else{result=result
      }
    
    if(input$gene=="genelist"){
      if(input$clust=="all"){
      g1=input$genelist1
      g2=input$genelist2
    }else if(input$clust=="clust"){
      g1=input$genelist1.1
      g2=input$genelist2.1
    }
      genes=read.table(g1$datapath,stringsAsFactors = F)#get complete gene list as string
      g1=as.vector(genes$V1)
      g1=tolower(g1)
      firstup <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        x
      }
      g1=firstup(g1)
      genes2=read.table(g2$datapath,stringsAsFactors = F)
      g2=as.vector(genes2$V1)
      g2=tolower(g2)
      g2=firstup(g2)
      result=result[(result$Receptor.ApprovedSymbol %in% g1) & (result$Ligand.ApprovedSymbol %in% g2),]
    }else{
      result=result
    }
   return(result)
  })

  #print data TABLE
  output$pairs_res = DT::renderDataTable({
    input$clust1
    input$clust2
    input$genelist1
    input$genelist2
    DT::datatable(datasetInput(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=FALSE,caption= "Result")
  })
  
  ###################################################
  ###################################################
  ############ run gage for kegg pathways ###########
  ###################################################
  ###################################################
  gageres <- reactive({
    #get receptor list and annotate it to get entrez ids
    tab=ligrecpairs()
    tab$pair=tab$Pair.Name
    tab= tab %>% separate(pair,c("lig","rec"),sep="_")
    res <- AnnotationDbi::select(org.Mm.eg.db, keys=as.character(tab$rec), columns=c("ENTREZID","SYMBOL"), keytype="SYMBOL")
    res <- subset(res,!duplicated(res$ENTREZID))
    final_tab<- left_join(tab,res, by=c("rec"="SYMBOL"))
    final_tab=unique(final_tab)
    final_tab$ENTREZID=paste0("ncbi-geneid:",final_tab$ENTREZID,sep="")
    
    #Convert entrez ids into kegg ids
    keggids=(keggConv("mmu", final_tab$ENTREZID))
    # keggids2=as.data.frame(keggids)
    # keggids2$name=names(keggids)
    # final_tab=left_join(final_tab,keggids2,by=c("ENTREZID"="name"))
    # final_tab=unique(final_tab)
    
    #Find pathways
    res=keggLink("pathway",keggids)
    res2=as.data.frame(res)
    res2$id=names(res)
    res3=res2[res2$res %in% unique(res2$res[duplicated(res2$res)]),]
    tab=as.data.frame(table(res3$res))
    
    
  })
  ###################################################
  ###################################################
  ################ PATHWAY ANALYSIS #################
  ###################################################
  ###################################################
  
}#end of server
