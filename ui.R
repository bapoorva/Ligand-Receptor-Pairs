library(shinydashboard)
#library(shinyIncubator)
library(shiny)
library(shinyBS)
library(plotly)
library(d3heatmap)
library(shinyjs)
library(rglwidget)

ui <- dashboardPage(
  dashboardHeader(title = "Ligand-Receptor Pairs",titleWidth = 300),
  dashboardSidebar(width = 300,
                   div(style="overflow-y: scroll"),
                   tags$head(tags$style(HTML(".sidebar { height: 250vh; overflow-y: auto; }" ))),
                   #uiOutput("projects"),
                   sidebarMenu(
                     menuItem('Compare datasets', tabName = 'compare', icon = icon('hand-o-right')),
                     menuItem("View TSNE", tabName = "tsne", icon = icon("hand-o-right")),
                     menuItem("scRNA data", tabName = "dashboard", icon = icon("hand-o-right")),
                     menuItem("Pathway Analysis", tabName = "pathway", icon = icon("hand-o-right"))
         )#end of sidebar menu

  ),#end dashboardSidebar
  
  
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    useShinyjs(),
    tabItems(
      tabItem(tabName = "dashboard",
              box(width = 10, status = "primary",solidHeader = TRUE,title = "Controls",
                  uiOutput("projects"),
                  radioButtons("clust","Select one", c("All clusters"="all","Select Cluster"="clust"),selected = "clust"),
                  radioButtons("gene","Select one", c("All genes"="allgene","Enter Genelist"="genelist"),selected = "allgene"),

                  conditionalPanel(
                    condition = "input.clust == 'all' && input.gene == 'genelist'" ,
                    uiOutput("list1"),
                    uiOutput("list2")
                  ),
                  conditionalPanel(
                    condition = "input.clust == 'clust' && input.gene == 'allgene'" ,
                    uiOutput("clust1"),
                    uiOutput("clust2")
                  ),
                  conditionalPanel(
                    condition = "input.clust == 'clust' && input.gene == 'genelist'" ,
                    uiOutput("list1.1"),
                    uiOutput("list2.1"),
                    uiOutput("clust1.1"),
                    uiOutput("clust2.1")
                  ),
                  fluidRow(
                    column(6,checkboxInput("checksource", label = "Check to select by source", value = FALSE)),
                    column(6,checkboxInput("checkevi", label = "Check to select by evidence", value = FALSE)),
                    conditionalPanel(
                      condition = "input.checksource ==true",
                      column(6,uiOutput('source'))
                    ),
                    conditionalPanel(
                      condition = "input.checkevi ==true",
                      column(6,uiOutput('evidence'))
                    )
                    
                    
                  )
                  ),

              box(
                width = 10, status = "primary",solidHeader = TRUE,
                title = "Ligand Receptor pairs",
                DT::dataTableOutput('pairs_res')
                  )#end of box
              ),#end of tabitem
      ######################################################################
      ######################################################################
      tabItem(tabName = "compare",
              box(width = 6, status = "primary",solidHeader = TRUE,title = "Ligand Selection Panel",
                  selectInput("ligand", "Select Experiment type",c('RNA-Seq' = "rna",'Single Cell' = "scrna", 'Microarray' = "microarray")),
                  uiOutput("ligprj"),
                  uiOutput("ligtype"),
                  conditionalPanel(
                    condition = "input.ligand == 'rna' | input.ligand == 'microarray'" ,
                    sliderInput("explig", label = "Set Expression threshold", min =6,max = 12, value = 6)),
                  conditionalPanel(
                    condition = "input.ligand == 'scrna'" ,
                    sliderInput("ligumi", label = "Set UMI threshold", min =1,max = 25, value = 1),
                    sliderInput("ligsamp", label = "Set Percent Samples", min =0,max = 100, value = 50)),
                 checkboxInput("liggene", label = "Upload Gene List", value = FALSE),
    
      conditionalPanel(
        condition = "input.liggene ==true",
        fileInput('liggeneli', 'Upload Receptor Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
      )),
              box(width = 6, status = "primary",solidHeader = TRUE,title = "Receptor Selection Panel",
                  selectInput("receptor", "Select Experiment type",c('RNA-Seq' = "rna",'Single Cell' = "scrna", 'Microarray' = "microarray")),
                  uiOutput("recprj"),
                  uiOutput("rectype"),
                  conditionalPanel(
                    condition = "input.receptor == 'rna' | input.receptor == 'microarray'" ,
                    sliderInput("exprec", label = "Set Expression threshold", min =6,max = 12, value = 6)),
                  conditionalPanel(
                    condition = "input.receptor == 'scrna'" ,
                    sliderInput("recumi", label = "Set UMI threshold", min =1,max = 25, value = 1),
                    sliderInput("recsamp", label = "Set Percent Samples", min =0,max = 100, value = 50)),
                  checkboxInput("recgene", label = "Upload Gene List", value = FALSE),
                  conditionalPanel(
                    condition = "input.recgene ==true",
                    fileInput('recgeneli', 'Upload Receptor Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
                  )),
              box(width = 12, status = "primary",solidHeader = TRUE,title = "Ligand-Receptor pairs",
                  DT::dataTableOutput('ligrecpairs'),uiOutput("dwldtab"))
      ),#end of tabitem
      ######################################################################
      ######################################################################
      tabItem(tabName = "tsne",
              box(width = 12,height = 10, status = "primary",solidHeader = TRUE,title = "TSNE Plot",
                  uiOutput("tsneprj"),
                  checkboxInput("seuratclus", label = "Check to view Seurat Clusters", value = FALSE),
                  uiOutput("imp_pdf")
              )
      ),#end of tabitem
      ######################################################################
      ######################################################################
      tabItem(tabName = "pathway",
              box(width = 12, status = "primary",solidHeader = TRUE,title = "Ligand Receptor Pairs",
                  DT::dataTableOutput('rec')
              ),#end of box
              box(width = 12, status = "primary",solidHeader = TRUE,title = "KEGG Pathways",
                  DT::dataTableOutput('Keggpaths')
              ),
              box(width = 12,height = 12, status = "primary",solidHeader = TRUE,title = "Pathview",
                  plotOutput("plots")
              )
      )#end of tabitem
      ######################################################################
      ######################################################################
    )#end of tabitems
  )#end of dashboardbosy
)#end of dashboard page


