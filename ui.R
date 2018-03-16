library(shinydashboard)
#library(shinyIncubator)
library(shiny)
library(shinyBS)
library(plotly)
library(d3heatmap)
library(shinyjs)
library(rglwidget)

ui <- dashboardPage(
  dashboardHeader(title = "Ligand-Receptor Pairs",titleWidth = 350),
  dashboardSidebar(width = 350,
                   div(style="overflow-y: scroll"),
                   tags$head(tags$style(HTML(".sidebar { height: 170vh; overflow-y: auto; }" ))),
                   #uiOutput("projects"),
                   sidebarMenu(
                     menuItem('Compare datasets', tabName = 'compare', icon = icon('hand-o-right')),
                     menuItem("View TSNE", tabName = "tsne", icon = icon("hand-o-right"))
                     #menuItem("scRNA data", tabName = "dashboard", icon = icon("hand-o-right"))
                     
                     
                     #          menuSubItem("PCA Plot", tabName = "dashboard"),
                     #          menuSubItem('Display Variances', tabName = 'var'),
                     #          menuSubItem('Show 3D plot', tabName = '3dplot'))
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
                  )
                  ),

              box(
                width = 10, status = "primary",solidHeader = TRUE,
                title = "Ligand Receptor pairs",
                DT::dataTableOutput('pairs_res')
                  )#end of box
              ),#end of tabitem
      ######################################################################
      tabItem(tabName = "compare",
              box(width = 6, status = "primary",solidHeader = TRUE,title = "Ligand Selection Panel",
                  selectInput("ligand", "Select Experiment type",c('RNA-Seq' = "rna",'Single Cell' = "scrna", 'Microarray' = "microarray")),
                  uiOutput("ligprj"),
                  uiOutput("ligtype"),
                  conditionalPanel(
                    condition = "input.ligand == 'rna' | input.ligand == 'microarray'" ,
                    sliderInput("explig", label = "Set Expression threshold", min =6,max = 12, value = 6)),
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
                  uiOutput("imp_pdf")
              )
      )#end of tabitem
      ######################################################################
    )#end of tabitems
  )#end of dashboardbosy
)#end of dashboard page


