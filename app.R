library(shiny)
library(igraph)
library(visNetwork)

seedrows=function(seed_graph,seedvec){
  vertexlist=unlist(igraph::vertex_attr(seed_graph))
  rowindex=which(is.element(vertexlist,seedvec))
}

cand_net=function(g,cand,s2boutput,indexa,indexb,symbnodes){
  cindex=which(s2boutput$s2btable$symbol==cand)
  seeda=indexa[which(s2boutput$seedmat1[cindex,]>0)]
  seedb=indexb[which(s2boutput$seedmat2[cindex,]>0)]
  seeds=union(seeda,seedb)
  sp1=igraph::distances(g,v=cindex,to=igraph::V(g))
  sp2=igraph::distances(g,v=seeds,to=igraph::V(g))
  sp1[sp1==Inf]=igraph::vcount(g)
  sp2[sp2==Inf]=igraph::vcount(g)
  d2seeds=sp1[seeds]
  m=matrix(d2seeds,nrow=length(seeds),ncol=length(sp1),byrow=F)
  sp1=matrix(sp1,nrow=length(seeds),ncol=length(sp1),byrow=TRUE)
  sumsp=sp1+sp2-m
  midnodes=which(colSums(sumsp==0)>0)
  nodesindex=union(midnodes,seeds)
  nodesindex=union(nodesindex,cindex)
  g=set.vertex.attribute(g, "name",value=symbnodes$symbol)
  candg=induced_subgraph(g,nodesindex)
  nodes=as_data_frame(candg,what="vertices")
  names(nodes)="label"
  nodes$id=(1:length(nodesindex))
  edgemat=as_edgelist(candg,names=F)
  edges=data.frame(from=edgemat[,1],to=edgemat[,2])
  #nodes$color="grey"
  #nodes$color[is.element(nodes$label,symbnodes[seeda,1])]="red"
  #nodes$color[is.element(nodes$label,symbnodes[seedb,1])]="blue"
  nodes$group="node"
  nodes$group[is.element(nodes$label,symbnodes[seeda,1])]="A"
  nodes$group[is.element(nodes$label,symbnodes[seedb,1])]="B"
  nodes$group[nodes$label==cand]="S2B"
  for (i in 1:length(nodesindex)){
    nodes$value[i]=s2boutput$s2btable$bcount[s2boutput$s2btable$symbol==nodes$label[i]]  
  }
  
  output=list(nodes=nodes,edges=edges)
  
  
}

S2B=function(seed_graph,index1,index2,nrep,nrep2,meandist,symbnodes,genbet){
  
  bt=subS2B(seed_graph,index1,index2,meandist)
  pbt=rep(0,igraph::gorder(seed_graph))
  nscore=rep(0,igraph::gorder(seed_graph))
  deglist=igraph::degree(seed_graph)
  if (nrep2>0){
    rbt_matrix2=matrix(nrow=length(bt$allcount),ncol=nrep2)
    for (i in 1:nrep2){
      rindex1=sample(igraph::gorder(seed_graph),length(index1),replace=FALSE)
      rindex2=sample(igraph::gorder(seed_graph),length(index2),replace=FALSE)
      rbt=justS2B(seed_graph,rindex1,rindex2,meandist)
      nscore[rbt$allcount<bt$allcount]=nscore[rbt$allcount<bt$allcount]+1
      rbt_matrix2[,i]=rbt$allcount
    }
    nscore=nscore/nrep2
  } else {
    rbt_matrix2=matrix()
  }
  if (nrep>0){
    rbt_matrix=matrix(nrow=length(bt$allcount),ncol=nrep)
    for (i in 1:nrep){
      rg=igraph::sample_degseq(deglist,method="vl")
      rbt=justS2B(rg,index1,index2,meandist)
      pbt[rbt$allcount<bt$allcount]=pbt[rbt$allcount<bt$allcount]+1
      rbt_matrix[,i]=rbt$allcount
    }
    pbt=pbt/nrep
  } else {
    rbt_matrix=matrix()
  }
  bigvertexlist=igraph::vertex_attr(seed_graph)
  allstat=data.frame(protein=bigvertexlist[[1]],symbol=symbnodes,S2B=bt$allcount,bcount=bt$allcount*bt$maxS2B,Gbtw=genbet,S2Brank=rank(-bt$allcount),Gbtwrank=rank(-genbet),sscore=pbt, escore=nscore) ####
  #s2btable=makes2btable(allstat,seed_graph,index1,index2)
  list(s2btable=allstat,seedmat1=bt$smat1,seedmat2=bt$smat2,maxS2B=bt$maxS2B)
}

subS2B_old=function(seed_graph,index1,index2,meandist){
  betweencount=rep(0,igraph::gorder(seed_graph)) #lista del tama??o de seed_graph
  seedmat1=matrix(data=0,nrow=igraph::gorder(seed_graph),ncol=length(index1))
  seedmat2=matrix(data=0,nrow=igraph::gorder(seed_graph),ncol=length(index2))
  sp1=igraph::distances(seed_graph,v=index1,to=igraph::V(seed_graph))
  sp2=igraph::distances(seed_graph,v=index2,to=igraph::V(seed_graph))
  sp1[sp1==Inf]=igraph::vcount(seed_graph)
  sp2[sp2==Inf]=igraph::vcount(seed_graph)
  sp=sp1[,index2]
  maxbc=sum(sp>0 & sp<meandist)
  betweensub=union(index1,index2)
  for (i in 1:length(index1)){
    for (j in 1:length(index2)){
      m=sp1[i,index2[j]]
      if (m<meandist){
        sumsp=sp1[i,]+sp2[j,]
        nodelist=which(sumsp==m) # aquellos que cumplan la condicion, seran nodos presentes en un sh_path y seran a??adidos a nodelist
        betweencount[nodelist]=betweencount[nodelist]+1
        seedmat1[nodelist,i]=1
        seedmat2[nodelist,j]=1
      }else {
        nodelist=c(index1[i],index2[j])
        betweencount[nodelist]=betweencount[nodelist]+1
      }
    }
  }
  betweencount[index1]=betweencount[index1]-length(index2)
  betweencount[index2]=betweencount[index2]-length(index1)
  betweencount[intersect(index1,index2)]=betweencount[intersect(index1,index2)]+1
  betweencount=betweencount/maxbc
  list(allcount=betweencount,smat1=seedmat1,smat2=seedmat2,maxS2B=maxbc)
}

subS2B=function(seed_graph,index1,index2,meandist){
  betweencount=rep(0,igraph::gorder(seed_graph)) #lista del tama??o de seed_graph
  seedmat1=matrix(data=0,nrow=igraph::gorder(seed_graph),ncol=length(index1))
  seedmat2=matrix(data=0,nrow=igraph::gorder(seed_graph),ncol=length(index2))
  sp1=igraph::distances(seed_graph,v=index1,to=igraph::V(seed_graph))
  sp2=igraph::distances(seed_graph,v=index2,to=igraph::V(seed_graph))
  sp1[sp1==Inf]=igraph::vcount(seed_graph)
  sp2[sp2==Inf]=igraph::vcount(seed_graph)
  sp=sp1[,index2]
  maxbc=sum(sp>0 & sp<meandist)
  betweensub=union(index1,index2)
  for (i in 1:length(index1)){
      m=sp1[i,index2]
      betweencount[index1[i]]=betweencount[index1[i]]+sum(m>=meandist)
      betweencount[index2[m>=meandist]]=betweencount[index2[m>=meandist]]+1
      m[m>=meandist]=0

      mmat=matrix(m,nrow=length(index2),ncol=igraph::vcount(seed_graph),byrow=F)
      sp1mat=matrix(sp1[i,],nrow=length(index2),ncol=igraph::vcount(seed_graph),byrow=TRUE)
      sumsp=sp1mat+sp2-mmat
      #nodelist=which(sumsp==m) # aquellos que cumplan la condicion, seran nodos presentes en un sh_path y seran a??adidos a nodelist
      betweencount=betweencount+colSums(sumsp==0)
      seedmat1[(colSums(sumsp==0)>0),i]=1
      #ij=which(sumsp==0,arr.ind=T)
      seedmat2=seedmat2+t((sumsp==0))
      #}else {
      #  nodelist=c(index1[i],index2[j])
      #  betweencount[nodelist]=betweencount[nodelist]+1
      #}
  }
  seedmat2=1*(seedmat2>0)
  betweencount[index1]=betweencount[index1]-length(index2)
  betweencount[index2]=betweencount[index2]-length(index1)
  betweencount[intersect(index1,index2)]=betweencount[intersect(index1,index2)]+1
  betweencount=betweencount/maxbc
  list(allcount=betweencount,smat1=seedmat1,smat2=seedmat2,maxS2B=maxbc)
}

justS2B=function(seed_graph,index1,index2,meandist){
  betweencount=rep(0,igraph::gorder(seed_graph)) #lista del tama??o de seed_graph
  sp1=igraph::distances(seed_graph,v=index1,to=igraph::V(seed_graph))
  sp2=igraph::distances(seed_graph,v=index2,to=igraph::V(seed_graph))
  sp1[sp1==Inf]=igraph::vcount(seed_graph)
  sp2[sp2==Inf]=igraph::vcount(seed_graph)
  sp=sp1[,index2]
  maxbc=sum(sp>0 & sp<meandist)
  for (i in 1:length(index1)){
    m=sp1[i,index2]
    betweencount[index1[i]]=betweencount[index1[i]]+sum(m>=meandist)
    betweencount[index2[m>=meandist]]=betweencount[index2[m>=meandist]]+1
    m[m>=meandist]=0
    
    mmat=matrix(m,nrow=length(index2),ncol=igraph::vcount(seed_graph),byrow=F)
    sp1mat=matrix(sp1[i,],nrow=length(index2),ncol=igraph::vcount(seed_graph),byrow=TRUE)
    sumsp=sp1mat+sp2-mmat
    betweencount=betweencount+colSums(sumsp==0)
    
  }
  betweencount[index1]=betweencount[index1]-length(index2)
  betweencount[index2]=betweencount[index2]-length(index1)
  betweencount[intersect(index1,index2)]=betweencount[intersect(index1,index2)]+1
  betweencount=betweencount/maxbc
  list(allcount=betweencount,maxS2B=maxbc)
}

fastS2B=function(seed_graph,index1,index2,meandist){
  betweencount=rep(0,igraph::gorder(seed_graph)) #lista del tama??o de seed_graph
  sp1=igraph::distances(seed_graph,v=index1,to=igraph::V(seed_graph))
  sp2=igraph::distances(seed_graph,v=index2,to=igraph::V(seed_graph))
  sp1[sp1==Inf]=igraph::vcount(seed_graph)
  sp2[sp2==Inf]=igraph::vcount(seed_graph)
  sp=sp1[,index2]
  not2=colSums(sp>=meandist)
  not1=rowSums(sp>=meandist)
  maxbc=sum(sp>0 & sp<meandist)
  sp1mat=sp1 %x% rep(1,length(index2))
  sp2mat=rep(1,length(index1)) %x% sp2
  m=matrix(as.vector(sp2[,index1]),nrow=length(index1)*length(index2),ncol=igraph::vcount(seed_graph),byrow=F)
  m[m>=meandist]=0
  sumsp=sp1mat+sp2mat-m
  betweencount=colSums(sumsp==0)
  #for (i in 1:length(index1)){
  #  m=sp1[i,index2]
  #  betweencount[index1[i]]=betweencount[index1[i]]+sum(m>=meandist)
  #  betweencount[index2[m>=meandist]]=betweencount[index2[m>=meandist]]+1
  #  m[m>=meandist]=0
  #  
  #  mmat=matrix(m,nrow=length(index2),ncol=igraph::vcount(seed_graph),byrow=F)
  #  sp1mat=matrix(sp1[i,],nrow=length(index2),ncol=igraph::vcount(seed_graph),byrow=TRUE)
  #  sumsp=sp1mat+sp2-mmat
  #  betweencount=betweencount+colSums(sumsp==0)
  #  
  #}
  betweencount[index1]=betweencount[index1]-length(index2)+not1
  betweencount[index2]=betweencount[index2]-length(index1)+not2
  betweencount[intersect(index1,index2)]=betweencount[intersect(index1,index2)]+1
  betweencount=betweencount/maxbc
  list(allcount=betweencount,maxS2B=maxbc)
}



s2bthreshold=function(s2bvec){
  n=length(s2bvec)
  x=max(s2bvec)*(1:n)/n
  y=sort(s2bvec)
  distvec=sqrt((x-max(s2bvec))^2+y^2)
  s2bt=y[which.min(distvec)]
  s2bt
}

load("apid_main.Rdata")
symbnodes=read.table("symbnodes.txt",header=T,stringsAsFactors = F)
meandist=3.857
#genbet=estimate_betweenness(apid_main,vids=V(apid_main), directed=FALSE, 3)
load("genbet.RData")

# Define UI for random distribution app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("S2B - Specific Specific Betweenness"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select a file ----
      fileInput("file1", "Protein list A",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
  
      # Input: Select a file ----
      fileInput("file2", "Protein list B",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      # Input: Numeric entry for number of obs to view ----
      numericInput(inputId = "nrep1",
                   label = "Number of network edge randomizations:",
                   value = 0, min = 0),
      
      # Input: Numeric entry for number of obs to view ----
      numericInput(inputId = "nrep2",
                   label = "Number of seed randomizations:",
                   value = 0, min = 0),
      
      # Input: Slider for the number of observations to generate ----
      sliderInput("spec1",
                  "Edge specificity filter:",
                  value = 0,
                  min = 0,
                  max = 1),
      
      # Input: Slider for the number of observations to generate ----
      sliderInput("spec2",
                  "Seed specificity filter:",
                  value = 0,
                  min = 0,
                  max = 1),
      
      # Input: Slider for the number of observations to generate ----
      #sliderInput("resid",
      #            "Global betweenness deviation:",
      #            value = 0,
      #            min = -1,
      #            max = 1),
      
      # Input: Slider for the number of observations to generate ----
      sliderInput("n",
                  "Number of top candidates:",
                  value = 50,
                  min = 1,
                  max = 100)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Summary",
                           tableOutput("summary"),
                           plotOutput("s2bthres"),
                           br(),
                           plotOutput("s2b_gbt")
                           ),
                  
                  tabPanel("Table", 
                           
                           # Output: HTML table with requested number of observations ----
                           tableOutput("view"),
                           htmlOutput("linkDownload"),
                           br()
                           ),
                  
                  tabPanel("Network", 
                           #numericInput(inputId = "candn",
                            #            label = "S2B candidate:",
                            #            value = 1, min = 1, max=100),
                           htmlOutput("choose_cand"),
                           htmlOutput("out_network")
                            
                            )
      )
      
    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output) {
  
  # Reactive expression
  indexA <- reactive({
    req(input$file1)
    uniprotA=read.table(input$file1$datapath, header=T,stringsAsFactors = F)
    indexA=seedrows(apid_main,uniprotA[,1])
  })
  
  indexB <- reactive({
    req(input$file2)
    uniprotB=read.table(input$file2$datapath, header=T,stringsAsFactors = F)
    indexB=seedrows(apid_main,uniprotB[,1])
  })
  
  indexA_I = reactive({
    indexA_I=setdiff(indexA(),indexB())
  })
  
  indexB_I = reactive({
    indexB_I=setdiff(indexB(),indexA())
  })

  s2b_output= reactive({
    s2b_output=S2B(apid_main,indexA_I(),indexB_I(),input$nrep2,input$nrep1,meandist,symbnodes,genbet)
  })
    
  #s2b_pre= reactive({
  #  s2b_pre=S2B(apid_main,indexA_I(),indexB_I(),0,0,meandist,symbnodes,genbet)
  #})
  
  #s2b_rep2= reactive({
  #  s2b_rep2=S2B(apid_main,indexA_I(),indexB_I(),input$nrep2,0,meandist,symbnodes,genbet)
  #})
  
  #s2b_rep1= reactive({
  #  s2b_rep1=S2B(apid_main,indexA_I(),indexB_I(),0,input$nrep1,meandist,symbnodes,genbet)
  #})
  

  #s2b_output=reactive({
  #  s2b_output=s2b_pre()
  #  s2b_output$s2btable$sscore=s2b_rep2()$s2btable$sscore
  #  s2b_output$s2btabel$escore=s2b_rep1()$s2btable$escore
  #})

  
  s2b_t=reactive({
    s2b_t=s2bthreshold(s2b_output()$s2btable$S2B)
  })
  
  s2b_tabs=reactive({
    s2b_tabs=s2b_t()*s2b_output()$maxS2B
  })
  
  s2b_table=reactive({
    s2b_table=s2b_output()$s2btable
    s2b_table=s2b_table[s2b_table$S2Brank<=input$n,]
    s2b_table=s2b_table[s2b_table$escore>=input$spec1,]
    s2b_table=s2b_table[s2b_table$sscore>=input$spec2,]
    s2b_table=s2b_table[order(s2b_table$S2Brank),]
  })
  
  output$choose_cand <- renderUI({
    req(s2b_table())
    selectInput("candn", "S2B candidate:", choices=as.character(s2b_table()$symbol), selected=s2b_table()$symbol[1])
  })
  
  output$out_network = renderUI({
    req(input$candn)
    visNetworkOutput("network")  
  })
  
  
  s2b_candnet=reactive({
    s2b_candnet=cand_net(apid_main,input$candn,s2b_output(),indexA_I(),indexB_I(),symbnodes)
  })
  

  s2b_sum=reactive({
    s2b_sum=data.frame(
      Property=c("Network","# Nodes", "# Edges", "Avg. Distance","# List A", "# List B", "Max S2B", "S2B threshold", "# above threshold"),
      Value=c("APID PPI", vcount(apid_main), ecount(apid_main), meandist, length(indexA()), length(indexB()), s2b_output()$maxS2B, s2b_t(), sum(s2b_output()$s2btable$S2B>s2b_t()) ),
      Unique=c("-", "-", "-", "-",length(indexA_I()), length(indexB_I()),"-","-","-" )
    )
  })
  
  # Show the first "n" observations ----
  output$view <- renderTable({
    s2b_table()
  })
  
  output$linkDownload=renderUI({
    req(s2b_table())
    downloadLink('downloadData', 'Download full S2B table')  
  })
  
  
  output$summary = renderTable({
    s2b_sum()
  })
  
  output$s2bthres=renderPlot({
    plot((1:length(symbnodes[,1])),sort(s2b_output()$s2btable$S2B),type="l",xlab="Nodes", ylab="S2B")
    points(c(1,length(symbnodes[,1])),c(s2b_t(),s2b_t()),type="l",col="red")},
    width = 400, height = 450
  )
  
  
  output$s2b_gbt=renderPlot(
    plot(log10(genbet),log10(s2b_output()$s2btable$S2B),xlab="log10(Global Betweenness)", ylab="log10(S2B)"),
    width = 400, height = 450
  )
  
  output$network = renderVisNetwork({
    req(s2b_candnet())
    nodesdf=s2b_candnet()$nodes
    edgesdf=s2b_candnet()$edges
    visNetwork(nodesdf, edgesdf)%>% 
      visInteraction(navigationButtons = TRUE)%>%
      visPhysics(stabilization = FALSE) %>%
      visGroups(groupname = "A", color = "#de2d26") %>%
      visGroups(groupname = "B", color = "#2c7fb8") %>%
      visGroups(groupname = "S2B", color = "#756bb1") %>%
      visGroups(groupname = "node", color = "grey") %>%
      visLegend()
  })
  

    output$downloadData <- downloadHandler(
      filename = "s2btable.txt",
      content = function(con) {
        write.table(s2b_output()$s2btable, con, sep="\t", row.names=F, quote=F)
      }
    )
      
  

  


  
}

# Create Shiny app ----
shinyApp(ui, server)