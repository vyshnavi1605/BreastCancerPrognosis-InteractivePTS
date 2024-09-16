library(shiny)
library(shinythemes)
library(ggfortify)
library(igraph)
library(flexclust)
library(shinyjs)
library(seqHMM)

dat = read.csv('/Users/vyshnavidoddi/Documents/R Studio/FYP/data/breastcancerwisconsin.csv', header = TRUE,
               col.names= c("Id","ClumpThickness"," CellSize"," CellShape","MarginalAdhesion","EpithelialCellSize","BareNuclei","BlandChromatin","NormalNucleoli","Mitoses","Class"), sep=",")
dat = dat[2:11]

# Define UI for app that draws a histogram ----
ui <- fluidPage(theme = shinytheme("superhero"),
                navbarPage("Cancer Prognosis!",
                           tabPanel("Home",
                                    tags$img(src='dna.png', style = "width: 100%; height: 600px;"),
                                    fluidRow(
                                      style = "margin-top: 20px; margin-bottom: 20px;",
                                      column(width = 4,
                                             style = "background-color: #2C74B3; height: 200px; border: 1px solid #DDDDDD;",
                                             h3("ABOUT"),
                                             p("This application will help you predict cancer prognosis. It explains how
                                               to analyse the models that are created with artificial intelligence
                                               to predict prognosis. The models used in this application are: Cluster plots,
                                               Pseudo Time Series and Hidden Markov Models.")
                                      ),
                                      column(width = 4,
                                             style = "background-color: #205295; height: 200px; border: 1px solid #CCCCCC;",
                                             h3("DATA"),
                                             p("This application uses the breast cancer dataset which is open source and can be
                                               found on the UCI Machine Learning Repository. The dataset consists of 8 variables about
                                               the cells where the values range from 1-10, and a class variable ranging from 2 to 4,
                                               where 2 suggests healthy patients and 4 suggests terminal stage.")
                                      ),
                                      column(width = 4,
                                             style = "background-color: #144272; height: 200px; border: 1px solid #BBBBBB;",
                                             h3("HOW TO USE"),
                                             p("To get the best use out of this application, go through all the tabs and read the descriptions
                                               provided and make use of the interactivity provided.")
                                      )
                                    )
                           ),
                           tabPanel("Hidden Markov Model", h4("Hidden Markov Model"),
                                    p("Hidden Markov Model can be used to analyse time series data as it allows for pattern and trend
                                      identification in complex datasets. As they are able to capture the underlying processes that generate 
                                      the data, they aid in decision making and predictions."),
                                    p("Pseduo time series points are extracted to set up a hidden markov model. Below 
                                      shows 10 different models. The red line through the blue points shows a time series of
                                      the disease."),
                                    actionButton("prev", "Previous plot"),
                                    actionButton("next_plot", "Next Plot"),
                                    plotOutput("pca_plot")
                           ),
                           tabPanel("Clustering", h4("Clusters"), h5("Here is an exmaple of the clusters created on the breast cancer dataset"),
                                    p("Clustering is a useful tool for exploring and analysing data as it aids the discovery of patterns and 
                                      relationships within complex datasets. This is done by grouping similar data points together into clusters, 
                                      based on their features or attributes."),
                                    fluidRow(
                                      column(width = 8,
                                             tags$div(
                                               tags$img(src='clusters.png', width = 900,height = 400),
                                               tags$img(src = "clustertable.png", width = 900,height = 300),
                                             ),
                                      ),
                                      column(width = 4,
                                             style = "background-color: #125D98; height: 525px; border: 1px solid #BBBBBB;",
                                             h3("Explanation"),
                                             tags$ul(
                                               tags$li("Clusters 4 and 8 indicate healthy patients, with class column value of 2 for both suggesting good health.", style = "margin-bottom:10px;"),
                                               tags$li("Clusters 5 and 3 are in the very early stages of cancer, with class values slightly higher than 2 and higher 
                                                       clump thickness and cell size compared to clusters 4 and 8.", style = "margin-bottom:10px;"),
                                               tags$li("Cluster 3 is closely connected to cluster 2, while cluster 5 has faint to no connections to other clusters, 
                                                       which is also seen in the clump thickness distribution (3.27 for cluster 5, 4.53 for cluster 3), indicating 
                                                       early stages of cancer that can progress.", style = "margin-bottom:10px;"),
                                               tags$li("Cluster 2 diverges into two clusters, 1 and 7, with similar class values but different cell sizes, 
                                                       suggesting different treatments for cancer present in both cases.", style = "margin-bottom:10px;"),
                                               tags$li("Clusters 6 and 10 diverge from cluster 1, indicating different treatment requirements for cancer 
                                                       in these clusters as well.", style = "margin-bottom:10px;"),
                                               tags$li("Cluster 9 seems to indicate a terminal stage, with most values being the highest and being the 
                                                       last cluster on the cluster plot, and a class value of 4.", style = "margin-bottom:10px;"),
                                               tags$li("The important clusters can be identified by the diverging stages.", style = "margin-bottom:10px;")
                                             ),
                                      ),
                                    ),
                                    h3("Examples"),
                                    sidebarLayout(
                                      sidebarPanel(
                                        numericInput("obs", "Enter the number of the cluster to see its breakdown", value = 1, min = 1, max = 10),
                                        verbatimTextOutput("value"),
                                        actionButton("regenerate", "Recluster", class = "btn-info"),
                                      ),
                                      mainPanel(
                                        plotOutput(outputId = "clusterPlot"),
                                        tableOutput("cluster")
                                      )
                                      
                                    ),
                           ),
                           tabPanel("Pseduo time series", h3("Pseduo time series"),
                                    p("A pseduo time series plot, plots trajectories through 
                                       data based on disease stage this creates a realistic disease progression. 
                                       They are a useful tool for modeling the temporal structure 
                                       of data that was not collected over time."), 
                                    p("We can analyse the trajectories created and really understand the disease progression. And identify the 
                                      the important trajectories by identifying the time points at which significant changes or transitions 
                                      occur in the underlying data."),
                                    p("To alter the number of repitions chage the slider accordingly and click submit, click recluster
                                      to recluster the clusters and plot them along with the pseduo time series."),
                                    # Sidebar layout with input and output definitions ----
                                    sidebarLayout(
                                      
                                      # Sidebar panel for inputs ----
                                      sidebarPanel(
                                        tags$h3("PTS attributes"),
                                        # Input: Slider for the number of bins ----
                                        sliderInput("repititions", "Repititions:",
                                                    min = 0, max = 100,
                                                    value = 100),
                                        actionButton("submitbutton", "Submit", class = "btn-info"),
                                        actionButton("recluster", "Recluster", class = "btn-info"),
                                      ),
                                      
                                      
                                      
                                      # Main panel for displaying outputs ----
                                      mainPanel(
                                        plotOutput(outputId = "distPlotPTS"),
                                      )
                                    )),
                           tabPanel("Bar plot", h4("Bar Plot"), p("Here is a cluster plot, using this plot below
                                                                     select the clusters depending on the path of the disease 
                                                                     you want to see the bar plot for."), 
                                    plotOutput(outputId = "clustPlot"),br(), 
                                    sidebarLayout(
                                      sidebarPanel(
                                        selectInput("column_select", "Select a column", choices = names(dat)),
                                        numericInput("number", "Enter number between 1 and 10:", min = 1, max = 10, value = 1),
                                        actionButton("add", "Add"),
                                        actionButton("clear", "Clear"),
                                        verbatimTextOutput("selected_numbers"),
                                      ),
                                      mainPanel(
                                        p("Bar plots can help identify patterns and trends in the data, as well as
                                          highlight differences and similarities between groups. The bar plot below can be 
                                          used to visualise the clusters for each feature in the dataset."),
                                        p("Using the cluster plot provided above, identify a potential course of the disease,
                                          and input the relevant cluster numbers into the field and press add. You can change the features
                                          of that you want to visualise using the dropdown provided."),
                                        plotOutput(outputId = "boxPlot")
                                      )),
                                    
                           ),
                           tabPanel("Density graph", h4("Density Graph"), plotOutput(outputId = "densityPlot"),
                                    p("The plot above shows us the distribution of the cell size over time. We can see
                                      hoe each clusters distribution changes over time, i.e as the cancer progresses.")
                           )
                ),
)

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  dat = read.csv('/Users/vyshnavidoddi/Documents/R Studio/FYP/data/breastcancerwisconsin.csv', header = TRUE,
                 col.names= c("Id","ClumpThickness"," CellSize"," CellShape","MarginalAdhesion","EpithelialCellSize","BareNuclei","BlandChromatin","NormalNucleoli","Mitoses","Class"), sep=",")
  dat$BareNuclei <- as.numeric(as.character(dat$BareNuclei))
  dat[ dat == "?" ] <- NA
  dat = na.omit(dat) #deletion of missing data
  data = dat[c(2,4:11)]
  clusterdata = dat[2:11]
  
  
  fullclass = dat[,3]#vector of classes where 0 is control
  class = 1+(fullclass>1)#converts multiclasses into 2: 1 = control, 2 = others
  
  Kfit <- kmeans(clusterdata, 10) # 3 cluster solution
  df = Kfit$centers
  
  plot_number <- reactiveVal(1) # initialize plot number to 1
  
  
  observeEvent(input$next_plot, {
    plot_number(plot_number() + 1)
  })
  
  observeEvent(input$prev, {
    plot_number(plot_number() - 1)
  })
  
  selectedValues1 <- reactive({
    as.numeric(input$checkGroup)
  })
  
  output$selectedValues <- renderText({
    paste("Selected Values:", paste(selectedValues1(), collapse = ", "))
  })
  
  
  selected_numbers <- reactiveVal(c())
  
  observeEvent(input$add, {
    if (input$number %in% 1:10 && !(input$number %in% selected_numbers())) {
      selected_numbers(c(selected_numbers(), input$number))
    }
  })
  
  observeEvent(input$clear, {
    selected_numbers(c())
  })
  
  output$selected_numbers <- renderPrint({
    selected_numbers()
  })
  
  obs <- reactive({
    validate(need(input$obs >= 1 & input$obs <= 10, "Input must be between 1 and 10."))
    input$obs
  })
  
  output$value <- renderText({
    obs()
  })
  
  
  banfrom = c(3,4);
  banto = c(4,3);
  
  #build full dist matrix
  dis <- dist(data,method="euclidean")
  dis=as.matrix(dis)
  
  sampsize=nrow(data);
  cons=dis;
  for (i in 1:sampsize)
  {
    for (j in 1:sampsize)
    {
      for (k in 1:length(banfrom))
      {
        if ((fullclass[i]==banfrom[k]) & (fullclass[j]==banto[k]))
        {
          cons[i,j]=999;
          cons[j,i]=999;
        }
      }
    }
  }
  
  pcadat=(pc.cl  <- princomp(data))
  
  clusters <- cclust(x = pcadat$scores[,1:2], k = 10, method = "kmeans", dist = "euclidean", simple = FALSE, save.data = TRUE)
  Kclusters <- clusters@cluster
  plot(clusters)
  
  
  distribution <- reactiveVal(
    cclust(x = pcadat$scores[,1:2], k = 10, method = "kmeans", dist = "euclidean", simple = FALSE, save.data = TRUE)
  )
  
  observeEvent(input$regenerate, {
    distribution(cclust(x = pcadat$scores[,1:2], k = 10, method = "kmeans", dist = "euclidean", simple = FALSE, save.data = TRUE))
  })
  
  cluster_means_reactive <- reactive({
    # get the cluster assignment for each observation
    cluster_assignment <- predict(distribution())
    # get the size of each cluster
    cluster_sizes <- table(cluster_assignment)
    # get the mean values for each variable in each cluster
    cluster_means <- aggregate(clusterdata, by = list(cluster_assignment), mean)
    # subset the cluster_means based on the selected observation number
    return(cluster_means[obs(), ])
  })
  
  #****PTS****
  sampsize = 50;#determines the length of the pts
  observeEvent(input$repititions, {
    nreps<<-input$repititions
  })
  endclass=2;#endclass is represented by 2 and startclass by 1
  
  pts = NULL
  for (i in 1:nreps)
  {
    #resample from data
    repeat
    {
      mysamp = sample(1:nrow(data), sampsize, replace=TRUE)
      if (is.element(1,class[mysamp]) & is.element(endclass,class[mysamp]))
      {
        break
      }
    }    
    dunesamp = data[mysamp,]
    classsamp = class(mysamp)
    
    repeat
    {
      startp = sample(1:sampsize,1)
      endp = sample(1:sampsize,1)
      if ((class[mysamp[startp]]==1) && (class[mysamp[endp]]==endclass))
      {
        break
      }
    }
    
    dissamp = cons[mysamp,mysamp]
    
    
    #build minimum spanning tree
    mode(dissamp) <- "numeric"
    g <- graph.adjacency(dissamp,weighted=TRUE,mode = "undirected")
    
    datamst=minimum.spanning.tree(g)#, weights = TRUE, algorithm = NULL) ALSO mst
    E(datamst)$weight
    
    #get shortest path from class 1 to class 2
    datshort=get.shortest.paths(datamst, startp,endp,mode="all",output = "both")
    
    #plot the trajectories of each mst
    ptsind = datshort$vpath
    pts[i] = list(mysamp[ptsind[[1]]])
  }
  
  mylist = list()
  for(i in 1:10){
    datatmp <- data[pts[[i]],]
    list_seq <- seqdef(datatmp)
    mylist <- append(mylist, list(list_seq))
  }
  
  set.seed(1010)
  init_hmm_bf2 <- build_hmm(
    observations = mylist,
    n_states = 10, left_right = TRUE, diag_c = 2)
  
  
  sim <- simulate_hmm(
    10,
    init_hmm_bf2$initial_probs,
    init_hmm_bf2$transition_probs,
    init_hmm_bf2$emission_probs,
    9)
  
  
  pcadat=(pc.cl  <- princomp(data))
  
  
  ###ggplots
  FULLPTS=data.frame()
  for (i in 1:nreps)
  {
    t=1:length(pts[[i]])
    c=fullclass[pts[[i]]]
    pcas=pcadat$scores[pts[[i]],1:2]
    IDPTS=cbind(t,c,pcas)
    
    FULLPTS=rbind(FULLPTS,IDPTS)
  }
  
  
  output$distPlotPTS <- renderPlot({
    
    if(input$recluster>0){
      Kfit <- kmeans(clusterdata, 10) # 3 cluster solution
      Kclusters = Kfit$cluster
      
      plot(pcadat$scores[,1],pcadat$scores[,2],col=Kclusters, pch=19)
      
      #impose plot of first n pts on pca plot
      for (i in 1:nreps)
      {
        lines(pcadat$scores[pts[[i]],1],pcadat$scores[pts[[i]],2], col=class[pts[[i]]])
      }
    }
    
    if(input$submitbutton>0){
      plot(pcadat$scores[,1],pcadat$scores[,2],col=Kclusters, pch=19)
      
      #impose plot of first n pts on pca plot
      for (i in 1:nreps)
      {
        lines(pcadat$scores[pts[[i]],1],pcadat$scores[pts[[i]],2], col=class[pts[[i]]])
      }
    }
    else{
      return("Server is ready for calculation.")
    }
  })

  output$clusterPlot <- renderPlot({
    plot(distribution())
  })
  
  output$cluster <- renderTable(
    {
      return(cluster_means_reactive())
    }
  )
  
  output$clusterData <- renderTable(
    {
      X <- Kfit$centers  
      means <- colMeans(df)
      sds <- apply(df,2,sd)
      return(means)
    }
  )
  
  output$pca_plot <- renderPlot({
    
    
    i <- plot_number()
    
    if (i < 1) {
      plot_number(1)
      i <- 1
    } else if (i > 10) {
      plot_number(10)
      i <- 10
    }
    
    
    simdatatmp <- sim$observations[[i]]
    csv_fname = "/Users/vyshnavidoddi/Documents/R Studio/FYP/data/timeseries.csv"
    write.table(simdatatmp, file = csv_fname, sep = ",",
                append = FALSE, quote = FALSE,
                col.names = FALSE, row.names = FALSE)
    timeseries <- read.csv("/Users/vyshnavidoddi/Documents/R Studio/FYP/data/timeseries.csv", header = TRUE,
                           col.names= c("ClumpThickness"," CellShape","MarginalAdhesion","EpithelialCellSize","BareNuclei","BlandChromatin","NormalNucleoli","Mitoses","Class"), sep=",")
    combined_data <- rbind(data, timeseries)
    pca_result <- prcomp(combined_data)
    pca_result=(pc.cl  <- princomp(combined_data))
    temp = rbind(matrix(1,682,1),matrix(2,10,1))
    plot(pca_result$scores[,1],pca_result$scores[,2],col=temp, pch=19)
    last_10_points = pca_result$scores[(dim(pca_result$scores)[1]-9):dim(pca_result$scores)[1],1]
    last_10_points_ = pca_result$scores[(dim(pca_result$scores)[1]-9):dim(pca_result$scores)[1],2]
    lines(last_10_points, last_10_points_, col="blue", lwd=2)
  })
  
  output$boxPlot <- renderPlot({
    
    clusterdata$cluster <- clusters@cluster
    # subset the data based on selected clusters
    subset_data <- subset(clusterdata, cluster %in% selected_numbers())
    
    # specify the order of levels based on selected clusters
    subset_data$cluster <- factor(subset_data$cluster, levels = selected_numbers())
    
    # plot boxplot
    ggplot(subset_data, aes(x = factor(cluster), y = !!as.symbol(input$column_select), fill = factor(cluster))) +
      geom_boxplot() +
      labs(x = "Cluster", y = input$column_select)
  })
  
  output$clustPlot <- renderPlot({
    plot(clusters)
  })
  
  output$densityPlot <- renderPlot({
    g <- ggplot(FULLPTS, aes(FULLPTS$t))
    g + geom_density(aes(fill=factor(FULLPTS$c)), alpha=0.8)
  })
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)