#Càrrega de llibreries 

library(shiny)
library(SRAdb)
library(DT)
library(shinycssloaders)
library(shinyjs)
library(dplyr)
library(ggplot2)

# Funció per a la descàrrega de dades i el control de qualitat
runQualityControl <- function(run_accession, desti) {
  
  desti <- "/home/miquela/RNASeqAnalysis/arxius_fastq_descarregats"
  
  # Indica el nom del fitxer per comprovar si ja existeix
  nom_fitxer <- paste0(run_accession, ".fastq")
  nom_fitxer_1 <- paste0(run_accession, "_1.fastq")
  nom_fitxer_2 <- paste0(run_accession, "_2.fastq")
  
  # Comprovem si el fitxer o fitxers ja existeix en el directori indicat
  if (!file.exists(file.path(desti, nom_fitxer)) || !file.exists(file.path(desti, nom_fitxer_1)) || !file.exists(file.path(desti, nom_fitxer_2))) {
    
    # Construïm la comanda per descarregar l'arxiu
    comanda <- paste("fasterq-dump --split-files --progress", run_accession, "-o", file.path(desti, nom_fitxer), sep=" ")
    
    # Executem la comanda amb system() i capturem la sortida
    system(comanda, intern=TRUE)
  }
  
  # Comprovar si ja s'ha executat abans el control de qualitat
  fitxer_resultats <- paste0(run_accession, "_fastqc.zip")
  fitxer_resultats_1 <- paste0(run_accession, "_1_fastqc.zip")
  fitxer_resultats_2 <- paste0(run_accession, "_2_fastqc.zip")
  
  
  # Comprovem si el fitxer o fitxers de control de qualitat ja existeix en el directori indicat
  if (!file.exists(file.path(desti, fitxer_resultats)) || !file.exists(file.path(desti, fitxer_resultats_1)) || !file.exists(file.path(desti, fitxer_resultats_2))) {
    
    # Construïm la comanda per fer el control de qualitat de l'arxiu
    comanda1 <- paste("fastqc", file.path(desti, nom_fitxer), "-o", desti, sep=" ")
    system(comanda1, intern=TRUE)
  }
  #Comanda per executar multiqc en el directori de desti
  comanda2 <- paste("python3 -m multiqc", desti, "-o", desti, "-f", sep=" ")
  system(comanda2, intern=TRUE)
}


server <- function(input, output) {
  
  # Definim la taula que es mostrarà automàticament quan s'obre l'aplicació a la pestanya de cercador d'estudis
  taula_consulta<-reactive({
    #if (input$demo){
     # sra_dbname <- file.path(system.file('extdata', package='SRAdb'), 'SRAmetadb_demo.sqlite')	
#      sra_con <- dbConnect(dbDriver("SQLite"), sra_dbname)
      
 #   }else{}
    
    sqlfile <- "/media/disc_sda4/gran_SRAmetadb.sqlite"
    #sqlfile <- "/home/miquela/RNASeqAnalysis/SRAmetadb.sqlite"
    
    sra_con <- dbConnect(SQLite(),sqlfile)
    
    #sra_dbname <- file.path(system.file('extdata', package='SRAdb'), 'SRAmetadb_demo.sqlite')	
    #sra_con <- dbConnect(dbDriver("SQLite"), sra_dbname)
    
    keyword <- input$keyword
    study_id <- input$studyID
    run_accession <- input$runAccession
    center_project_name <- input$specie
    run_date <- input$year
    
    # Comença la consulta
    consulta <- "SELECT study.study_id, study.study_accession, study.study_title, study.study_type, study.study_abstract, study.study_description, study.center_name, study.center_project_name, run.run_accession, run.run_date FROM study JOIN run ON study.submission_accession = run.submission_accession"
    
    # Creem una llista amb les condicions
    condicions <- list()
    
    if (study_id != "") {
      condicions <- c(condicions, paste0('study_ID = "', study_id, '"'))
    }
    
    if (keyword != "") {
      condicions <- c(condicions, paste0('study_abstract like "%', keyword, '%"'))
    }
    
    if (run_accession != "") {
      condicions <- c(condicions, paste0('run.run_accession = "', run_accession, '"'))
    }
    
    if (center_project_name != "") {
      condicions <- c(condicions, paste0('center_project_name like "%', center_project_name, '%"'))
    }
    
    if (run_date != "") {
      condicions <- c(condicions, paste0('run.run_date like "%', run_date, '%"'))
    }
    
    # Afegim les condicions a la consulta
    if (length(condicions) > 0) {
      consulta <- paste(consulta, "where", paste(condicions, collapse = " and "))
    }
    
    # Executem la consulta
    taula <- dbGetQuery(sra_con, consulta)
    taula
  })
  
  # Creem un histograma del nombre de runs per any
  output$run_histogram <- renderPlot({
    run_data <- taula_consulta()
    run_data$year <- format(as.Date(run_data$run_date), "%Y")
    
    run_counts <- table(run_data$year) #Crea una taula del nombre de runs per any
    
    ggplot(data.frame(year = names(run_counts), count = as.integer(run_counts)), aes(x = year, y = count)) +
      geom_col(fill = "#009081") +
      labs(title = "Run number per year", x = "Year", y = "Number of runs") +
      theme_minimal()
    })
  
  
  # Creem un histograma del nombre d'estudis per any
  output$study_histogram <- renderPlot({
    study_data <- taula_consulta()
    
    # Agrupem els experiments segons study_id
    study_data <- study_data %>%
      distinct(study_ID, .keep_all = TRUE)
    
    study_data$year <- format(as.Date(study_data$run_date), "%Y")
    
    study_counts <- study_data %>%
      group_by(year) %>%
      summarise(count = n(), .groups = 'drop')
    
    ggplot(study_counts, aes(x = year, y = count)) +
      geom_col(fill = "#009081") +
      labs(title = "Study number per year", x = "Year", y = "Number of studies") +
      theme_minimal()
  })
  
  
  # Creem una variable reactiva per emmagatzemar els estudis seleccionats
  selected_studies <- reactiveVal()
  
  output$taula <- DT::renderDT({
    datatable(
      data.frame(taula_consulta()),
      extensions = 'Buttons',
      options = list(
        dom = 'Blfrtip',
        buttons = list(
          list(
            extend = 'csv',
            text = 'Download as CSV'
          ),
          list(
            extend = 'excel',
            text = 'Download as Excel'
          ),
          list(
            extend = 'pdf',
            text = 'Download as PDF'
          )
        ),
        pageLength = 5,
        lengthMenu = list(c(5, 10, 15, 20, -1), c('5', '10', '15', '20', 'All')),
        # Afegim l'opció de selecció de múltiples files
        selection = 'multiple'
      )
    )
  })
  
  # Quan l'usuari selecciona una nova fila, actualitzem la variable reactiva
  observeEvent(input$taula_rows_selected, {
    selected_studies(taula_consulta()[input$taula_rows_selected, ])
  })
  
  output$selected_studies <- DT::renderDT({
    selected_data <- selected_studies()
    if (!is.null(selected_data)) {
      datatable(selected_data)
    }
  })
  
  observeEvent(input$runQC, {
    # Mostra l'spinner
    shinyjs::show("loading")
    
    selected_data <- selected_studies()
    
    if (!is.null(input$qcRunAccession) && input$qcRunAccession != "") {
      run_accession <- input$qcRunAccession
    } else if (!is.null(selected_data)) {
      run_accession <- selected_data$run_accession
    } else {
      # Amaga l'spinner si no hi ha cap run_accession introduït ni cap fila seleccionada
      shinyjs::hide("loading")
      return()
    }
    
    runQualityControl(run_accession, desti)
    
    output$qc_report <- renderUI({
      # Amaga l'spinner quan l'informe està llest per ser mostrat
      shinyjs::hide("loading")
      includeHTML("/home/miquela/RNASeqAnalysis/arxius_fastq_descarregats/multiqc_report.html")
    })
  })
}