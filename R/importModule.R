#' input UI
#'
#' @return A shiny tagList object that contains the input UI components
#' @rdname INTERNAL_inputUI
#' @keywords internal
#'
#' @importFrom shiny fluidRow NS actionButton icon uiOutput helpText tags div
#' @importFrom shinydashboardPlus box
#' @importFrom htmltools tagList h2
#' @importFrom shinyBS bsTooltip


importUI <- function(id="import")
{
  fluidRow(
    column(width=12,
           tags$head(tags$script(HTML(
             "$(document).on('shiny:connected', function() {
                $('.btn-file').removeClass('btn-default').addClass('btn-primary');
              });"
           ))),
           div(
             list(
               tags$label("Software", `for`="software"),
               selectInput(NS(id, "software"), label = NULL,
                           choices = c("DIA-NN"      = "diann",
                                       "DIA-NN (protein groups)" = "diann_pg",
                                       "Spectronaut" = "spectronaut",
                                       "Max-Quant (peptides)" = "maxquant",
                                       "Max-Quant (protein groups)" = "maxquant_pg",
                                       "Other" = "other"),
                           selected = "diann"
                           ),
               
               helpText(id = "tooltip_select_software",
                        "Specify the software used for identification. DIA-NN is set as default.")
             )
           ),
           div(
             list(tags$label("Input file", `for`="pe"),
                  fileInput(inputId=NS(id,"pe"), label=NULL, multiple = FALSE, accept = NULL, width = NULL),
                    helpText(
                      id="tooltip_input",
                      "Select the input file."
                    )
             )
           ),
           div(
             list(
               tags$label("Input file preview"),
               DT::dataTableOutput(NS(id, "pePreview"),)
             )
           ),
           uiOutput(NS(id, "columnSelectors")),
           uiOutput(NS(id, "quantColSelector")),
           uiOutput(NS(id, "nameInput")),
           div(
             list(
               tags$label("QFeatures", `for`="buildQFeatures"),
               actionButton(NS(id, "buildQFeatures"), "Build the QFeatures object", class = "btn-success"),
               helpText(id = "tooltip_buildQFeatures",
                        "The button generate a QFeatures object.")
             )
           ),
           div(
             list(
               tags$label("Export sample names for the annotation file", `for`="printed_annot"),
               downloadButton(NS(id, "printed_annot"), "Export sample names", class = "btn-primary"),
               helpText(id = "tooltip_printed_annot",
                        "The button generate a csv file with the samples' names of the input file, which the user can modify to create an annotation file.")
             )
           ),
           div(
             list(
               tags$label("Annotation file", `for`="annot"),
               fileInput(inputId=NS(id,"annot"), label=NULL, multiple = FALSE, accept = NULL, width = NULL),
               actionButton(NS(id, "addAnnot"), "Add annotation", class = "btn-success"),
               helpText(id = "tooltip_annotation",
                        "Upload the annotation file and click 'Add annotation'. The first column must match the sample names exactly. Grouping variables (e.g. condition) should be written as strings (e.g. 'D2', 'D4') so they are treated as factors. Numeric columns are kept as numeric.")
             )
           ),
           div(
             list(
               tags$label("Annotation file preview"),
               DT::dataTableOutput(NS(id, "annotPreview"),)
             )
           ),
           div(
             list(
               tags$label("QFeatures summary"),
               verbatimTextOutput(NS(id, "qfeaturesSummary"))
             )
           ),
           
          div(
            list(
              br(),
              h4("How do I cite MSqRob?"),
              h4("msqrob2 is free to use and  open source.
            When making use of msqrob2, we would appreciate it if you could cite our two papers."),
              h4("1. Sticker A, Goeminne L, Martens L, Clement L (2020). Robust Summarization and Inference in Proteome-wide Label-free Quantification. Molecular & Cellular Proteomics, 19(7), 1209-1219. doi: 10.1074/mcp.ra119.001624"),
              h4("2. Goeminne L, Gevaert K, Clement L (2016). Peptide-level Robust Ridge Regression Improves Estimation, Sensitivity, and Specificity in Data-dependent Quantitative Label-free Shotgun Proteomics. Molecular & Cellular Proteomics, 15(2), 657-668. doi: 10.1074/mcp.m115.055897")
              )
            )
    ) #end column
  )
}

#' Server for input tab
#'
#' @param id module id
#' @param variables global reactive values object to share objects across modules
#' @return list of reactive inputs
#' @rdname INTERNAL_inputServer
#' @keywords internal
#'
#' @importFrom shiny moduleServer updateSelectInput observeEvent eventReactive is.reactive showNotification removeNotification renderUI
#' @importFrom shinymeta metaReactive
#' @importFrom MultiAssayExperiment getWithColData
#' @importFrom DT datatable
#' @importFrom arrow read_parquet
#'
importServer <- function(id="import", variables){
  moduleServer(
    id,
    function(input,output,session){
      #extract datapath of inputfile
      peDatapath <- metaReactive({getDataPath(..(input$pe$datapath))})

      #import the input file
      observeEvent({input$pe},{
        peOut <- try(arrow::read_parquet(file=peDatapath()))
        if (inherits(peOut, "try-error")) {
          peOut <- try(data.table::fread(file=peDatapath(), check.names = TRUE, integer64 = "double"))

          if (inherits(peOut, "try-error")){
            showNotification("Upload proper input file",id="noProperFile",type="error",duration=NULL,closeButton=FALSE)
          }

        } else {
          removeNotification(id="File Not Valid")}
        variables$pe <- peOut
        variables$rawFilePath <- input$pe$datapath
        variables$rawFileName <- input$pe$name
      })
      
      
      
      # Track the software selection each of the two blocks below last rendered
      # for, so a software change can force its defaults even when the old pick
      # happens to also be a valid column name in the (possibly unchanged) file.
      prevSoftwareForColumns <- reactiveVal(NULL)
      prevSoftwareForQuant   <- reactiveVal(NULL)

      # select columns (feature ID + run column only — must NOT depend on
      # input$runCol/input$fnames themselves, or picking a value here would
      # immediately rebuild and reset this very block)
      output$columnSelectors <- renderUI({
        req(variables$pe)

        cols <- colnames(variables$pe)

        cfg <- switch(input$software,
                      diann = list(fnames = "Precursor.Id",   runCol = "Run"),
                      diann_pg = list(fnames = "Protein.Group", runCol = NULL),
                      spectronaut = list(fnames = "EG_PrecursorId", runCol = "R_FileName"),
                      maxquant = list(fnames = "Sequence",     runCol = NULL),
                      maxquant_pg = list(fnames = "Protein.IDs", runCol = NULL),
                      other = list(fnames = NULL, runCol = NULL)
        )

        # Preserve the user's current picks across re-renders (e.g. a new file
        # for the same software) as long as they're still valid columns in this
        # file; but force the software's own default whenever the software
        # selection itself just changed, even if the old pick happens to also
        # be a valid column name in this file. isolate() so that reading
        # input$fnames/input$runCol here does not itself retrigger this block.
        softwareChanged <- isolate(!identical(prevSoftwareForColumns(), input$software))

        fnamesSelected <- isolate({
          if (!softwareChanged && !is.null(input$fnames) && input$fnames %in% cols) {
            input$fnames
          } else {
            cfg$fnames
          }
        })
        runColSelected <- isolate({
          if (!softwareChanged && !is.null(input$runCol) && input$runCol %in% c("None", cols)) {
            input$runCol
          } else if (!is.null(cfg$runCol)) {
            cfg$runCol
          } else {
            "None"
          }
        })

        # Only treat the software change as "consumed" once its default was
        # actually realisable against this file's columns (i.e. is among the
        # choices being rendered below). If the software was switched before
        # the matching file was uploaded, the default can't be selected yet
        # (the widget would silently fall back to its first choice) — so keep
        # retrying on the next re-render instead of losing track of the change.
        appliedCleanly <- isolate({
          (is.null(cfg$fnames) || cfg$fnames %in% cols) &&
            (is.null(cfg$runCol) || cfg$runCol %in% cols)
        })
        if (!softwareChanged || appliedCleanly) isolate(prevSoftwareForColumns(input$software))

        list(
          tags$label("Column selection"),
          selectizeInput(NS(id, "fnames"), "Feature ID column",
                         choices = cols, selected = fnamesSelected),
          selectizeInput(NS(id, "runCol"), "Run column",
                         choices = c("None", cols), selected = runColSelected)
        )
      })

      # intensity column(s): long (quantCol) vs wide (quantCols) depends on
      # whether a run column is selected, kept separate from columnSelectors
      # so toggling it doesn't rebuild fnames/runCol
      output$quantColSelector <- renderUI({
        req(variables$pe)

        cols <- colnames(variables$pe)

        cfg <- switch(input$software,
                      diann = list(quantCol = "Precursor.Quantity"),
                      diann_pg = list(quantCol = setdiff(cols, c("Protein.Group", "Protein.Ids",
                                                                 "Protein.Names", "Genes",
                                                                 "First.Protein.Description",
                                                                 "N.Sequences",
                                                                 "N.Proteotypic.Sequences"))),
                      spectronaut = list(quantCol = "FG_MS2RawQuantity"),
                      maxquant = list(quantCol = grep("Intensity.", cols, value = TRUE)),
                      maxquant_pg = list(quantCol = grep("^LFQ.intensity", cols, value = TRUE)),
                      other = list(quantCol = NULL)
        )

        softwareChanged <- isolate(!identical(prevSoftwareForQuant(), input$software))

        isLong <- !is.null(input$runCol) && input$runCol != "None"

        # As above: only consume the software change once cfg$quantCol is
        # actually realisable against this file's columns; otherwise retry on
        # the next re-render (e.g. once the matching file is uploaded).
        appliedCleanly <- isolate(is.null(cfg$quantCol) || all(cfg$quantCol %in% cols))
        if (!softwareChanged || appliedCleanly) isolate(prevSoftwareForQuant(input$software))

        if (isLong) {
          selectizeInput(NS(id, "quantCol"), "Intensity column",
                         choices = cols,
                         selected = if (!softwareChanged && !is.null(input$quantCol) && input$quantCol %in% cols) input$quantCol else cfg$quantCol)
        } else {
          selectizeInput(NS(id, "quantCols"), "Intensity columns",
                         choices = cols,
                         selected = if (!softwareChanged && !is.null(input$quantCols) && all(input$quantCols %in% cols)) input$quantCols else cfg$quantCol,
                         multiple = TRUE)
        }
      })
      
      # preview raw input file
      output$pePreview <- DT::renderDataTable({
        req(variables$pe)
      
        DT::datatable(
          head(variables$pe, 5),  
          options = list(scrollX = TRUE)
        )
      })
      
      
      #name input summarized experiment, only showing for wide-format imports
      output$nameInput <- renderUI({
        if (input$software %in% c("maxquant", "other", "diann_pg", "maxquant_pg")) {
          div(list(
            tags$label("Set name"),
            textInput(NS(id, "name"), label = NULL, value = "quants_initial"),
            helpText("Assign a name to the summarized experiment.")
          ))
        }
      })
      
      
      # build qfeatures
      qfeatures <- eventReactive(input$buildQFeatures, {
        req(variables$pe)
        req(input$fnames)
        
        
        isLong <- !is.null(input$runCol) && input$runCol != "None"

        pe <- if (isLong) {
          req(input$runCol, input$quantCol)
          QFeatures::readQFeatures(
            assayData = variables$pe,
            fnames    = input$fnames,
            runCol    = input$runCol,
            quantCol  = input$quantCol,
            name      = input$name
          )
        } else {
          req(input$quantCols)
          QFeatures::readQFeatures(
            assayData = variables$pe,
            fnames    = input$fnames,
            quantCols = which(colnames(variables$pe) %in% input$quantCols),
            name      = input$name
          )
        }
        
        pe
       
      })
      
      # store in variables
      observe({
        req(qfeatures())
        variables$qfeatures        <- qfeatures()
        variables$qfeatures_import <- qfeatures()
        variables$software <- input$software
      })
      
      #extract datapath of annotation file
      annotDatapath <- metaReactive({getDataPath(..(input$annot$datapath))})
      
      #import the annotation file
      observeEvent({input$annot},{
        annotOut <- try(data.table::fread(file=annotDatapath(),check.names = TRUE, integer64 = "double",stringsAsFactors = TRUE))
        if (inherits(annotOut, "try-error")) {
          showNotification("Upload proper annotation file",id="noProperAnnotFile",type="error",duration=NULL,closeButton=FALSE)
        }
        else {
          removeNotification(id="File not valid. Rownames must match the sample names of the QFeatures object")}
        variables$annot <- annotOut
        variables$annotFilePath <- input$annot$datapath
        variables$annotFileName <- input$annot$name
      }
      )
      
      
      observe({
        req(variables$qfeatures)
        
        variables$annot_tmp <- data.frame(
        sampleName = rownames(colData(variables$qfeatures))
      )}
      )
      
      # write in csv the runCol with sample annotation
      output$printed_annot <- downloadHandler(
        filename = "annotation.tsv",
        content  = function(file) {
          
          write.table(variables$annot_tmp, file, sep = ";", row.names = FALSE)  
          
        },
        contentType = "text/csv"
      )
      
      # add annotation table to the QFeatures
      observeEvent(input$addAnnot, {
        req(variables$qfeatures)
        req(variables$annot)

        annot       <- as.data.frame(variables$annot)
        annotSamples <- annot[, 1]
        qf           <- variables$qfeatures
        sampleNames  <- rownames(colData(qf))

        missing <- setdiff(sampleNames, annotSamples)
        if (length(missing) > 0) {
          showNotification(
            paste("Annotation does not match sample names. Missing:", paste(missing, collapse = ", ")),
            type = "error", duration = NULL
          )
          return()
        }

        rownames(annot) <- annot[, 1]
        annot <- annot[sampleNames, , drop = FALSE]

        # convert character columns to factors for modelling
        annot <- as.data.frame(lapply(annot, function(col) {
          if (is.character(col)) as.factor(col) else col
        }))
        rownames(annot) <- sampleNames

        SummarizedExperiment::colData(qf) <- S4Vectors::DataFrame(annot)

        variables$qfeatures        <- qf
        variables$qfeatures_import <- qf
        showNotification("Annotation added successfully", type = "message")
      })
      
      # preview annotation input file
      output$annotPreview <- DT::renderDataTable({
        req(variables$qfeatures)
        req(variables$annot)
        
        DT::datatable(
          head(variables$annot, 5),  
          options = list(scrollX = TRUE)
        )
      })
      
      
      
      output$qfeaturesSummary <- renderPrint({
        req(variables$qfeatures)
        variables$qfeatures
      })
      
      return(
        list(
          qfeatures = reactive(variables$qfeatures),
          software  = reactive(variables$software),
          fnames    = reactive(input$fnames),
          runCol    = reactive(input$runCol),
          quantCol  = reactive(input$quantCol),
          quantCols = reactive(input$quantCols),
          name      = reactive(input$name)
        )
      )
    }
  )
}

