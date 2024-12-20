


teReadLib <- function(libPath, libIdentifier = NULL){#libpath in the format ("./outDir/species/strain/lib.fa")
  if(file.exists(libPath)){
    cat("Loading TE library:", libPath, "\n")
    teLib <- Biostrings::readDNAStringSet(libPath)
    species <- libPath %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 3) %>% unlist()
    strain <- libPath %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 4) %>% unlist()
    teLib@metadata$species <- species
    teLib@metadata$strain <- strain
    teLib@ranges@NAMES <- gsub(pattern = "#DNA", replacement = "#TIR", teLib@ranges@NAMES)
    teLib@ranges@NAMES <- gsub(pattern = "#Unknown", replacement = "#Unclassified", teLib@ranges@NAMES)
    teLib@ranges@NAMES <- gsub(pattern = "?", replacement = "", teLib@ranges@NAMES, fixed = TRUE)
    teLib@metadata$libID <- libIdentifier
    if(length(teLib) > 0){#prevent error when there are no sequences in the fasta file
      repType <- teLib@ranges@NAMES %>%
        strsplit("#") %>% lapply("[[", 2) %>% unlist() %>% #separate header by the # and select the second part (eliminate rnd-X...)
        strsplit(" ") %>% lapply("[[", 1) %>% unlist() #separate the second part of the header by the space and take the repetitive family type
      teOrder <- repType %>% strsplit("/", fixed = TRUE) %>% lapply("[[", 1) %>% unlist()
      teType <- ifelse(!is.na(match(teOrder, teClassification$Order))
      , yes = teClassification$AppName[match(teOrder, teClassification$Order)]
      , no = repType)
      teLib@ranges@metadata$repetitiveType <- teType
      cat(libPath, " successfully loaded\n")
    } else {
      cat("There are no sequences in the selected library: ", libPath, " -------\n")
      teLib <- NULL
    }
  } else {
    cat(libPath, " does not exist\n")
    teLib <- NULL
  }
  teLib #return the library
}

tePlotLib <- function(teLibList, libType = "Please identify the library plot!!"){
  if(class(teLibList) == "list"){
    # cat("This is executed")
    inputDataComb <- lapply(teLibList, function(teLib){
      if(length(teLib) >= 1){
        inputData <- data.frame(
          element = teLib@ranges@metadata$repetitiveType
          , width = teLib@ranges@width
          , species = teLib@metadata$species, strain = teLib@metadata$strain
          , libID = teLib@metadata$libID
        )
      }else{
        cat("There are empty libraries", teLib, "\n")
        inputData <- data.frame(element = character(), width = numeric(), species = character(), strain = character(), libID = character())
      }
      inputData  
    }) %>% do.call(rbind, .)
    # cat("Data for plotting is generated")
    p <- inputDataComb %>%
      ggplot(aes(x = element)) + theme_classic() +
      geom_bar(aes(fill = libID), position = "dodge") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "inside", legend.position.inside = c(0.9,1)
            ,plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
      labs(x = "", y = "Count", title = paste0(unique(inputDataComb$species), "--", unique(inputDataComb$strain))) +
      scale_fill_discrete(name = "Library type")
    ggplotly(p) %>% layout(bargap = 0.3, legend = list(orientation = 'h', x = 0.1, y = 1))
  } else {
    teLib <- teLibs
    inputData <- data.frame(
      element = teLib@ranges@metadata$repetitiveType
      , width = teLib@ranges@width
      , species = teLib@metadata$species, strain = teLib@metadata$strain
    )
    p2 <- inputData %>%
      ggplot(aes(x = element)) + theme_classic() +
      geom_bar(aes(fill = element), position = "dodge", show.legend = FALSE) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "bottom"
            ,plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
      labs(x = "", y = "Count", title = paste0(unique(inputData$species), "--", unique(inputData$strain)), subtitle = libType)
    p2
  }
}

teConsDivergencePlot <- function(seqName, blastFiles = NULL, leftTrim = 1, rightTrim = 1){
  # teSeqId <- teSeq %>% strsplit("#", fixed = TRUE) %>% lapply("[[", 1) %>% unlist()
  blast <- tryCatch({
    blast = read.table(file = grep("trimming", x = blastFiles, value = TRUE, invert = TRUE)
                       , col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
    blast$cons_len <- max(blast$qend)
    blast$is_full_length <- blast$length >= FL_thresh*unique(blast$cons_len)
    nFlcNoTrim <- sum(blast$is_full_length)
    blast = read.table(file = grep("trimming", x = blastFiles, value = TRUE)
                       , col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
    blast$pTitle <- paste(seqName, "(trimmed)")
    blast$nFlcNoTrim <- nFlcNoTrim
    blast
  }, error = function(e){
    blast = read.table(file = grep("trimming", x = blastFiles, value = TRUE, invert = TRUE)
                       , col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
    blast$pTitle <- seqName
    blast
  })
  # blast = read.table(file = blastFiles
  #                    , col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) #https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
  blast$cons_len <- max(blast$qend)
  blast$is_full_length <- blast$length >= FL_thresh*unique(blast$cons_len)
  blast$full_length <- ifelse(blast$is_full_length, yes = "FLF", "No FLF")
  #subset blast noFLF to prevent slow app loading
  if(nrow(blast) > 2000) {
    tmpBlast <- blast %>% subset(full_length == "No FLF")
    subBlast <- tmpBlast[sample(1:nrow(tmpBlast), size = min(nrow(tmpBlast), 1000)),]
    tmpBlast <- blast %>% subset(full_length == "FLF")
    tmpBlast <- tmpBlast[sample(1:nrow(tmpBlast), size = min(nrow(tmpBlast), 1500)),]
    subBlast <- rbind(tmpBlast, subBlast)
  } else {
    subBlast <- blast
  }
  p1 <- subBlast %>% arrange(desc(full_length)) %>%
    ggplot() + theme_classic() +
    geom_segment(aes(x = qstart, xend = qend, y = 100 - pident, colour = full_length)) +
    geom_vline(xintercept = leftTrim, lty = 2) +
    geom_vline(xintercept = rightTrim, lty = 2) +
    # geom_label(aes(x = 1, y = 1, label = paste0("FLC no trim:", unique(nFlcNoTrim)))) +
    # annotate(x = 10, y = 1, geom = "label", label = paste0("FLC no trim: ", unique(nFlcNoTrim)))+
    labs(x = "TE consensus (bp)", y = "divergence to consensus (%)"
         , title = paste0("TE: ", unique(blast$pTitle))
         , subtitle = paste0("consensus size:", unique(blast$cons_len)
                             , "bp; fragments: ", nrow(blast)
                             , "; full length: ", sum(blast$is_full_length), "(>=", FL_thresh*unique(blast$cons_len), "bp)"
                             , "\n", paste0("FLC no trim: ", unique(nFlcNoTrim)))) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    scale_color_manual(values = c("FLF" = "red", "No FLF" = "grey"), name = "Full length")
  p1
}


teConsCoveragePlot <- function(seqName, blastFiles = NULL, leftTrim = 1, rightTrim = 1, source = "A"){
  # teSeq <- input$teManualSeq %>% strsplit("#", fixed = TRUE) %>% lapply("[[", 1) %>% unlist()
  coverage <- tryCatch({
    blastTrim <- read.table(file = grep("trimming", blastFiles, value = TRUE)
                       , col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
    trimPos <- grep("trimming", blastFiles, value = TRUE) %>%
      strsplit("trimming") %>% lapply("[[", 1) %>% unlist() %>%
      paste0(., "trimming/", seqName, ".fa_trimPositions.tsv") %>%
      read.table(header = TRUE, sep = "\t")
    # blastTrim$qstart <- blastTrim$qstart+trimPos$leftTrim-1
    # blastTrim$qend <- blastTrim$qend+trimPos$rightTrim-1
    blastTrim$cons_len <- max(blastTrim$qend)
    coverageTrim <- lapply(seq_along(blastTrim$qseqid), function(i){
      seqCoverage <- c(
        rep(FALSE, blastTrim$qstart[i]-1) #inicio de la secuencia sin cubrir (parte izquierda)
        , rep(TRUE,abs(blastTrim$qend[i]-blastTrim$qstart[i])+1) #secuencia cubierta
        , rep(FALSE, blastTrim$cons_len[i]-blastTrim$qend[i]) #secuencia final sin cubrir (parte derecha)
      )#location of each sequence in consensus sequence
    }) %>% do.call(rbind, .) %>% colSums() %>%
      data.frame(seqPosition = 1:length(.), Coverage = ., Trim = "Trimmed") # coverage of the consensus
    coverageTrim$seqPosition <- coverageTrim$seqPosition+trimPos$leftTrim-1
    leftTrim <- leftTrim + trimPos$leftTrim-1
    rightTrim <- rightTrim + trimPos$leftTrim-1
    
    blast = read.table(file = grep("trimming", blastFiles, value = TRUE, invert = TRUE)
                       , col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) #https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
    blast$cons_len <- max(blast$qend)
    coverage <- lapply(seq_along(blast$qseqid), function(i){
      seqCoverage <- c(
        rep(FALSE, blast$qstart[i]-1) #inicio de la secuencia sin cubrir (parte izquierda)
        , rep(TRUE,abs(blast$qend[i]-blast$qstart[i])+1) #secuencia cubierta
        , rep(FALSE, blast$cons_len[i]-blast$qend[i]) #secuencia final sin cubrir (parte derecha)
      )#location of each sequence in consensus sequence
    }) %>% do.call(rbind, .) %>% colSums() %>%
      data.frame(seqPosition = 1:length(.), Coverage = ., Trim = "Not trimmed") # coverage of the consensus
    rbind(coverage, coverageTrim) %>% cbind(pTitle = "\nNot trimmed and Trimmed")
  }, error = function(e){
    blast = read.table(file = grep("trimming", blastFiles, value = TRUE, invert = TRUE)
                       , col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) #https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
    blast$cons_len <- max(blast$qend)
    coverage <- lapply(seq_along(blast$qseqid), function(i){
      seqCoverage <- c(
        rep(FALSE, blast$qstart[i]-1) #inicio de la secuencia sin cubrir (parte izquierda)
        , rep(TRUE,abs(blast$qend[i]-blast$qstart[i])+1) #secuencia cubierta
        , rep(FALSE, blast$cons_len[i]-blast$qend[i]) #secuencia final sin cubrir (parte derecha)
      )#location of each sequence in consensus sequence
    }) %>% do.call(rbind, .) %>% colSums() %>%
      data.frame(seqPosition = 1:length(.), Coverage = ., Trim = "Not trimmed", pTitle = "") # coverage of the consensus
  })
  
  # blast$is_full_length <- blast$length >= FL_thresh*unique(blast$cons_len)
  # blast$full_length <- ifelse(blast$is_full_length, yes = "FLF", "No FLF")
  # coverage <- lapply(seq_along(blast$qseqid), function(i){
  #   seqCoverage <- c(
  #     rep(FALSE, blast$qstart[i]-1) #inicio de la secuencia sin cubrir (parte izquierda)
  #     , rep(TRUE,abs(blast$qend[i]-blast$qstart[i])+1) #secuencia cubierta
  #     , rep(FALSE, blast$cons_len[i]-blast$qend[i]) #secuencia final sin cubrir (parte derecha)
  #   )#location of each sequence in consensus sequence
  # }) %>% do.call(rbind, .) %>% colSums() %>%
  #   data.frame(seqPosition = 1:length(.), Coverage = .) # coverage of the consensus
  (covPlot <- coverage %>%
      ggplot(aes(x = seqPosition, y = Coverage)) + theme_classic() +
      geom_line(aes(colour = Trim)) +
      geom_vline(xintercept = leftTrim, lty = 2) +
      geom_vline(xintercept = rightTrim, lty = 2) +
      labs(x = "TE consensus (bp)", y = "Coverage (bp)", title = paste0("TE consensus genomic coverage", unique(coverage$pTitle)), subtitle = seqName) +
      scale_color_manual(name = "Seq type", values = c("Trimmed" = "black", "Not trimmed" = "lightgrey"))
  )
  p <- plotly::ggplotly(covPlot, source = source) %>%
    layout(legend = list(
      orientation = "h",    # Horizontal
      x = 0.5,              # Centrado horizontalmente
      y = -0.2,             # Debajo del gráfico
      xanchor = "center",   # Anclado al centro
      yanchor = "top"       # Anclado arriba
    ))
  #   # Registrar el evento para este grafico
  event_register(p, "plotly_click")
}

teConsSelfPlot <- function(teSeq, blastFiles = NULL, leftTrim = 1, rightTrim = 1){
  # teSeq <- input$teManualSeq %>% strsplit("#", fixed = TRUE) %>% lapply("[[", 1) %>% unlist()
  selfBlast <- tryCatch({
    blast = read.table(file = grep("trimming", blastFiles, value = TRUE)
                       , col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
    blast$pTitle <- " (trimmed)"
    blast
  }, error = function(e){
    blast = read.table(file = grep("trimming", blastFiles, value = TRUE, invert = TRUE)
                       , col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
    blast$pTitle <- ""
    blast
  })
  # selfBlast = read.table(file = blastFiles
  #                        , col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) #https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
  ###########################-
  ## dot-plot (bottom left)##
  ###########################-

  plot(x = 1, type = "n", xlim = c(0,selfBlast$qend[1]),ylim = c(0,selfBlast$qend[1]), col = "white",
       main = paste0("TE consensus self dotplot\n", unique(selfBlast$pTitle)),
       ylab = paste(as.character(selfBlast$qseqid[1]), "(bp)", sep = " "),
       xlab = paste(as.character(selfBlast$qseqid[1]), "(bp)", sep = " ")
  )
  abline(v = leftTrim, lty = 2)
  abline(v = rightTrim, lty = 2)
  for(i in 1:length(selfBlast$qseqid)){
    if(selfBlast$send[i] > selfBlast$sstart[i]){
      segments(x0 = selfBlast$qstart[i], x1 = selfBlast$qend[i], y0 = selfBlast$sstart[i], y1 = selfBlast$send[i], col = "black", lwd = 1.5)
    } else {
      segments(x0 = selfBlast$qstart[i], x1 = selfBlast$qend[i], y0 = selfBlast$sstart[i], y1 = selfBlast$send[i], col = "#009E73", lwd = 1.5)
    }
    # if orientation
  } # for each segment end
}

# seqName <- "rnd-1_family-20_unconfirmed"
# teSeq <- seqName
# blastFile <- fileNames$te_aidSelfBnFiles <- list.files(fileNames$teManualOutDir, pattern = glob2rx(paste0(seqName, ".fa*self.blastn*")), recursive = TRUE, full.names = TRUE) %>% tail(1)
teConsStructurePlot <- function(seqName, blastFiles = NULL, leftTrim = 1, rightTrim = 1){
  # teSeq <- input$teManualSeq %>% strsplit("#", fixed = TRUE) %>% lapply("[[", 1) %>% unlist()
  selfBlast <- tryCatch({
    blastFile <- grep("trimming", blastFiles, value = TRUE)
    blast = read.table(file = blastFile
                       , col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
    blast$pTitle <- " (trimmed)"
    blast
  }, error = function(e){
    blastFile <- grep("trimming", blastFiles, value = TRUE, invert = TRUE)
    blast = read.table(file = blastFile
                       , col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
    blast$pTitle <- ""
    blast
  })
  # selfBlast = read.table(file = blastFile
  #                        , col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) #https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
  selfBlast=selfBlast[order(selfBlast$qstart, decreasing = F),]
  #####################################-
  ## Annotation graph (bottom right) ##
  #####################################-
  
  orfFile <- gsub(pattern = ".fa.self.blastn.out", replacement = ".fa.orftetable"
                  , blastFile)
  
  # test if there are orf detected; will later store in orf if TRUE
  ORFs <- tryCatch({
    orfs <- read.table(orfFile, sep = "\t", fill = NA)
    
    if (orfs$V1[1] == seqName) {
      cat("ORF table generated by MCHelper and not by TE-aid, skipping first line ----\n")
      orfs <- read.table(orfFile, sep = "\t", skip = 1)
    } else {
      cat("ORF table generated by TE-aid ----\n")
      orfs <- read.table(orfFile)
    }
    
    orfs
  }, error = function(e) {
    # message("Error: ", e$message)
    data.frame() # Devolver un data.frame vacío en caso de error
  })
  

  ## Arrows layer ##

  plot(x = 1, type = "n", xlim = c(0,selfBlast$qend[1]),ylim = c(-max(nrow(ORFs),10),length(selfBlast$qseqid)), col = "white", yaxt="n",
       main = paste0("TE consensus structure and protein hits\n", unique(selfBlast$pTitle)),
       xlab = paste(as.character(selfBlast$qseqid[1]), "(bp)", sep = " "),
       ylab = ""
  )
  abline(v = leftTrim, lty = 2)
  abline(v = rightTrim, lty = 2)
  for(i in seq(1:length(selfBlast$qseqid))){
    if(selfBlast$sstart[i] > selfBlast$qstart[i]){
      arrows(x0 = selfBlast$qstart[i], x1 = selfBlast$qend[i], y0 = i, y1 = i, col = rainbow(length(selfBlast$qseqid))[i], lwd = 3, length = 0.1)
      arrows(x0 = selfBlast$sstart[i], x1 = selfBlast$send[i], y0 = i, y1 = i, col = rainbow(length(selfBlast$qseqid))[i], lwd = 3, length = 0.1)
    } # if to draw (filter)
  } # for each segment end

  ## Orfs layers ##
  
  
  tryCatch({
    if(nrow(ORFs) < 1){ # if orf table is empty
      text(paste("no orf >",os," bp detected", sep=""), x=selfBlast$length[1]/2, y=-5, cex = 2) #os is just the minimum length for ORF
    } else { # if ORF table, plot ORFs
      for(i in seq(1:length(ORFs$V1))){
        if(ORFs$V1[i] < ORFs$V2[i]){ # checking ORF orientation
          rect(xleft = ORFs$V1[i], xright = ORFs$V2[i], # draw a + ORF
               ybottom = -i-0.15, ytop = -i+0.15, lwd = 1, border = "black")
        } else {
          rect(xleft = ORFs$V1[i], xright = ORFs$V2[i], # draw a - ORF
               ybottom = -i-0.15, ytop = -i+0.15, lwd = 1, border = "red")
        } # orientation
        
        ## TE protein hits (blastp) ##
        rect(xleft = ORFs$V5[i], xright = ORFs$V6[i],
             ybottom = -i-0.15, ytop = -i+0.15, lwd = 1, col = as.character(paste("#",ORFs$V8[i], sep="")), border = "white") # draw colored rectangle same way as orf
        text(paste(ORFs$V3[i], ORFs$V4[i]), x = (min(ORFs$V1[i],ORFs$V2[i])+max(ORFs$V1[i],ORFs$V2[i]))/2, y = -i+0.15, pos = 3) # print hit name
        
      } # for each segment
      #names(ORFs)<-c("orf.start", "orf.end", "hit.TE.prot", "TE.Class", "hit.start", "hit.end", "strand", "color")
      #print(ORFs)
    } # if orf present plot orfs and prot hits
  }
  , error = function(e) {
    # text(paste("no orf >",os," bp detected", sep=""), x=selfBlast$length[1]/2, y=-5, cex = 2)#os is just the minimum length for ORF
    cat("no orf to plot ----\n")
  }
  # , finally = cat("Hola")
  )
}




create_progress_bar <- function(status) {
  steps <- c("Not started", "MCHelper automatic running", "MCHelper MI for curated sequences", "MI started", "MI finished", 
             "Non 80 MI started", "70 models", "Incomplete Models", "Step 9")
  
  current_step <- steps[status]
  
  # Crear un elemento HTML para la barra de progreso
  sprintf(
    '<div class="progress-container" style="width: 100%%; cursor: pointer;" title="%s" data-step="%s" onclick="Shiny.setInputValue(\'progress_clicked\', this.getAttribute(\'data-step\'), {priority: \'event\'})">
       <div class="progress-bar" style="width: %.0f%%; background-color: #007bff;">%d</div>
     </div>',
    current_step, current_step, (status / length(steps)) * 100, status
  )
}

getMchStatus <- function(species, strain, rawLib = "", outDir = "./MCHelper"){
  mchDir <- paste(outDir, species, strain, sep = "/")
  teRepeatLib <- teReadLib(paste("RM2_output", species, strain, rawLib, sep = "/"))
  
  if(!dir.exists(mchDir)) {
    dir.create(mchDir, recursive = TRUE);cat("Output folder created:", mchDir, "\n")
  }else{cat("Output dir already existing:", mchDir, "\n")}
  
  teCleanFile <- paste0(mchDir, "/", strain, "-clean_families.fa")
  #create clean families if they do not exist
  if(file.exists(teCleanFile)){
    cat("Clean lib already existing ---- \n")
  } else {
    ### Preparing input
    #we need to exclude some TE families as MCHelper does not manage them
    teCleanLib <- teRepeatLib[!grepl(pattern = "Satellite|Simple_repeat|tRNA|rRNA|Retroposon|snRNA|scRNA", x = names(teRepeatLib))]
    Biostrings::writeXStringSet(teCleanLib, filepath = teCleanFile)
  }
  
  status <- 1
  if(file.exists(paste0(mchDir, "/mchelper_automatic.log"))){
    status <- 2 #MCHelper automatic module has been executed
  } 
  if (file.exists(paste0(mchDir, "/mchelper_manual.log"))){
    logFile <- readLines(paste0(mchDir, "/mchelper_manual.log"))
    if(logFile[length(logFile)] == "MCHelper with TE-Aid successfully run"){
      #mchelper manual inspection finished, time to review manually
      status <- 3
    }
  }
  if (file.exists(paste0(mchDir, "/1_MI_MCH/tmp_curated_sequences_NR.fa"))) {
    #manual review started
    status <- 4
    tmpCuratedLib <- paste0(mchDir, "/1_MI_MCH/tmp_curated_sequences_NR.fa")
    if (length(Biostrings::readDNAStringSet(tmpCuratedLib)) == 0) {
      status <- 5
    }
  }
  if (file.exists(paste0(mchDir, "/1_MI_MCH/2_MI_MCH/tmp_non808080.fa"))){
    #non808080 manual review started
    status <- 6
    tmpNon80Lib <- paste0(mchDir, "/1_MI_MCH/2_MI_MCH/tmp_non808080.fa")
    if(length(Biostrings::readDNAStringSet(tmpNon80Lib)) == 0){
      status <- 7
    }
  }
  if (file.exists(paste0(mchDir, "/1_MI_MCH/2_MI_MCH/blast707070.log"))){
    status <- 8
  }
  if (file.exists(paste0(mchDir, "/1_MI_MCH/3_MI_MCH/tmp_incomplete_non808080.fa"))){
    #non808080 manual review started
    # status <- 8
    tmpIncNon80Lib <- paste0(mchDir, "/1_MI_MCH/3_MI_MCH/tmp_incomplete_non808080.fa")
    if(length(Biostrings::readDNAStringSet(tmpIncNon80Lib)) == 0){
      status <- 9
    }
  }
  status
}
# lapply(1:nrow(data), function(i) getMchStatus(data$Species[i], data$Strain[i], rawLib = data$Raw_lib[i], outDir = "./MCHelper")) %>% unlist()
# i <- 2
# getMchStatus(data$Species[i], data$Strain[i], rawLib = data$Raw_lib[i], outDir = "./MCHelper")
# 
# getMchStatus(species = "C.elegans", strain = "test")



# deNovoClass <- read.table("MCHelper/Anopheles/gambiae/1_MI_MCH/denovoLibTEs_PC.classif", sep = "\t", header = TRUE)
# seqName <- "rnd-6_family-30"
teClassOutput <- function(deNovoClass, seqName){
  teInfo <- deNovoClass %>% subset(Seq_name == seqName) %>% select(Seq_name, class, order, sFamily,  other, strand)
  
  teCoding <- tryCatch({
    teCoding <- deNovoClass %>% subset(Seq_name == seqName) %>% 
      pull(coding) %>% strsplit("profiles:", fixed = TRUE) %>% lapply("[[", 2) %>% unlist() %>%
      gsub(pattern = ")", replacement = "", x = .) %>%
      strsplit(",", fixed = TRUE) %>% unlist() %>% #trimws() %>%
      strsplit(":", fixed = TRUE) %>% lapply(., trimws) %>%
      do.call(rbind.data.frame, .)
    colnames(teCoding) <- c("Coding", "Evalue")
    # teCoding$Evalue <- as.numeric(teCoding$Evalue)
    teCoding$Evalue <- strsplit(teCoding$Evalue, split = "e") %>% 
      lapply(function(vec) paste(substr(vec[1], start = 1, stop = 2), vec[2], sep = "e")) %>%
      unlist()
    teCoding
  }
    , error = function(e) teCoding <- data.frame()
  )
  
  teStruc <- tryCatch({
    teStruc <- deNovoClass %>% subset(Seq_name == seqName) %>% 
      pull(struct) %>% strsplit("struct=(", fixed = TRUE) %>% lapply("[[", 2) %>% unlist() %>%
      gsub(pattern = ")", replacement = "", x = .) %>%
      strsplit(";", fixed = TRUE)%>% lapply(., trimws) %>%
      do.call(rbind.data.frame, .) %>% t()
    # colnames(teStruc) <- "Structure"
    # teStruc %>% select(Structure)
  }, error = function(e) teStruc <- data.frame())
  return(list(
    teInfo = teInfo
    , teCoding = teCoding
    , teStruc = teStruc
  ))
}

# teClassOutput(deNovoClass, seqName)

teCurationStats <- function(teLibList){
  libNames <- c("teRepeatLib" = "Raw lib", "teCleanLib" = "Clean lib"
                , "teAutoCuratedLib" = "Curated seqs", "teCompleteModelsLib" = "Complete models"
                , "teIncompleteModelsLib" = "Incomplete models"
                , "teAutoCuratedTmp" = "Remaining to MI"
                , "te80Lib" = "808080 lib", "teNon80TmpLib" = "Remaining to MI"
                , "te80RecovLib" = "808080 recovered"
                , "teStandbyLib" = "Stand by lib"
                , "te70FiltLib" = "Filtered 707070"
                , "teIncompleteRecovLib" = "Incomplete 808080 recovered"
                , "teIncompleteNon80Tmp" = "Remaining to MI"
                , "teFinalNR" = "Final lib"
                )
  
  teStats <- sapply(teLibList, function(lib){
    list(
      nSeqs = names(lib) %>% length()
      , nComplete = names(lib)[!grepl(names(lib), pattern = "_inc|_unconfirmed")] %>% length()
      , nIncomplete = names(lib)[grepl(names(lib), pattern = "_inc|_unconfirmed")] %>% length()
    ) %>% do.call(cbind, .)
  }, USE.NAMES = TRUE, simplify = FALSE) %>% do.call(rbind.data.frame, .)
  rownames(teStats) <- libNames[rownames(teStats)]
  teStats
}

# names(teLibs)

# teCurationStats(teLibs)

# lib1 <- teReadLib("./MCHelper/Drosophila/monieri/MCH_final/final_curated_NR.fa", libIdentifier = "Adrian")
# lib2 <- teReadLib("../D.monieri/Drosophila/monieri/final_curated_NR.fa", libIdentifier = "Marta")
# lib3 <- teReadLib("../D.monieri/Drosophila/monieri/D.monieri.TE_library.fasta", libIdentifier = "Simon")
# lib3 <- teLib
# lib3@metadata$libID <- "Simon"
# lib3@metadata$species <- "Drosophila"
# lib3@metadata$strain <- "monieri"
# lib3
# tePlotLib(list(lib3, lib2, lib1))
# 
# table(names(lib3) %in% names(lib2))
# table(names(lib1) %in% names(lib2))
# 
# # blastn -q seq[i]
# 
# adriBlast <- data.table::fread("../D.monieri/Drosophila/monieri/adri_vs_Simon.blast.out", sep = "\t")
# 
# unique(adriBlast$V1)
# 
# martaBlast <- data.table::fread("../D.monieri/Drosophila/monieri/marta_vs_Simon.blast.out")
# martaBlast$V1 %>% unique()
# 
# 
# setdiff( names(lib1), adriBlast$V1)
# 
# adriBlast2 <- data.table::fread("../D.monieri/Drosophila/monieri/adri_vs_Marta.blast.out", sep = "\t")
# 
# adriBlast2$V1 %>% unique()
# 
# setdiff(names(lib1), adriBlast2$V1) %>% unique()
# 
# intersect(names(lib3), names(lib2))
# setdiff(martaBlast$V1, adriBlast$V1) %>% unique()

# martaBlast$V1 %>% unique()
# martaBlast$V2 %>% unique() #recupera 116 de Simon
# 
# #secuencias de simon que marta recupera
# 
# #secuencias de simon que marta no recupera
# martaNotSimon <- lib3[!names(lib3) %in% unique(martaBlast$V2)] # las secuencias que Marta no recupera de simon
# 
# 
# adriBlast$V1 %>% unique()
# adriBlast$V2 %>% unique() # recupero 113 de Simon
# adriNotSimon <- lib3[!names(lib3) %in% unique(adriBlast$V2)] #estas son las 14 de simon que no me quedo yo
# 
# intersect(names(adriNotSimon),
# names(martaNotSimon))
# 
# 
# lib2[martaBlast$V1 %>% unique()] %>%
#   Biostrings::writeXStringSet(filepath = "../D.monieri/Drosophila/monieri/final_curated_NR_80Simon.fa")
# # makeblastdb -in final_curated_NR_80Simon.fa -dbtype "nucl"
# 
# cat("blastn -query final_curated_NR_Adri.fa -db final_curated_NR_80Simon.fa -outfmt 6 -qcov_hsp_perc 80 -perc_identity 80 -max_hsps 1 -out adri_vs_Marta.blast.out")
# 
# adriMalSeqs <- setdiff(martaBlast$V2,adriBlast$V2) #secuencias de Simon que marta recupera y yo no
# lib3[adriMalSeqs] %>% 
#   Biostrings::writeXStringSet(filepath = "../D.monieri/Drosophila/monieri/adriNotIncluded.fa")
# # conda run -n MCHelper python3 mchelper-ats/MCHelper.py 
# " conda run -n MCHelper --no-capture-output python3 ../../../TEammo_app-dev/mchelper-ats/MCHelper.py -r T --input_type fasta -l adriNotIncluded.fa -g ../../../TEammo_app_v2.1/0_raw/Drosophila/monieri/GCA_035047585.1.fna -o adriNotIncluded_MCH -t 5"
# adriMalSeqs2 <- martaBlast %>% subset(V2 %in% adriMalSeqs)
# adriMalSeqs2$V1 %in% names(lib1) #el gypsy yo también lo tengo pero no me ha hecho homología 80 con el de Simon
# lib2[adriMalSeqs2$V1] %>%
#   Biostrings::writeXStringSet(filepath = "../D.monieri/Drosophila/monieri/adriNotIncluded2.fa")
# 
# "conda run -n MCHelper --no-capture-output python3 ../../../TEammo_app-dev/mchelper-ats/MCHelper.py -r T --input_type fasta -l adriNotIncluded2.fa -g ../../../TEammo_app_v2.1/0_raw/Drosophila/monieri/GCA_035047585.1.fna -o adriNotIncluded2_MCH -t 5"
# martaMalSeqs <- setdiff(adriBlast$V2, martaBlast$V2) #secuencias de Simon que Adri recupera y Marta no
# martaMalSeqs2 <- adriBlast %>% subset(V2 %in% martaMalSeqs)
# 
# 
# 
# 
# blastFinal <- data.table::fread("HELIANO/C.elegans/N2/hel-NR_vs_final-curated-NR.blast")
# blastAll_80 <- data.table::fread("HELIANO/C.elegans/N2/hel-NR_vs_allDB-merged.blast")
# blastAll_Final <- data.table::fread("HELIANO/C.elegans/N2/hel-NR_vs_allDB-final.blast")
# blastAll <- data.table::fread("HELIANO/C.elegans/N2/hel-NR_vs_allDB.blast")
# ggvenn::ggvenn(list(
#   Final = blastFinal$V1
#   , all80 = blastAll_80$V1
#   , allFinal = blastAll_Final$V1
#   , all = blastAll$V1
#   )) +
#   labs(title = "Query Helitrons recovered in each database") + 
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# ggvenn::ggvenn(list(
#   Final = blastFinal$V2
#   , all80 = blastAll_80$V2
#   , allFinal = blastAll_Final$V2
#   , all = blastAll$V2
# )) +
#   labs(title = "Database Helitrons recovered in each database") + 
#   theme(plot.title = element_text(hjust = 0.5))
# 


















