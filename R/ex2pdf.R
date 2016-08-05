#' Generation of Exercise Sheets and Exams in PDF Format
#' 
#' @description Automatic generation of exercise sheets and exams in PDF format,
#' inspired by the \code{exams} package.
#'
#' @param sheets      [\code{list}]\cr
#'                    A list with the following mandatory entries:
#'                    \itemize{
#'                      \item{\code{exercises}}{ [\code{list}] with file names 
#'                        to look for in \code{poolDir}. This is done 
#'                        recursively. Don't give the absolute path or relative 
#'                        path, just the file name. You can omit the extension. 
#'                        If a list entry is a \code{character} vector of file 
#'                        names, they get sampled according to \code{nsamp}.}
#'                      \item{\code{nsamp}}{ [\code{integer}] giving the number
#'                        of exercises to sample, if the corresponding entry of
#'                        \code{exercises} is a \code{character}. As most as 
#'                        long as \code{exercises}. Gets recycled if shorter.}
#'                    }
#'                    Additional entries are available during brewing of the
#'                    \code{templates} files. If the list entries are named, the
#'                    names are appended to \code{names} to generate unique file
#'                    names.
#' @param poolDir     [\code{string}]\cr
#'                    Path of the directory storing the exercise files.
#' @param templates   [\code{character}]\cr
#'                    Path of template files. Each file must follow the 
#'                    \code{\link[brew]{brew}} specifications. Each template will be 
#'                    brewed, knited and compiled into one PDF.
#' @param outDir      [\code{string}]\cr
#'                    Path of the output directory. Default is the current 
#'                    working directory.
#' @param names       [\code{character}]\cr
#'                    Specifies the naming scheme of each template in the output
#'                    directory. 
#' @param n           [\code{integer}]\cr
#'                    Specifies the number of copies to be generated per 
#'                    template. If \code{n}>1, the numer of copy is appended to 
#'                    the output file name.
#' @param compiler    [\code{list}]\cr
#'                    Specifies the LaTeX compiler and its arguments. Must have
#'                    the following entries:
#'                    \itemize{
#'                      \item{\code{command}}{ [\code{string}] name of the 
#'                      executable of the LaTeX compiler to be called.}
#'                      \item{\code{args}}{ [\code{character}] vector of 
#'                      arguments passed to the LaTeX compiler.}
#'                    }
#' @param clean       [\code{character} or \code{flag}]\cr
#'                    If \code{TRUE} all files but PDF files in \code{outDir} 
#'                    are removed after compilation. If a \code{character} 
#'                    vector, only files with the given extensions are removed
#'                    after compilation. If \code{FALSE} no clean up is done.
#' @param ...         Any further arguments passed to the environment of 
#'                    \code{\link[brew]{brew}}. This is helpful to pass static
#'                    elements such as e-mail adresses or semester identifier to
#'                    all \code{\link[brew]{brew}} calls.
#'
#' @return Invisible \code{character} with path of generated PDF files.
#' @export
#' 
#' @author Damian Gola, \email{gola@@imbs.uni-luebeck.de}
#'
#' @examples
#' \dontrun{
#' if(require(ggplot2)) {
#' sheets <- list(w1 = list(exercises = list("sum1", "sum2"), 
#'                          nsamp = c(1, 1),
#'                          date = "yyyy-mm-dd"),
#'                w2 = list(exercises = list(c("boxplot1", "regression1")),  
#'                          nsamp = 1,
#'                          date = "yyyy-mm-dd"))
#'
#' ex2pdf(sheets = sheets,
#'        poolDir = system.file("pool", package = "imbs"),
#'        templates = c(system.file("templates/exercise.brew", package = "imbs"),
#'                      system.file("templates/solution.brew", package = "imbs")),
#'        outDir = file.path(getwd(), "ex2pdfTest"),
#'        names = c("Exercise", "Solution"),
#'        n = 2,
#'        compiler = list(command = "pdflatex",
#'                        args = c("-interaction=batchmode")),
#'        year = 2016,
#'        email = "mail@example.net")
#'  }
#' }
ex2pdf <- function(sheets,
                   poolDir,
                   templates,
                   outDir = getwd(), names, n = 1,
                   compiler = list(command = "latexmk",
                                   args = c("-pdf", "-interaction=batchmode")),
                   clean = TRUE,
                   ...) {
  
  checkmate::assertList(sheets)
  
  if(all(sapply(sheets, is.list))) {
    # check if sheets is a list of valid sheet definitions
    invisible(sapply(sheets, checkSheet))
  } else {
    # check if sheets is a valid sheet definition
    checkSheet(sheets)
    sheets <- list(sheets)
  }
  
  checkmate::assertDirectory(poolDir)
  checkmate::assertFile(templates)
  if(!checkmate::testDirectory(outDir)) {
    warning("Directory '", outDir, "' does not exist. Creating now...")
    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  }
  checkmate::assertDirectory(outDir)
  checkmate::assertCharacter(names, len = length(templates))
  checkmate::assertList(compiler, "character", min.len = 2)
  checkmate::assertSubset(names(compiler), c("command", "args"))
  if(checkmate::testLogical(clean)) {
    checkmate::assertFlag(clean)
    cleanExt <- character()
  } else {
    checkmate::assertCharacter(clean, any.missing = FALSE, null.ok = FALSE, 
                               min.len = 1, min.chars = 1)
    cleanExt <- clean
    clean <- TRUE
  }
  
  # save the current working directory
  wd <- getwd()
  
  # get additional parameters for the template
  brewEnv <- list(...)
  
  # iterate over all sheets
  pdfs <- sapply(seq_along(sheets), function(nr, ...) {
    sheet <- sheets[[nr]]
    brewEnv <- c(brewEnv, sheet)
    
    # iterate over all versions
    sapply(1:n, function(version, ...) {
      # sample the exercises
      brewEnv$exercises <- unlist(mapply(ex = sheet$exercises, size = sheet$nsamp,
                                         FUN = function(ex, size) {
                                           sapply(sample(x = ex, size = min(length(ex), size)),
                                                  function(file) {
                                                    list.files(poolDir,
                                                               pattern = paste0("^", file),
                                                               full.names = TRUE,
                                                               recursive = TRUE)
                                                  })
                                         }))
      
      # create an environment for brewing
      brewEnv <- list2env(brewEnv)
      
      # brew the templates
      rnws <- sapply(templates, function(tmp, ...) {
        # change the working directory
        setwd(dirname(tmp))
        
        # create rnw name
        rnw <- sub(pattern = tools::file_ext(tmp),
                   replacement = "Rnw", x = tmp)
        
        # brew
        brew::brew(tmp,
                   output = rnw,
                   ...)
        rnw
      }, envir = brewEnv)
      
      # knit the brewed Rnw files
      setwd(outDir)
      texs <- mapply(rnw = rnws, name = names,
                     FUN = function(rnw, name) {
                       # build the file name
                       tex <- tools::file_path_sans_ext(name)
                       if(!is.null(names(sheets)[[nr]])) {
                         # sheet is not named
                         tex <- c(tex, names(sheets)[[nr]])
                       }
                       if(n > 1) {
                         # only one version requested
                         tex <- c(tex, version)
                       }
                       tex <- file.path(outDir,
                                        paste(paste(tex, collapse = "_"),
                                              "tex", sep = "."))
                       
                       # ensure that the output directory exists
                       dir.create(dirname(tex),
                                  showWarnings = FALSE,
                                  recursive = TRUE)
                       
                       # knit
                       knitr::knit(rnw, output = tex)
                     })
      
      # compile the tex files to get pdfs
      sapply(texs, FUN = function(tex) {
        message("Compiling ", tex, "...")
        sysOut <- system2(command = compiler$command,
                          args = c(compiler$args, tex),
                          stdout = TRUE,
                          stderr = TRUE)
        
        if(!is.null(attr(sysOut, "status"))) {
          # error during compilation
          log <- sub(pattern = tools::file_ext(tex),
                     replacement = "log", x = tex)
          # extract errors from log file and throw an error
          stop(paste(c(grep("error", sysOut, ignore.case = TRUE, value = TRUE),
                       paste("See", log, "for more information.")),
                     collapse = "\n"))
        }
        message("Done!")
      })
      
      # clean up
      if(clean) {
        rmFiles <- if(length(cleanExt)) {
          # find all files with specified file extension
          pattern <- paste0("[", 
                 paste(basename(tools::file_path_sans_ext(texs)), 
                       collapse = "|"), 
                 "]\\.[", paste(cleanExt, collapse = "|"), "]*$")
          list.files(path = outDir, 
                     pattern = pattern,
                     full.names = TRUE,
                     recursive = TRUE,
                     ignore.case = TRUE)
        } else {
          # find all files but PDF files
          grep(pattern = ".*(?<!pdf)$", 
               list.files(path = outDir, 
                          full.names = TRUE,
                          recursive = TRUE),
               perl = TRUE,
               value = TRUE)
        }
        unlink(rmFiles)
      }
      
      paste(tools::file_path_sans_ext(texs), "pdf", sep = ".")
    }, ...)
    
  }, brewEnv, cleanExt)
  
  names(pdfs) <- names(sheets)
  
  # set the working directory back
  setwd(wd)
  
  invisible(pdfs)
  
}

checkSheet <- function(sheet) {
  checkmate::assertList(sheet, min.len = 2, names = "named")
  checkmate::assertSubset(c("exercises", "date", "nsamp"), choices = names(sheet))
  checkmate::assertList(sheet$exercises, types = "character", min.len = 1)
  checkmate::assertIntegerish(sheet$nsamp, max.len = length(sheet$exercises))
}