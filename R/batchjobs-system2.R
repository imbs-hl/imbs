
#' Wrapper around typically \code{batchMap} setup
#'
#' @title update_system2_job
#' @param reg.name            [\code{string}]\cr
#'                            Name of registry to load or create.
#' @param reg.dir             [\code{string}]\cr
#'                            Path where files regarding the registry / jobs should be saved.
#' @param work.dir            [\code{string}]\cr
#'                            Working directory for R process when experiment is executed.
#' @param sys.cmd             [\code{string}]\cr
#'                            Name of script to run through \code{\link{system2}}.
#' @param vec.cmd.args        [\code{character}] 
#'                            Vector to loop over.
#' @param more.cmd.args       [\code{list}]\cr
#'                            Static arguments, handed to more.args of \code{\link[BatchJobs]{batchMap}}.
#' @param seed                [\code{int}]\cr
#'                            Start seed for experiments. The first experiment in the registry will use this seed, for the subsequent ones the seed is incremented by 1.
#' @param cpus                [\code{int}]\cr
#'                            Number of CPUs to use for each job.
#' @param memory              [\code{int}]\cr 
#'                            Memory to use for each job in MB. Default is 5000.
#' @param walltime            [\code{int}]\cr
#'                            Maximum runtime for each job in minutes. Default is 59.
#' @param partition           [\code{string}]\cr
#'                            The partition the jobs are sent to.
#' @param env                 [\code{character}]\cr 
#'                            Directly passed to  \code{\link{system2}}. Character vector of name=value strings to set environment variables.
#' @param update              [\code{flag}]\cr
#'                            \code{TRUE} deletes the registry to force re-execution of jobs, default is \code{FALSE}, i.e. to load the registry and resume jobs not done.
#'
#' @return Nothing, throws an error if not all jobs are finished 
#' @export
#' @import BBmisc checkmate
#' @importFrom BatchJobs makeRegistry
#' @importFrom BatchJobs batchMap
#' @importFrom BatchJobs findNotDone
#' @importFrom BatchJobs submitJobs
#' @importFrom BatchJobs waitForJobs
update_system2_job <- function(reg.name,
                               reg.dir,
                               work.dir,
                               sys.cmd,
                               vec.cmd.args,
                               more.cmd.args,
                               seed,
                               cpus = 1,
                               memory = 5000,
                               walltime = 60,
                               partition = "fast",
                               env = character(),
                               update = FALSE){
  
  assertions <- checkmate::makeAssertCollection()
  
  assertString(reg.name, add = assertions)
  assertPathForOutput(reg.dir, overwrite = TRUE, add = assertions)
  assertPathForOutput(work.dir, overwrite = TRUE, add = assertions)
  assertCommand(sys.cmd, add = assertions)
  assertCharacter(vec.cmd.args, add = assertions)
  assertList(more.cmd.args, add = assertions)
  assertInt(seed, lower = 1, add = assertions)
  assertInt(cpus, lower = 1, add = assertions)
  assertInt(memory, lower = 1000, add = assertions)
  assertChoice(partition, choices = c("fast", "batch", "prio"), add = assertions)
  assertCharacter(env, add = assertions)
  assertFlag(update, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  
  # Delete registry if update
  if (update) {
    reg_abs <- file.path(reg.dir, reg.name)
    if ( file.exists(reg_abs)){
      unlink(reg_abs)
    }
  }
  
  # create or load reg
  reg <- BatchJobs::makeRegistry(id = reg.name,
                      file.dir = file.path(reg.dir, reg.name),
                      work.dir = work.dir,
                      #      src.files = batchjobs_args_file,
                      seed = seed)
  
  # Add jobs to map, if reg is empty
  if (length(findJobs(reg))) {
    ids <- findJobs(reg)
  } else {
    # build job map
    message('Build job map')
    BatchJobs::batchMap(reg,
             fun = function(cmd, vec_args, single_args){
               system2(command = cmd,
                       args = c(vec_args,
                                single_args),
                       env = env)
             },
             vec.cmd.args,
             more.args = list(cmd = sys.cmd,
                              single_args = more.cmd.args)
    )
  }
  
  # submit unfinished jobs, i.e. for first run: all
  unfinished_jobs <- BatchJobs::findNotDone(reg)
  if (length(unfinished_jobs) > 0){
    message(length(unfinished_jobs), ' jobs found, (re)submitting')
    BatchJobs::submitJobs(reg,
               ids = BBmisc::chunk(BatchJobs::findNotDone(reg), 
                                   n.chunks = 1),
               resources = list(ntasks = 1, 
                                ncpus = cpus, 
                                memory = memory,
                                walltime = walltime,
                                partition = partition),
               chunks.as.arrayjobs = TRUE)
  }
  # wait for jobs to finish
  wait <- BatchJobs::waitForJobs(reg)
  if (!wait){
    stop('Jobs for registry ', reg.name, 'not completed')
    
  }
}
