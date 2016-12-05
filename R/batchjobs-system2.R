#' Wrapper around typically batchMap() setup
#'
#' @title update_system2_job
#' @param reg_name [char] Name of registry to load
#' @param sys_cmd  [char] Name of script to run through \code{\link{system2()}}
#' @param vec_cmd_args [char] Vector to loop over
#' @param more_cmd_args [char] Static arguments, handed to more.args of \code{\link{BatchJobs::batchMap}}
#' @param update [logical] TRUE deletes the registry to force re-execution of jobs, default is FALSE
#'
#' @return Nothing, throws an error if not all jobs are finished 
#' @export
update_system2_job <- function(reg_name, 
                               sys_cmd,
                               vec_cmd_args,
                               more_cmd_args,
                               update = FALSE){
  library(BatchJobs)
  # Delete registry if update
  if (update) {
    reg_abs <- file.path(reg_dir, reg_name)
    if ( file.exists(reg_abs)){
      unlink(reg_abs)
    }
  }
  
  # create or load reg
  reg <- makeRegistry(id = reg_name,
                      file.dir = file.path(reg_dir, reg_name),
                      work.dir = work_dir,
                      #      src.files = batchjobs_args_file,
                      seed = seed)
  
  # Add jobs to map, if reg is empty
  if (length(findJobs(reg))) {
    ids <- findJobs(reg)
  } else {
    # build job map
    message('Build job map')
    batchMap(reg,
             fun = function(cmd, vec_args, single_args){
               system2(command = cmd,
                       args = c(vec_args,
                                single_args))
             },
             vec_cmd_args,
             more.args = list(cmd = sys_cmd,
                              single_args = more_cmd_args)
    )
  }
  
  # submit unfinished jobs, i.e. for first run: all
  unfinished_jobs <- findNotDone(reg)
  if (length(unfinished_jobs) > 0){
    message(length(unfinished_jobs), ' jobs found, (re)submitting')
    submitJobs(reg,
               ids = BBmisc::chunk(BatchJobs::findNotDone(reg), 
                                   n.chunks = 1),
               resources = list(ntasks = 1, 
                                ncpus = 1, 
                                memory = memory,
                                walltime = walltime,
                                partition = partition),
               chunks.as.arrayjobs = TRUE)
  }
  # wait for jobs to finish
  wait <- waitForJobs(reg)
  if (!wait){
    stop('Jobs for registry ', reg_name, 'not completed')
    
  }
}
