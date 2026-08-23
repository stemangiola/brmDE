#' Lazy evaluation notice for unevaluated pipelines
#'
#' A `brmDE` pipeline is a targets graph. Building it (`brmDE()`,
#' `estimate()`, `hypothesis()`, `adjust()`) does not run it. Printing it, or
#' calling [tt_evaluate()], does. Assignment never prints, so a user who
#' writes `pipeline <- se |> brmDE() |> estimate()` would otherwise see
#' nothing and think the call failed.
#'
#' The notice cannot be printed at construction time. A pipe is one
#' expression, and mid-pipe more verbs may still be coming. It also cannot
#' be printed unconditionally after the expression: auto-print runs
#' `print.brmDE_hpc()` *before* top-level task callbacks, and that print
#' already evaluates. Saying "ready to be evaluated" after a result table
#' would be a lie.
#'
#' So the notice is deferred with [addTaskCallback()], and a small table of
#' stores records which pipelines were built and which were then run. After
#' the expression, if any store is still waiting, one message is printed.
#'
#' The two objects below are session state, not helpers, which is why they
#' sit at the top of this file: every function here reads or writes them,
#' and they must exist for the life of the package, not the life of a call.
#'
#' @name lazy_evaluation_message
#' @keywords internal
#' @noRd
NULL

#' Name of the top-level task callback that prints the ready notice.
#'
#' Registration is idempotent: `schedule_pipeline_ready_notice()` skips a
#' second callback of this name, and the callback returns `FALSE` so R drops
#' it after one expression.
#'
#' @keywords internal
#' @noRd
PIPELINE_READY_CALLBACK <- "brmDE-pipeline-ready"

#' Stores of pipelines built but not yet run.
#'
#' Keys are [brmDE()] `store` paths, one per pipeline, so evaluating one
#' object does not silence another left waiting in the same expression.
#' `parent = emptyenv()` keeps lookups from falling through to the session.
#'
#' @keywords internal
#' @noRd
pipeline_notice <- new.env(parent = emptyenv())

#' Tell the user the pipeline graph is waiting.
#'
#' @return `NULL`, invisibly, after [message()].
#'
#' @keywords internal
#' @noRd
pipeline_ready_message <- function() {
  message(
    "brmDE says: The pipeline is ready to be evaluated lazily: print the object, ",
    "or call tt_evaluate(), to run it."
  )
}

#' Whether `store` can be an environment name in `pipeline_notice`.
#'
#' `as_brmde_hpc()` reads `x$initialisation$store`, which is `NULL` on a
#' hand-built object. `pipeline_notice[[NULL]]` would error, and the notice
#' must not break a pipeline.
#'
#' @param store Candidate store path.
#'
#' @return `TRUE` for a single non-empty, non-missing character string.
#'
#' @keywords internal
#' @noRd
is_store_key <- function(store) {
  length(store) == 1L && is.character(store) && !is.na(store) && nzchar(store)
}

#' Record that a pipeline at `store` has been built and not yet run.
#'
#' @param store Store directory from [brmDE()].
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
pipeline_pending_add <- function(store) {
  if (is_store_key(store)) {
    pipeline_notice[[store]] <- TRUE
  }
  invisible()
}

#' Record that the pipeline at `store` has already been evaluated.
#'
#' Called from `tt_evaluate.brmDE_hpc()` so auto-print, which evaluates
#' before the task callback, does not then say the object is still waiting.
#'
#' @param store Store directory from [brmDE()].
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
mark_pipeline_evaluated <- function(store) {
  if (is_store_key(store) && exists(store, envir = pipeline_notice, inherits = FALSE)) {
    rm(list = store, envir = pipeline_notice)
  }
  invisible()
}

#' Whether any pipeline is still waiting, then forget them all.
#'
#' One top-level expression reports at most once. `all.names = TRUE` so a
#' store path that begins with a dot is not skipped by [ls()].
#'
#' @return `TRUE` if at least one store was pending.
#'
#' @keywords internal
#' @noRd
pipeline_pending_flush <- function() {
  stores <- ls(pipeline_notice, all.names = TRUE)
  rm(list = stores, envir = pipeline_notice)
  length(stores) > 0L
}

#' Defer the ready notice until the current top-level expression finishes.
#'
#' Assignment does not print, so nothing would otherwise say the graph is
#' waiting. A pipe of `brmDE() |> estimate() |> hypothesis()` is one
#' expression: the notice cannot fire mid-pipe, and auto-print runs the
#' pipeline before [addTaskCallback()] callbacks. Non-interactive sessions
#' skip this entirely, so batch runs, tests, and knitr never accumulate keys
#' with no callback to flush them.
#'
#' @param store Store directory from [brmDE()].
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
schedule_pipeline_ready_notice <- function(store) {
  if (!interactive()) {
    return(invisible())
  }
  pipeline_pending_add(store)
  if (PIPELINE_READY_CALLBACK %in% getTaskCallbackNames()) {
    return(invisible())
  }
  addTaskCallback(
    function(...) {
      if (pipeline_pending_flush()) {
        pipeline_ready_message()
      }
      FALSE
    },
    name = PIPELINE_READY_CALLBACK
  )
  invisible()
}
