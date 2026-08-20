#' Check and install cmdstanr and CmdStan
#'
#' Same pattern as sccomp: ensure `cmdstanr` is available and that CmdStan
#' itself is installed. Used when `backend = "cmdstanr"`.
#'
#' @importFrom instantiate stan_cmdstan_exists
#' @noRd
check_and_install_cmdstanr <- function() {
  rlang::check_installed(
    pkg = "cmdstanr",
    reason = paste(
      "The {cmdstanr} package is required to run Stan models with",
      "backend = \"cmdstanr\". Install it with",
      "install.packages(\"cmdstanr\",",
      "repos = c(\"https://stan-dev.r-universe.dev\", getOption(\"repos\")))."
    ),
    action = function(...) {
      utils::install.packages(
        ...,
        repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org")
      )
    }
  )

  if (!stan_cmdstan_exists()) {
    stop(
      "CmdStan is required to proceed.\n\n",
      "Install it with:\n",
      "  cmdstanr::check_cmdstan_toolchain(fix = TRUE)\n",
      "  cmdstanr::install_cmdstan()\n",
      "See https://mc-stan.org/users/interfaces/cmdstan",
      call. = FALSE
    )
  }
  invisible(NULL)
}
