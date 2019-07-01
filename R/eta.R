#' Functions related to eta parameter used in optim and kkt checks
#'
#' @description Used for gradient of eta. Currently being passed to optim in
#'   \code{\link{lmmlasso}} and used in \code{\link{kkt_check}}
#' @seealso \code{\link{logliklasso}}, \code{\link{kkt_check}}, \code{\link{lmmlasso}}
#' @inheritParams logliklasso
#' @inheritParams kkt_check
gr_eta_lasso_fullrank <- function(eta, sigma2, beta, eigenvalues, x, y, nt, myweights) {
  di <- 1 / myweights + eta * (eigenvalues - 1 / myweights)

  (1 / 2) * sum(((eigenvalues - 1) / di) * (1 - (((y - x %*% beta)^2) / (sigma2 * di))))
}

gr_e_lasso_fullrank <- function(e, sigma2, beta, eigenvalues, x, y, nt, myweights) {
  eta = e2eta(e)
  gr_eta_lasso_fullrank(eta, sigma2, beta, eigenvalues, x, y, nt, myweights) * (-eta)
}




#' @rdname gr_eta_lasso_fullrank
fn_eta_lasso_fullrank <- function(eta, sigma2, beta, eigenvalues, x, y, nt, myweights) {

  # this is based on the negative log-lik

  di <- 1 / myweights + eta * (eigenvalues - 1 / myweights)

  (nt / 2) * log(2 * pi) +
    (nt / 2) * log(sigma2) +
    0.5 * sum(log(di)) +
    (1 / (2 * sigma2)) * sum((y - x %*% beta)^2 / di)
}

fn_e_lasso_fullrank <- function(e, sigma2, beta, eigenvalues, x, y, nt, myweights) {
  eta = e2eta(e)
  fn_eta_lasso_fullrank(eta, sigma2, beta, eigenvalues, x, y, nt, myweights)
}

e2eta = function(e) {
  exp(-e)
}

eta2e = function(eta) {
  -log(eta)
}
