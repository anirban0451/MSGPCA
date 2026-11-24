#' Select a normalized Pareto ideal-point solution
#' @param df Candidate data frame; objectives must be finite numeric columns.
#' @param lambda1_col,lambda2_col Names of smoothing columns.
#' @param loss_col,roughness_col Names of unweighted reconstruction and roughness columns.
#' @return List of selected parameters, row_index and annotated results.
#' @details Min-max ranges use all supplied eligible rows. Constant ranges map to
#'   zero. Strict dominance requires a strict improvement in one objective.
#'   Exact distance ties select the first row in input order. Selection is grid-dependent.
#' @examples
#' d <- data.frame(s=1:3,l=1:3,R=c(1,2,4),P=c(4,2,1))
#' find_optimal_pareto_points(d,'s','l','R','P')
#' @export
find_optimal_pareto_points <- function(df,lambda1_col,lambda2_col,loss_col,roughness_col) {
  if (!is.data.frame(df) || !nrow(df) || !all(c(lambda1_col,lambda2_col,loss_col,roughness_col)%in%names(df))) stop('Invalid candidate data')
  R <- df[[loss_col]]; P <- df[[roughness_col]]
  if (!is.numeric(R)||!is.numeric(P)||any(!is.finite(c(R,P)))) stop('Objectives must be finite numeric values')
  scale01 <- function(x) if (max(x)>min(x)) (x-min(x))/(max(x)-min(x)) else rep(0,length(x))
  df$normalized_reconstruction <- scale01(R); df$normalized_roughness <- scale01(P)
  df$distance <- sqrt(df$normalized_reconstruction^2+df$normalized_roughness^2)
  df$pareto <- vapply(seq_along(R),function(i) !any(R<=R[i]&P<=P[i]&(R<R[i]|P<P[i])),TRUE)
  idx <- which.min(df$distance)
  list(optimal_lambda1=df[[lambda1_col]][idx],optimal_lambda2=df[[lambda2_col]][idx],
       reconstruction_loss=R[idx],spatial_roughness=P[idx],row_index=idx,results=df)
}
