% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_group_fused_multi.R
\name{predict.groupFusedMulti}
\alias{predict.groupFusedMulti}
\title{Prediction for Group Fused Multi Lasso}
\usage{
\method{predict}{groupFusedMulti}(
  object,
  newx = NULL,
  s = NULL,
  sf = NULL,
  type = c("link", "response", "coefficients", "nonzero", "class", "nvars"),
  ...
)
}
\arguments{
\item{object}{fitted groupFusedMulti object}

\item{newx}{new matrix for predictions}

\item{s}{lambda value for the predictions. defaults to all values computed in the groupFusedMulti object}

\item{sf}{fused lasso tuning parameter value for the predictions. defaults to all values computed in the groupFusedMulti object}

\item{type}{type of predictions to be made.}

\item{...}{parameters to be passed to groupFusedMulti}
}
\value{
predictions or coefficients
}
\description{
Prediction for Group Fused Multi Lasso
}
