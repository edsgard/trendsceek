% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.R
\name{add_markdist_step}
\alias{add_markdist_step}
\title{Add step-gradient mark distributions to a point pattern}
\usage{
add_markdist_step(pp, low_marks, high_marks, step_border = 0.5)
}
\arguments{
\item{pp}{A point-pattern.}

\item{low_marks}{A numeric or numeric vector specifying the lower value of the step mark distribution. If a vector is given 'n' mark distributions are added where 'n' is the length of the vector.}

\item{high_marks}{A numeric or numeric vector specifying the upper value of the step mark distribution. If a vector is given 'n' mark distributions are added where 'n' is the length of the vector.}

\item{step_border}{A numeric specifying the relative x-position within the point-pattern window where the values of the marks will change from low to high.}
}
\value{
A point-pattern with added mark distributions.
}
\description{
\code{add_markdist_step} adds step-gradient mark distributions to a point pattern.
}
\examples{
low_expr = c(10, 10)
high_expr = c(15, 20)
pp = sim_pois(100)
pp = add_markdist_step(pp, low_expr, high_expr)

}

