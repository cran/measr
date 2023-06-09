% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{create_profiles}
\alias{create_profiles}
\title{Generate mastery profiles}
\usage{
create_profiles(attributes)
}
\arguments{
\item{attributes}{Positive integer. The number of attributes being measured.}
}
\value{
A \link[tibble:tibble-package]{tibble} with all possible attribute
mastery profiles. Each row is a profile, and each column indicates whether
the attribute in that column was mastered (1) or not mastered (0). Thus,
the tibble will have \code{2^attributes} rows, and \code{attributes} columns.
}
\description{
Given the number of attributes, generate all possible patterns of attribute
mastery.
}
\examples{
create_profiles(3L)
create_profiles(5)
}
