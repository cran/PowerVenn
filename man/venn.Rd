\name{venn}
\alias{venn}
\title{Class to Contain Multi-Set Intersections}
\description{
  This object contains data regarding the intersections between multiple sets.
}
\details{
This is an object created by the \code{PowerVenn} function.}
\value{
 \item{x}{a list of sets from input.}
 \item{set.names}{names of the sets. If the input sets do not have names, they will be automatically named as SetX where X is an integer from 1 to the total number of sets.}
 \item{set.sizes}{a vector of set sizes.}
 \item{n}{background population size.}
 \item{overlap.sizes}{a named vector of intersection sizes. Each intersection component is named as a character strings of 0 and 1, where a value of 1 in the \code{i}th position of the string indicates the intersection is involved with the \code{i}th set; 0 otherwise. E.g., string '000101' indicates that the intersection is an overlap between the 4th and 6th sets.}
 \item{P.value}{a vector of p values for the intersections.}
}
\keyword{classes}
\seealso{\code{\link{PowerVenn}}, \code{\link{summary.venn}}, \code{\link{plot.venn}}}
