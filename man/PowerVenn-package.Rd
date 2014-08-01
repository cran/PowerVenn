\name{PowerVenn}
\alias{PowerVenn}
\docType{package}
\title{Calculate Intersections Between Multiple Sets}
\description{
  This function calculates intersection sizes between multiple sets and optionally performs statistical tests of the intersection.
}
\usage{
PowerVenn(x,n=NULL,...)
}
\arguments{
 \item{x}{list; a collection of sets.}
 \item{n}{integer, background population size. Required for computing the statistical significance of overlap.}
 \item{...}{additional arguments (not implemented).}
}
\details{
This function calculates intersection sizes between multiple sets and, if background population size \code{n} is specified, performs statistical tests of the intersections.
}
\value{
An object of class \code{venn}.
}
\author{
Minghui Wang <minghui.wang@mssm.edu>
}
\examples{
##not run###
#fake data
#x=list(S1=letters[1:20],S2=letters[10:26],S3=sample(letters,10),
# S4=sample(letters,10))
#obj=PowerVenn(x,n=26)
#summary(obj)
#plot(obj)
}
\seealso{\code{\link{venn}}, \code{\link{MSET}}, \code{\link{cpsets}}}
