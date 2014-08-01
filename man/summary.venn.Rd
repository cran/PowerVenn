\name{summary.venn}
\alias{summary.venn}
\docType{package}
\title{Summarize Venn Ojbect}
\description{
  This function outputs summary statistics of a venn object.
}
\usage{
\method{summary}{venn}(object,...)
}
\arguments{
 \item{object}{a \code{venn} object.} 
 \item{...}{additional arguments (not implemented).}
}
\value{
A list:
 \item{Barcode}{a vector of 0/1 character strings, representing the set composition of each intersection.}
 \item{otab}{a vector of observed intersection size between any combination of sets.}
 \item{etab}{a vector of expected intersection size between any combination of sets if background population size is specified.}
 \item{set.names}{set names.}
 \item{set.sizes}{set sizes.}
 \item{n}{background population size.}
 \item{P.value}{upper tail p value for each intersection if background population size n is specified.}
 \item{Table}{a data.frame containing otab, etab, fold change, p value and the overlap elements.}
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
}
\seealso{\code{\link{venn}}}
