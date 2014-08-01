\name{MSET}
\alias{MSET}
\docType{package}
\title{Multi-Set Overlap Test}
\description{
  Calculate the probability of overlap between multiple sets.
}
\usage{
MSET(x,n,lower.tail=TRUE,log.p=FALSE)
}
\arguments{
 \item{x}{list; a collection of sets.} 
 \item{n}{integer; background population size.}
 \item{lower.tail}{logical; if TRUE, probability is \code{P[overlap < m]}, otherwise, 
	\code{P[overlap >= m]}, where \code{m} is the number of elements overlap between all sets.}
 \item{log.p}{logical; if TRUE, probability p is given as log(p).}
}
\value{
Probability.
}
\author{
Minghui Wang <minghui.wang@mssm.edu>
}
\examples{
##not run###
#fake data
#x=list(S1=letters[1:20],S2=letters[10:26],S3=sample(letters,10),
# S4=sample(letters,10))
#MSET(x,26,FALSE)
}
