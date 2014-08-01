\name{plot.venn}
\alias{plot.venn}
\docType{package}
\title{Draw Multi-Set Intersections}
\description{
  This function draws intersections between multiple sets.
}
\usage{
\method{plot}{venn}(x,Layout=c('circos','landscape'),which.intersections=0,
     keep.empty.intersections=TRUE,sort.by.size=FALSE,
     log.scale=FALSE,label.by.set.names=NULL,...)
}
\arguments{
 \item{x}{a \code{venn} object.} 
 \item{Layout}{layout for plotting.}
 \item{which.intersections}{integer; which intersections are chosen for plotting. E.g., when \code{which.intersections=2}, only those pairwise intersections will be plotted. By default, \code{which.intersections=0}, all intersections are chosen.}
 \item{keep.empty.intersections}{logical; if FALSE, empty intersection(s) will be discarded to save plotting space.}
 \item{sort.by.size}{logical; if TRUE, intersections will be sorted by size in decreasing order. By default, intersections are sorted by set.}
 \item{log.scale}{logical; whether to plot with log transformed intersection sizes.}
 \item{label.by.set.names}{logical; if TRUE, sets will be labelled by their names, otherwise labelled by numeric id. By default NULL, the function will use set names in \code{landscape} layout, and use numeric id in \code{circos} layout.}
 \item{...}{additional arguments for plot function (not fully implemented).}
}
\value{
No return.
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
#plot(obj)
}
\seealso{\code{\link{venn}}}
