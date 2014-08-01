#Utility to do set analysis
#Author: Minghui Wang, minghui.wang@mssm.edu
#Date: 20 July, 2014

#list all possible intersections
intersectElements=function(x){
#return Venn diagram entry sizes
#x: a list of sets
	if(!is.list(x)) stop('Input x must be list\n')
	nL=length(x)
	if(nL<2) stop('Input x should have at least two entries\n')
	x=lapply(x,unique)
	allE=unique(unlist(x))
	dat=matrix(0,nrow=length(allE),ncol=nL)
	for(i in 1:nL){
		dat[,i]=0+(allE %in% x[[i]])
	}
	barcode=apply(dat,1,function(r) paste(r,collapse=''))
	data.frame(Entry=allE,barcode=barcode,stringsAsFactors=FALSE)
}
#enumerate all overlap and non-overlap set sizes
listExclusiveIntersect=function(x){
#return mutual exclusive intersection sizes
#x: a list of sets
	intersects=intersectElements(x)
	nL=length(x)
	barcodes=mkBarcode(nL)
	otab=rep(0,length(barcodes))
	names(otab)=barcodes
	tab=table(intersects$barcode)
	otab[names(tab)]=tab
	otab
}
inclusiveIntersect=function(x){
#x, an object generated from function exclusiveIntersect
#return inclusive subset sizes
	otab=x
	otab[]=0
	rel=barcodeRelation(names(x))
	for(i in 1:length(x)){
		otab[rel[i,]]=otab[rel[i,]]+x[i]
	}
	otab
}
barcodeRelation=function(barcode){
#barcode, character strings of 0/1
#return matrix, row row specifies whether the row entry is a subset of the column entries.
	n=length(barcode)
	Mat=matrix(FALSE,n,n)
	Mat=sapply(barcode,function(a)
		sapply(barcode,function(b) is.subSet(a,b))
	)
	t(Mat)
}
is.subSet=function(a,b){
#a,b barcode
#return TRUE if a is b's subset; FALSE otherwise. Eg, a='00011' is a subset of b='00001', while a='00111' is not a subset of b='10001'.
	a1=strsplit(a,'')[[1]] == '1'
	b1=strsplit(b,'')[[1]] == '1'
	if(length(a1) != length(b1)) stop('Input a and b have different numbers of characters\n')
	all(a1[which(b1)])
}
#use barcode (character strings of 0/1) to denote the overlap/non-overlap sets
mkBarcode=function(n){
#n, number of sets
	if(n<1) stop('n must be positive integer')
	barcode=paste(rep('0',n),collapse='')
	for(i in 1:n){
		barcode0=barcode
		substr(barcode0,i,i)='1'
		barcode=c(barcode,barcode0)
	}
	sort(barcode[-1])
}
#reverse barcode
deBarcode <- function(barcode,grp){
	sapply(barcode,function(b){
		s=grp[strsplit(b,'')[[1]] == '1']
		s=paste('Intersect(',paste(s,collapse=','),')',sep='')
		s
	})
}
