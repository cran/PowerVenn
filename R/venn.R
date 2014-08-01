#Compute overlap test and visualize intersections between multiple sets
#Author: Minghui Wang
#minghui.wang@mssm.edu
#
#setGeneric("PowerVenn", function(x, n=NULL,...)
#           standardGeneric("PowerVenn"),signature="x")
PowerVenn<-function(x,n=NULL,...){
#x, a list of sets
#n, integer, background population size
#..., additional arguments, not implemented

	if(!is.list(x)) stop('Input must be a list\n')
	if(is.null(names(x))) names(x)=paste('Set',1:length(x),sep='')
	set.names=names(x)
	if(any(set.names=='')) stop('Please specify names for each list entry\n')
	obj=list()
	obj$x=x
	obj$set.names=set.names
	obj$set.sizes=sapply(x,function(x) length(unique(x)))
	obj$n=n
	obj$overlap.sizes=inclusiveIntersect(listExclusiveIntersect(x))
	obj$P.value=NULL
	if(!is.null(n)){
		if(any(obj$set.sizes>n)) stop('Background population size should not be smaller than set size\n')
		obj$P.value=sapply(1:length(obj$overlap.sizes),function(i){
			which.set=which(strsplit(names(obj$overlap.sizes)[i],'')[[1]]=='1')
			if(length(which.set)==1) return(NA)
			cpsets(max(obj$overlap.sizes[i]-1,0),obj$set.sizes[which.set],n,lower.tail=FALSE)
		})
		names(obj$P.value)=names(obj$overlap.sizes)
	}
	class(obj)='venn'
	obj
}
#setMethod("PowerVenn", signature=c(x="list"), PowerVenn.list)
print.venn=function(x,...){
	cat('A venn object\n')
}
summary.venn=function(object,...){
	otab=object$overlap.sizes
	etab=rep(NA,length(otab))
	if(!is.null(object$n)){
		for(i in 1:length(otab)){
			s=strsplit(names(otab)[i],'')[[1]] == '1'
			if(sum(s)==1) next
			etab[i]=object$n*do.call('prod',as.list(object$set.sizes[s]/object$n))
		}
	}
	Barcode=names(otab)
	names(otab)=deBarcode(Barcode,object$set.names)
	res=list(Barcode=Barcode,otab=otab,etab=etab,set.names=object$set.names,set.sizes=object$set.sizes,n=object$n,P.value=object$P.value)
	#find intersections
	el=intersectElements(object$x)
	bc=strsplit(el$barcode,'')
	intersects=sapply(Barcode,function(d){
		id=which(strsplit(d,'')[[1]]=='1')
		od=sapply(bc,function(b) ifelse(all(b[id]=='1'),TRUE,FALSE))
		paste(el[od,1],collapse=', ')
	})
	if(is.null(object$n)){
		res$Table=data.frame(Intersections=names(otab),Observed.Overlap=otab,Intersects=intersects)
	}else{
		res$Table=data.frame(Intersections=names(otab),Observed.Overlap=otab,Expected.Overlap=etab,Fold.Change=otab/etab,P.value=object$P.value,Intersects=intersects,stringsAsFactors=FALSE)
		rownames(res$Table)=Barcode
	}
	class(res)='summary.venn'
	res
}
print.summary.venn=function(x,...){
	cat('A venn object with',length(x$set.names),'sets:',x$set.names,'\n')
	if(!is.null(x$n)) cat('Background size:',x$n,'\n')
	cat('Summary of intersections:\n')
	x$Table$Intersects=sapply(x$Table$Intersects,function(d){
		if(nchar(d)>20) d=paste(substr(d,0,20),' ...',sep='')
		d
	})
	print(x$Table)
}
#
plot.venn=function(x,Layout=c('circos','landscape'),which.intersections=0,keep.empty.intersections=TRUE,sort.by.size=FALSE,log.scale=FALSE,label.by.set.names=NULL,...){
#keep.empty.intersections, whether to retain empty intersections in the plot
	Layout <- match.arg(Layout)
	if(Layout=='circos'){
		return(plot.venn.circos(x=x,which.intersections=which.intersections,keep.empty.intersections=keep.empty.intersections,sort.by.size=sort.by.size,log.scale=log.scale,label.by.set.names=label.by.set.names,...))
	}else if(Layout=='landscape'){
		return(plot.venn.landscape(x=x,which.intersections=which.intersections,keep.empty.intersections=keep.empty.intersections,sort.by.size=sort.by.size,log.scale=log.scale,label.by.set.names=label.by.set.names,...))
	}else{
		stop('Invalid Layout\n')
	}
}
plot.venn.landscape=function(x,which.intersections=0,keep.empty.intersections=TRUE,sort.by.size=FALSE,log.scale=FALSE,label.by.set.names=NULL,...){
	if(is.null(label.by.set.names)) label.by.set.names=TRUE
	Args=list(...)
	cex=ifelse(is.null(Args$cex),0.8,Args$cex)
	#set color scheme
	if(is.null(Args$heatmapColor)){
		heatmapColor = c('#eeeeee',rev(heat.colors(100)[1:50]))
	}
	nColors=length(heatmapColor)-1
	params=getPlotParams(x,nColors,which.intersections=which.intersections,keep.empty.intersections=keep.empty.intersections,sort.by.size=sort.by.size,log.scale=log.scale)
	otab=params$otab
	otab0=params$otab0
	nO=length(otab)
	nSet=length(x$set.sizes) #number of sets
	cid=params$cid
	mlogp=params$mlogp

	#start plotting
	#start a canvas
	grid.newpage()
	vp <- viewport(x=0.53,y=0.5,width=0.9, height=0.9)
	pushViewport(vp)
	###grid.rect(x=0.5, y=0.5, width=1,height=1,gp=gpar(fill="gray80"))
	if(nSet==1){
		grid.circle(x=0.5, y=0.5, r=0.3)
		grid.text(otab0[1],0.5,0.5,gp=gpar(cex=cex*2.5))
		return(invisible())
	}
	if(nSet==2){
	#to be modified
		grid.circle(x=0.4, y=0.5, r=0.2)		
		grid.circle(x=0.6, y=0.5, r=0.2)
		grid.text(otab0['10'],0.3,0.5,gp=gpar(cex=cex*2.5))
		grid.text(otab0['11'],0.5,0.5,gp=gpar(cex=cex*2.5))
		grid.text(otab0['01'],0.7,0.5,gp=gpar(cex=cex*2.5))
		return(invisible())
	}
	#sub canvas 1
	vp1 <- viewport(x=0.5, y=0.60, width=0.95, height=0.8)
	pushViewport(vp1)
	
	yLen=axisLen(max(otab0,na.rm=T))
	w=1/nO
	h=0.9/yLen
	#plot intersections
	for(i in 1:nO){
		posx=w*(i-1)+w/2
	#	grid.rect(x=posx,y=0.05,width=w*0.8,height=h*otab[i],just=c('center','bottom'),gp=gpar(fill=heatmapColor[cid[i]]))
	#	grid.text(otab0[i],posx,0.03,rot=45,gp=gpar(cex=cex),just=c('right','top')) #size
		grid.rect(x=posx,y=0.0,width=w*0.8,height=h*otab[i],just=c('center','bottom'),gp=gpar(fill=heatmapColor[cid[i]]))
	}
	upViewport()

	#y axis
	vp1 <- viewport(x=0.05, y=0.60, width=0.1, height=0.8)
	pushViewport(vp1)
	grid.yaxis(at=seq(0,0.9,length.out=6),label=seq(0,yLen,length.out=6))
	upViewport()
	
	#color scale
	if((!is.null(x$n)) & (! is.null(mlogp))){
		vp11 <- viewport(x=0.9, y=0.95, width=0.2, height=0.1)
		pushViewport(vp11)
		grid.text('-Log10(P value)',0.5,0.75,just=c('center','bottom'),gp=gpar(cex=cex))
		wc=1/length(heatmapColor)
		for(i in 2:length(heatmapColor)){
			grid.polygon(c((i-1)*wc,i*wc,i*wc,(i-1)*wc),c(0.45,0.45,0.7,0.7),gp=gpar(fill=heatmapColor[i],col=NA))
		}
		t1=floor(min(mlogp,na.rm=T));t2=ceiling(max(mlogp,na.rm=T))
		grid.text(t1,0,0.3,just=c('center','top'),gp=gpar(cex=cex))
		grid.lines(x = c(wc, wc),y = c(0.45, 0.35))
		grid.text(t2,1,0.3,just=c('center','top'),gp=gpar(cex=cex))
		grid.lines(x = c(1-wc, 1-wc),y = c(0.45, 0.35))
		t3=(t1+t2)/2
		if(t3-t1>2){
			grid.text(as.integer(t3),0.5,0.3,just=c('center','top'),gp=gpar(cex=cex))
			grid.lines(x = c(0.5, 0.5),y = c(0.45, 0.35))
		}
		upViewport()
	}

	#sub canvas 2
	vp2 <- viewport(x=0.5, y=0.1, width=0.95, height=0.2)
	pushViewport(vp2)
	h=0.9/nSet
	for(i in 1:nO){
		posx=w*(i-1)+w/2
		for(j in 1:nSet){
			vpJ <- viewport(x=posx, y=0.01+(j-0.5)*h, width=w*0.8, height=h*0.75)
			pushViewport(vpJ)
			if(substr(names(otab[i]),j,j)=='1'){
				grid.circle(0.5,0.5,0.5,gp=gpar(fill='#0000aa'))
			}else{
				grid.circle(0.5,0.5,0.5)
			}
			upViewport()
		}
	}
	for(j in 1:nSet){
		grid.text(ifelse(label.by.set.names,x$set.names[j],j), 0, 0.015+(j-0.5)*h,just=c('right','center'))
		grid.text(x$set.sizes[j], 1, 0.015+(j-0.5)*h,just=c('left','center'))
	}
	upViewport()
	return(invisible())
}
plot.venn.circos=function(x,which.intersections=0,keep.empty.intersections=TRUE,sort.by.size=FALSE,log.scale=FALSE,label.by.set.names=NULL,...){
	if(is.null(label.by.set.names)) label.by.set.names=FALSE
	Args=list(...)
	cex=ifelse(is.null(Args$cex),0.8,Args$cex)
	#set color scheme
	if(is.null(Args$heatmapColor)){
		heatmapColor = c('#eeeeee',rev(heat.colors(100)[1:50]))
	}
	nColors=length(heatmapColor)-1
	params=getPlotParams(x,nColors,which.intersections=which.intersections,keep.empty.intersections=keep.empty.intersections,sort.by.size=sort.by.size,log.scale=log.scale)
	otab=params$otab
	otab0=params$otab0
	nO=length(otab)
	nSet=length(x$set.sizes) #number of sets
	cid=params$cid
	mlogp=params$mlogp

	# set graph layout parameters
	width.sets=0.3
	width.intersections=0.5-width.sets
	track.offset=2 #number of phantom tracks in the middle
	track.width=width.sets/(nSet+track.offset)
	bar.width.unit=width.intersections/max(otab,na.rm=T)
	gap.within.track=0.1 #ratio of gap width over block width on the same track
	gap.between.track=0.1 #ratio of gap width over track width

	#start a canvas
	grid.newpage()

	#Plot tracks
	origin=c(0.5,0.5)
	vp <- viewport(x=origin[1],y=origin[2],width=0.95, height=0.95)
	pushViewport(vp)
	degreeUnit=2*pi/nO
	degreeStart=(c(1:nO)-1)*degreeUnit
	degreeEnd=(c(1:nO))*degreeUnit
	degree.gap=max(degreeUnit*gap.within.track,2*pi/360)

	fill.col=rep(c('#dddddd','#999999'),length.out=nSet)
	#plot tracks
	#for(i in 1:nSet){
	#	XY1=sapply(seq(0,2*pi,length.out=360), function(deg) getXY(origin,(i+track.offset-1)*track.width,deg))
	#	XY2=sapply(seq(0,2*pi,length.out=360), function(deg) getXY(origin,(i+track.offset)*track.width-track.width*gap.between.track,deg))
	#	pos.x=c(XY1[1,],rev(XY2[1,]))
	#	pos.y=c(XY1[2,],rev(XY2[2,]))
	#	grid.polygon(pos.x, pos.y,gp=gpar(col=1))
	#}
	#plot overlap
	for(i in 1:nO){
		which.set=strsplit(names(otab)[i],'')[[1]]=='1'
		for(j in 1:nSet){
			XY1=sapply(seq(degreeStart[i],degreeEnd[i]-degree.gap,length.out=40), function(deg) getXY(origin,(j+track.offset-1)*track.width,deg))
			XY2=sapply(seq(degreeStart[i],degreeEnd[i]-degree.gap,length.out=40), function(deg) getXY(origin,(j+track.offset)*track.width-track.width*gap.between.track,deg))
			pos.x <- c(XY1[1,],rev(XY2[1,]))
			pos.y <- c(XY1[2,],rev(XY2[2,]))
			grid.polygon(pos.x, pos.y,gp=gpar(col = 'black',fill=ifelse(which.set[j],'yellow','#eeeeee'))) #
		}
		#bar of intersection
		XY1=sapply(seq(degreeStart[i],degreeEnd[i]-degree.gap,length.out=40), function(deg) getXY(origin,width.sets,deg))
		XY2=sapply(seq(degreeStart[i],degreeEnd[i]-degree.gap,length.out=40), function(deg) getXY(origin,width.sets+bar.width.unit*otab[i],deg))
		pos.x=c(XY1[1,],rev(XY2[1,]))
		pos.y=c(XY1[2,],rev(XY2[2,]))
		grid.polygon(pos.x, pos.y,gp=gpar(fill = heatmapColor[cid[i]],col=1)) #bar width is proportional to the intersection size
		#text intersection size
		XY3=sapply(seq(degreeStart[i],degreeEnd[i]-degree.gap,length.out=4), function(deg) getXY(origin,width.sets+0.01+bar.width.unit*otab[i],deg))
		grid.text(otab0[i],mean(XY3[1,]),mean(XY3[2,]),rot=degreeStart[i],gp=gpar(cex=cex))
	}
	#track number
	for(j in 1:nSet){
		XY1=sapply(seq(degreeStart[1],degreeEnd[1]-degree.gap,length.out=40), function(deg) getXY(origin,(j+track.offset-1)*track.width,deg))
		XY2=sapply(seq(degreeStart[1],degreeEnd[1]-degree.gap,length.out=40), function(deg) getXY(origin,(j+track.offset)*track.width,deg))
		grid.text(ifelse(label.by.set.names,x$set.names[j],j),(XY1[1,1]+XY2[1,1])/2,mean(XY1[2,]),just=c('center'),gp=gpar(cex=cex))
	}
	upViewport()
	#color scale
	if((!is.null(x$n)) & (! is.null(mlogp))){
		vp11 <- viewport(x=0.85, y=0.95, width=0.2, height=0.1)
		pushViewport(vp11)
		grid.text('-Log10(P value)',0.5,0.75,just=c('center','bottom'),gp=gpar(cex=cex))
		wc=1/length(heatmapColor)
		for(i in 2:length(heatmapColor)){
			grid.polygon(c((i-1)*wc,i*wc,i*wc,(i-1)*wc),c(0.45,0.45,0.7,0.7),gp=gpar(fill=heatmapColor[i],col=NA))
		}
		t1=floor(min(mlogp,na.rm=T));t2=ceiling(max(mlogp,na.rm=T))
		grid.text(t1,0,0.3,just=c('center','top'),gp=gpar(cex=cex))
		grid.lines(x = c(wc, wc),y = c(0.45, 0.35))
		grid.text(t2,1,0.3,just=c('center','top'),gp=gpar(cex=cex))
		grid.lines(x = c(1-wc/2, 1-wc/2),y = c(0.45, 0.35))
		t3=(t1+t2)/2
		if(t3-t1>2){
			grid.text(as.integer(t3),0.5,0.3,just=c('center','top'),gp=gpar(cex=cex))
			grid.lines(x = c(0.5, 0.5),y = c(0.45, 0.35))
		}
		upViewport()
	}
	return(invisible())
}
getXY=function(origin,radius,degree){
	X=radius*cos(degree)
	Y=radius*sin(degree)
	origin+c(X,Y)
}
getPlotParams=function(x,nColors=50,which.intersections=0,keep.empty.intersections=TRUE,sort.by.size=FALSE,log.scale=FALSE){
	otab=x$overlap.sizes
	if(sort.by.size==TRUE){
		otab.order=order(otab,decreasing=TRUE)
		otab=otab[otab.order]
	}
	if(keep.empty.intersections==FALSE){
		otab.kept=otab>0
		otab=otab[otab.kept]
	}else{
		otab.kept=rep(T,length(otab))
	}
	if(which.intersections>0){
		kpt=sapply(names(otab),function(d) sum(strsplit(d,'')[[1]] == '1') == which.intersections)
		if(sum(kpt)<2) stop('Too few intersections left for plotting after applying filter with which.intersections\n')
		otab.kept[otab.kept]=kpt
		otab=otab[kpt]
	}
	otab0=otab
	if(log.scale==TRUE) otab=log(otab+1)
	nO=length(otab)
	cid=rep(1,nO) #color gradient id
	mlogp=NULL
	if((!is.null(x$n)) & nO>0){
		mlogp=-log(x$P.value,base=10)
		mlogp[mlogp == Inf] = max(320,mlogp[mlogp < Inf],na.rm=T)
		if(sort.by.size==TRUE) mlogp=mlogp[otab.order]
		mlogp=mlogp[otab.kept]
		cid=ceiling(nColors*mlogp/max(c(mlogp,1e-10),na.rm=T))+1
		cid[is.na(cid)]=1
		if(all(is.na(mlogp))) mlogp=NULL
	}
	return(list(otab=otab,otab0=otab0,cid=cid,mlogp=mlogp))
}
axisLen=function(Max,Min=0,n.ticks=6){
	u=ceiling(Max/(n.ticks-1))
	if(u %% 5 != 0){
		u= ((u+5) %/% 5)*5
	}
	u*(n.ticks-1)
}
