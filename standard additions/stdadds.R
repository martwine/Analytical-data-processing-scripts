# R program to calculate concentrations for raw data using calibration by standard additions. See README.txt for summary of .dat and .method files which are required inputs.
# Martin Johnson 11/2010 martin.johnson@uea.ac.uk
# Please use and distribute freely. But please acknowledge the originator (me).

#dataname gives name of subdirectory, .method and .data files
standard_additions<-function(dataname){
	
	data_filename<-paste(dataname,"/",dataname,".data",sep="")
	method_filename<-paste(dataname,"/",dataname,".method",sep="")

	rawdata<-read.table(data_filename,header=TRUE,fill=TRUE)

	#filter out lines flagged as bad (rowflag==0)
	rawdata<-subset(rawdata,rawdata$rowflag==1)
	
	#calculate a mean from all 'good' instrument readings (bad ones should be strings not numerics e.g. D13.23 not 13.23)
	rawdata$mean<-apply(rawdata,1,function(row){mean(as.numeric(row[5:14]),na.rm=T)})

	methoddata<-read.table(method_filename,header=TRUE)
	sample_vol<-methoddata$sample_volume #voulme in l
	stock_concn<-methoddata$stock_concn #concn in mol/l
	
	# if there's a tracer is assumed to be a member of one calibration (normally the largest addition to a particular set)
	# these variables used to get tracer which is member of calibration
	# the other (repeat) tracer measurements will be identified by type=-333
	tracer_name<-as.character(methoddata$tracer_name) #name of sample used for tracer
	tracer_type<-methoddata$tracer_type #vol of standard addition for tracer instance

	#which 'type' to be used as baseline?
	bl_method<-methoddata$baseline_method

	#start a blank list for the calibration data
	calib_out<-list()
	
	#for each calibration 'set' in the file...
	for(seti in unique(rawdata$set)){
		
		rows<-rawdata[which(rawdata$set==seti),]
		
		#calculate a baseline value if required 
		# (some intruments effectivel provide baseline-corrected raw data so this isn't needed
		# - if so, set baseline_method to -7777 in method file
		if(bl_method==-7777){baseline<-0}else{
			baselinerows<-rows[which(rows$type==bl_method),]
			baseline<-mean(baselinerows$mean)}
		#get the rows which are actua standard additions
		additions_rows<-rows[which(rows$type>=0),]
		
		# add a column giving the concentration due to the addition of standard
		#  addition vol is in ul so need to conver to l by *1e-6
		additions_rows$add_concn<-(additions_rows$type*1e-6*stock_concn)/(sample_vol+additions_rows$type*1e-6)

		#calculate a regression line and get the gradient and intercept
		regression<-lm(additions_rows$mean-baseline~additions_rows$add_concn)
		gradient<-regression[[1]][[2]]
		y_intercept<-regression[[1]][[1]]

		# calculate the concentration in the sample
		xE<-y_intercept/gradient

		# do some statistics to get the undertainty 
		#  - see http://www.uea.ac.uk/~e356/#[[Error%20analysis%20in%20standard%20additions]]
		n<-length(additions_rows$mean)
		Syoverx<-summary(regression)[[6]]
		SxE<-(Syoverx/gradient)*sqrt((1/n)+((mean(additions_rows$mean)^2)/(sum(sum((additions_rows$add_concn-mean(additions_rows$add_concn))^2))*gradient^2)))
		#95% confidence interval
		u_xE<-SxE*qt(0.975,n-2)

		#add the calibration stats for each set to the calib list
		templist<-list()
		templist["sample_id"]<-as.character(additions_rows$Desc[1])
		templist["extrapolated concentration"]<-xE
		templist$gradient<-gradient
		templist$confidence<-u_xE
		templist$baseline<-baseline
		key<-paste(dataname,"_set",seti,sep="")
		calib_out[[key]]<-templist	

		#plot calibration data and fit line
		filename<-paste(dataname,"/",key,".png",sep="")
		linedatax<-c(-xE,1.3*max(additions_rows$add_concn))
		linedatay<-(gradient*linedatax)+y_intercept		
		png(filename)
			plot(linedatax,linedatay,type="l",xlab="concentration from standard addition",ylab="instrument response",main=key,xaxs="i",yaxs="i")			
			points(additions_rows$add_concn,additions_rows$mean-baseline)
			abline(v=0)
		dev.off()	
	}

	#now for all 'sets'


	#add in a hacky extra column to correct tracer values with the correct baseline set
	# can't figure out a way to get the value of the previous row in apply in the concentration
	# calculation, otherwise this wouldn't be necessary
	bl_set_column_hack<-NULL
	for (index in 1:nrow(rawdata)){
		row<-rawdata[index,]
		ifelse(row[,"type"]==-333&&!index==1, newcolval<-rawdata[(index-1),"set"], newcolval<-0)
		bl_set_column_hack<-c(bl_set_column_hack,newcolval)
	}
	#ugly ugly ugly *must get better at R*
	rawdata<-cbind(rawdata,bl_set_column_hack)
	

	# get all the non-baseline data  (in a  roundabout way in case there are multiple baseline types)
	data_out<-subset(rawdata,rawdata$type==-999|rawdata$type==-333|rawdata$type>=0,select=c("Desc","set","type","mean","bl_set_column_hack"))
	
	#select only single element of calibration - the lowest in concentration (probably 0 added unless that row discarded with roflag) 
	# to represent sample that's been calibrated
	# don't discard calibration rows which are equal to tracers
	data_out<-data_out[-which((duplicated(data_out$Desc)&data_out$type>=0)&!(data_out$Desc==tracer_name&data_out$type==tracer_type)),]



	#calculate concentration - by value-baseline/gradient for unknowns and tracers, 
	# and by getting relevant extrpolated concentrations from calib_out for directly calibrated samples
	data_out$concentration<-apply(data_out,1,function(row){
		key<-paste(dataname,"_set",row[2],sep="")
		ifelse(as.numeric(row[3])==-333&!row[[5]]==0,bl_key<-paste(dataname,"_set",row[5],sep=""),bl_key<-key)	
		ifelse((as.numeric(row[3])>=0)&&!(as.character(row[1])==tracer_name&&as.numeric(row[3])==tracer_type),calib_out[[key]][[2]],(as.numeric(row[4])-calib_out[[bl_key]][[5]])/calib_out[[key]][[3]])
		})

	#apply uncertainty - at the moment this is just the extrapolation uncertainty for direclty calibrated - need to do better stats on unknowns
	data_out$u<-apply(data_out,1,function(row){
		key<-paste(dataname,"_set",row[2],sep="")
		calib_out[[key]][[4]]
	})

	#change type of tracer-equivalent calibration point to -333 to be clearly a tracer in the output
	data_out$type<-apply(data_out,1,function(row){
		ifelse(as.character(row[1])==tracer_name&&as.numeric(row[3])==tracer_type,row[3]<--333,row[3])
	})	

		
	data_out<-subset(data_out,select=c("Desc","set","type","mean","concentration","u"))
	#write out the calibration data and the results file
	dput(calib_out,file=paste(dataname,"/",dataname,".calib",sep=""))
	write.table(as.data.frame(data_out),file=paste(dataname,"/",dataname,".concn",sep=""),sep=",")

}

