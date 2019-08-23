# Functions for BAT 2.1
# By Mans Thulin
# mans@statistikkonsult.com

############

# Format time data

# Switch to next day if machine runs for more than 24 hours
fixTime<-function(timeString)
{
    hour<-as.numeric(substring(timeString,1,2))-24
    if(hour<10){hourString<-paste(0,hour,sep="")} else { hourString<-as.character(hour) }
    return(hourString)
}

# Adjust measurements to control for increasing OD in blank wells
adjustBlanks<-function(OD,adjustAlt,blankID)
{
    totalMeanBlank<-0
    diffTimes<-OD$Time
    
    # Alternative 1: adjust using the mean of all blanks
    if(adjustAlt==1)
    {
        for(a in blankID)
        {
            totalMeanBlank<-totalMeanBlank+mean(eval(parse(text=paste("OD$",as.character(a),sep=""))))
        }
        totalMeanBlank<-totalMeanBlank/length(blankID)
    }
    
    # Alternative 2: adjust using only the blanks at the same timepoint
    if(adjustAlt==2)
    {	
        matrixOfBlanks<-matrix(NA,dim(OD)[1],length(blankID))
        for(a in 1:length(blankID))
        {
            matrixOfBlanks[,a]<-eval(parse(text=paste("OD$",as.character(blankID[a]),sep="")))
        }
        totalMeanBlank<-matrix(c(diffTimes,rep(rowMeans(matrixOfBlanks),dim(OD)[2]-1)),dim(OD))
    }

    OD<-OD-totalMeanBlank
    OD$Time<-diffTimes
    OD[OD<0]<-0.0001
    return(OD)
}

referencePlotter<-function(namesList,OD,maskingLower,maskingUpper,calculatedValues,ncols,maskType=1,maskingLowerX=1,maskingUpperX=2)
{
    x<-OD$Time
    par(mfrow=c(1,ncols),cex=1.1)
    
    # Step through all reference wells:
    for(i in 1:length(namesList)){
        
        # Get the name of the well:
        whichWell<-namesList[i]
        
        # Mask:
        if(maskType==1)
        {
            # Vertical masking interval:
            ODmasked<-eval(parse(text=paste("subset(OD,",maskingLower,"<OD$",as.character(whichWell),")",sep="")))
            ODmasked<-eval(parse(text=paste("subset(ODmasked,",maskingUpper,">ODmasked$",as.character(whichWell),")",sep="")))
        } else
        {
            # Horizontal masking interval:
            ODmasked<-eval(parse(text=paste("subset(OD,",maskingLowerX,"<=Time)",sep="")))
            ODmasked<-eval(parse(text=paste("subset(ODmasked,",maskingUpperX,">=Time)",sep="")))
        }
        
        ###################################################
        
        # Fit curve
        maskedy<-log(eval(parse(text=paste("ODmasked$",as.character(whichWell),sep=""))))
 #       if(nrow(ODmasked)>0) { m<-lm(maskedy~ODmasked$Time) }
        
        Rnow<-round(cor(ODmasked$Time,maskedy),4)
        pointsUsedNow<-dim(ODmasked)[1]
        
        # Save plot for diagnostics
        y<-log(eval(parse(text=paste("OD$",as.character(whichWell),sep=""))))
        
        plot(x,y,type="l",lwd=2,main=paste("R=",as.character(Rnow),"in",as.character(whichWell),sep=" "),xlab="Time",ylab="ln(OD)",sub=paste(pointsUsedNow,"points used for estimate",sep=" "))
        grid()
        points(x,y,type="l",lwd=3)
        if(nrow(ODmasked)>0) { points(ODmasked$Time,maskedy,col=2,type="l",lwd=2) }
        if(maskType==1) { axis(2,at=log(c(maskingLower,maskingUpper)),col=4,tck=1,labels=rep("",2),lwd=2,lty=3) } else { axis(1,at=c(maskingLowerX,maskingUpperX),col=4,tck=1,labels=rep("",2),lwd=2,lty=3) }

     }
}

###################################################

referencePlotter2<-function(refNamesList,namesList,OD,maskingLower,maskingUpper,numRep,calculatedValues)
{
    fullNamesList<-c(refNamesList,namesList)
    x<-OD$Time
    
    # Step through all reference wells:
    for(i in 1:length(fullNamesList)){
        
        if(i>nrow(calculatedValues)) { calculatedValues<-rbind(calculatedValues,rep(NA,nrow(calculatedValues))) }
        
        # Get the name of the well:
        whichWell<-fullNamesList[i]
        
        # Add well name to data frame:
        calculatedValues$Well[i]<-whichWell
        
        # Mask:
        ODmasked<-eval(parse(text=paste("subset(OD,",maskingLower,"<OD$",as.character(whichWell),")",sep="")))
        ODmasked<-eval(parse(text=paste("subset(ODmasked,",maskingUpper,">ODmasked$",as.character(whichWell),")",sep="")))
        
        ###################################################
        
        # Fit curve
        if(nrow(ODmasked)>0)
        { 
            maskedy<-log(eval(parse(text=paste("ODmasked$",as.character(whichWell),sep=""))))
            m<-lm(maskedy~ODmasked$Time)
            calculatedValues$fittedValue[i]<-m$coefficients[2]
            calculatedValues$R[i]<-cor(ODmasked$Time,maskedy)
            
            calculatedValues$pointsUsed[i]<-dim(ODmasked)[1]
            
            # Calculate doubling time
            calculatedValues$doubTime[i]<-log(2)/m$coefficients[2]
            
            # Add masking interval to data frame
            calculatedValues$mask1[i]<-maskingLower
            calculatedValues$mask2[i]<-maskingUpper
            calculatedValues$masktype[i]<-1
        } else { calculatedValues$R[i]<-10; calculatedValues$warningMessages[i]<- "No growth detected - the well appears to be empty or blank." }
        
    }
    
    # Relative growth rate
    refMean<-mean(calculatedValues$fittedValue[1:length(refNamesList)])
    calculatedValues$growthRate<-calculatedValues$fittedValue/refMean
    
    # Group averages
    calculatedValues$groupMeanDoubTime<-calculatedValues$groupMeanGrowthRate<-calculatedValues$groupGrowthRateSD<-" "
    calculatedValues$groupMeanDoubTime[1:length(refNamesList)]<-mean(calculatedValues$doubTime[1:length(refNamesList)])
    calculatedValues$groupMeanGrowthRate[1:length(refNamesList)]<-1
    calculatedValues$groupGrowthRateSD[1:length(refNamesList)]<-sd(calculatedValues$growthRate[1:length(refNamesList)])
    
    j<-length(refNamesList)
    for(i in 1:length(namesList))
    {
        if((i+j) %% numRep == 0)
        {
            calculatedValues$groupMeanDoubTime[i+j]<-mean(calculatedValues$doubTime[(i+j-numRep+1):(i+j)],na.rm=TRUE)
            calculatedValues$groupMeanGrowthRate[i+j]<-mean(calculatedValues$growthRate[(i+j-numRep+1):(i+j)],na.rm=TRUE)
            calculatedValues$groupGrowthRateSD[i+j]<-sd(calculatedValues$growthRate[(i+j-numRep+1):(i+j)],na.rm=TRUE)
        }
    }
    return(calculatedValues)
}

###################################################

referencePlotter3<-function(whichWell,refNamesList,namesList,OD,maskingLower,maskingUpper,numRep,calculatedValues,maskType=1,maskingLowerX=1,maskingUpperX=2)
{
    x<-OD$Time
    
        i<-which(calculatedValues$Well==whichWell)
    
        # Mask:
        if(maskType<1.5)
        {
            # Vertical masking interval:
            ODmasked<-eval(parse(text=paste("subset(OD,",maskingLower,"<OD$",as.character(whichWell),")",sep="")))
            ODmasked<-eval(parse(text=paste("subset(ODmasked,",maskingUpper,">ODmasked$",as.character(whichWell),")",sep="")))
        } else
        {
            # Horizontal masking interval:
            ODmasked<-eval(parse(text=paste("subset(OD,",maskingLowerX,"<=Time)",sep="")))
            ODmasked<-eval(parse(text=paste("subset(ODmasked,",maskingUpperX,">=Time)",sep="")))
        }        
        ###################################################
        
        # Fit curve
        maskedy<-log(eval(parse(text=paste("ODmasked$",as.character(whichWell),sep=""))))
        m<-lm(maskedy~ODmasked$Time)
        calculatedValues$fittedValue[i]<-m$coefficients[2]
        calculatedValues$R[i]<-cor(ODmasked$Time,maskedy)
        
        calculatedValues$pointsUsed[i]<-dim(ODmasked)[1]
        
        # Calculate doubling time
        calculatedValues$doubTime[i]<-log(2)/m$coefficients[2]
        
        # Add masking interval to data frame
        if(maskType==1)
        { calculatedValues$mask1[i]<-maskingLower; calculatedValues$mask2[i]<-maskingUpper } else
        { calculatedValues$mask1[i]<-maskingLowerX; calculatedValues$mask2[i]<-maskingUpperX; calculatedValues$masktype[i]<-2 }    
        
    # Relative growth rate
    refMean<-mean(calculatedValues$fittedValue[1:length(refNamesList)])
    calculatedValues$growthRate<-calculatedValues$fittedValue/refMean
    
    # Group averages
    calculatedValues$groupMeanDoubTime[1:length(refNamesList)]<-mean(calculatedValues$doubTime[1:length(refNamesList)])
    calculatedValues$groupMeanGrowthRate[1:length(refNamesList)]<-1
    calculatedValues$groupGrowthRateSD[1:length(refNamesList)]<-sd(calculatedValues$growthRate[1:length(refNamesList)])
    
    j<-length(refNamesList)
    for(k in 1:length(namesList))
    {
        if((k+j) %% numRep == 0)
        {
            calculatedValues$groupMeanDoubTime[k+j]<-mean(calculatedValues$doubTime[(k+j-numRep+1):(k+j)])
            calculatedValues$groupMeanGrowthRate[k+j]<-mean(calculatedValues$growthRate[(k+j-numRep+1):(k+j)])
            calculatedValues$groupGrowthRateSD[k+j]<-sd(calculatedValues$growthRate[(k+j-numRep+1):(k+j)])
        }
    }
    return(calculatedValues)
}

#######################

finalPlotter<-function(namesList,OD,calculatedValues)
{
    x<-OD$Time
    
    # Step through all wells:
    for(i in 1:length(namesList)){
        
        # Get the name of the well:
        whichWell<-namesList[i]
        
        maskingLower<-calculatedValues$mask1[i]
        maskingUpper<-calculatedValues$mask2[i]
        
        if(is.na(maskingLower)) { maskingLower<- 10}
        if(is.na(maskingUpper)) { maskingUpper<- 11}
        
        # Mask:
        if(calculatedValues$masktype[i]<1.5)
        {
            # Vertical masking interval:
            ODmasked<-eval(parse(text=paste("subset(OD,",maskingLower,"<OD$",as.character(whichWell),")",sep="")))
            ODmasked<-eval(parse(text=paste("subset(ODmasked,",maskingUpper,">ODmasked$",as.character(whichWell),")",sep="")))
        } else
        {
            # Horizontal masking interval:
            ODmasked<-eval(parse(text=paste("subset(OD,",maskingLower,"<=Time)",sep="")))
            ODmasked<-eval(parse(text=paste("subset(ODmasked,",maskingUpper,">=Time)",sep="")))
        }             
        ###################################################
        
        # Fit curve
        maskedy<-log(eval(parse(text=paste("ODmasked$",as.character(whichWell),sep=""))))
        # if(nrow(ODmasked)>0) { m<-lm(maskedy~ODmasked$Time) }
        
        Rnow<-round(cor(ODmasked$Time,maskedy),4)
        pointsUsedNow<-dim(ODmasked)[1]
        
        # Save plot for diagnostics
        y<-log(eval(parse(text=paste("OD$",as.character(whichWell),sep=""))))
        
        plot(x,y,type="l",lwd=2,main=paste("R=",as.character(Rnow),"in",as.character(whichWell),sep=" "),xlab="Time",ylab="ln(OD)",sub=paste(pointsUsedNow,"points used for estimate",sep=" "))
        grid()
        points(x,y,type="l",lwd=3)
        if(nrow(ODmasked)>0) { points(ODmasked$Time,maskedy,col=2,type="l",lwd=2) }
        #axis(2,at=log(c(maskingLower,maskingUpper)),col=4,tck=0.025,labels=rep("",2),lwd=2)
        
    }
}