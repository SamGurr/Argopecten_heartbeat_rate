#This algorithm was developed by Ian P. Dwyer and Samuel J. Gurr
#in January 2017 for the analysis of data produced by oxygen sensors
#around Long Island, New York. Criteria below define the duration
#and character of hypoxic events and offer a series of options for
#analysis and reporting of these events.


#User-defined variables:

Thresholds<-c(.99,1.99,2.99,3.49,4.79); #Choose any number of thresholds for hypoxia (DO<=Threshold)
ThresholdLabels<-c("Severe","Moderately Severe", "Moderate","Marginal", "Sub-Optimal"); #Choose labels for these thresholds
#Length of the labels list must be the same as the list of thresolds, and must be in the same order. Format labels as "Label"

MinEvent<-59.9; #defined in minutes, the minimum duration of hypoxia to be considered an event 
#(slightly shorter than interval because exactly 1 hr events were being discarded for some reason.)
Interval<-60; #defined in minutes, the maximum distance (<=) between two measured hypoxic points that are part of same event
timeDigits<-5; #specifies the number of decimal places to which the time data are accurate.
#For reference, 5 decimal places is a resolution of slightly smaller than seconds.

#Define some desired outputs here as booleans. Don't worry about it for now.
#Might exclude this entirely. Not sure.

#USER SHOULD NOT TOUCH ANYTHING BELOW THIS LINE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Imports the data.
#Note that the following functions expect the columns to be 
#TIME, TIME_NUM_FORMAT, and DO in that order.
filename<-as.character(file.choose());
data<-data.frame(read.csv(filename,header=TRUE));

#Sorts the data by time in ascending order before anything else happens.
data<-data[order(data[[2]]),];

#Separates out the elements into their own variables for easy use.
timeString<-as.character(data[[1]]);
timeNum<-data[[2]];
DO<-(data[[3]]);

#Creates a character vector to store what will eventually be the report.
dataReport<-character();

#Creates variables(s) that will eventually be output to (a) file(s) (hopefully if I get this right).

type<-character();
startTime<-character();
endTime<-character();
duration<-numeric();
eventDO<-numeric();
startNum<-numeric();
endNum<-numeric();

avgAvgDO<-numeric();
avgDur<-numeric();

#Defines the minimum event and maximum distance of measurements in terms of time code (new unit=days).
minEventNum<-round(MinEvent/24/60,digits=timeDigits);#Rounded to the same level of precision as the time data
intervalNum<-round(Interval/24/60,digits=timeDigits);#Rounded to the same level of precision as the time data

#Here's where the analysis begins in earnest:
#Defines and initializes variables to hold the starting and ending row numbers of each run detected
runStart<-0;#must start at 0 in order to work.
runEnd<-0;#must start at 0 in order to work.

#Defines a vector for keeping track of the number of runs of each type.
runCounts<-numeric();

#Loop through for each threshold.
for(j in 1:length(Thresholds)){
  #Defines and initializes a variable for counting runs.
  NumberOfRuns<-0;
  #Defines and vector for keeping track of DOs.
  AvgDOs<-numeric();
  #defines a vector to keep track of durations
  durs<-numeric();

#The loop that iterates through each row of the data and finds runs:
for (i in 1:length(timeNum)){
  #if we're not currently in a run...
  if(runEnd==0){
    #if we've found a starting point, set the start and end of the run to that time.
    if(DO[i]<=Thresholds[j]){
      runStart<-i;
      runEnd<-i;
    };
  }#LINE TERMINATOR WILL MESS THIS UP
  #if we're already in a run...
  else{
    #if we haven't exceeded the maximum allowable gap between hypoxic points...
    if(round(timeNum[i]-timeNum[runEnd],digits=timeDigits)<=intervalNum){
      #and we've just found another applicable point, reset the end point to this new point.
      if(DO[i]<=Thresholds[j]){
        runEnd<-i;
      };
      #If we hit the end of the file during a run that meets our criteria.
      if(i==length(timeNum)&&(round(timeNum[runEnd]-timeNum[runStart],digits=timeDigits)>=minEventNum)){
        #DO SHIT HERE FOR ANALYSIS
        NumberOfRuns<-NumberOfRuns+1;
        avgDO<-mean(DO[runStart:runEnd]);
        avgDOround<-round(avgDO,digits=2);
        dataReport[length(dataReport)+1]<-(paste("Found",ThresholdLabels[j],"Event #",NumberOfRuns," Avg DO:",
        avgDOround," Start/End:",timeString[runStart],"-",timeString[runEnd],"(Cut Off By End Of File)"));
        #for the output reports
        startTime[length(startTime)+1]<-timeString[runStart];#stores the start time
        endTime[length(endTime)+1]<-NA;#stores the end time
        startNum[length(startNum)+1]<-round(timeNum[runStart],digits=timeDigits);#stores the start num
        endNum[length(endNum)+1]<-NA;#stores the end num
        #Actually we don't want to update DO for incomplete events.
        #AvgDOs[length(AvgDOs)+1]<-avgDO; #stores the avg DO for each event
        duration[length(duration)+1]<-NA;
        eventDO[length(eventDO)+1]<-NA;
        type[length(type)+1]<-ThresholdLabels[j];
        #Actually we don't want to update durations for averaging because this doesn't count
        #durs[length(durs)+1]<-(timeNum[runEnd]-timeNum[runStart])*24*60; #stores event durations
        
        #set these back to 0 for the next file, or else shit will go badly if I try to automate.
        runStart<-0;
        runEnd<-0;
      };
    }#LINE TERMINANTOR WILL MESS THIS UP
    #if we've exceeded the maximum allowable gap and not found an applicable point...
    else{
      #if the run we've found is at least the minimum length we've defined... 
     if(round(timeNum[runEnd]-timeNum[runStart],digits=timeDigits)>=minEventNum){
       #DO SHIT HERE FOR ANALYSIS
       NumberOfRuns<-NumberOfRuns+1;
       avgDO<-mean(DO[runStart:runEnd]);
       avgDOround<-round(avgDO,digits=2);
       #Contingency for if run starts at beginning of file
       if(runStart==1){dataReport[length(dataReport)+1]<-(
         paste("Found",ThresholdLabels[j],"Event #",NumberOfRuns," Avg DO:",
         avgDOround," Start/End:",timeString[runStart],"(Cut off by beginning of file)","-",timeString[runEnd]));
       #for the output reports
       startTime[length(startTime)+1]<-NA;#stores the start time
       endTime[length(endTime)+1]<-timeString[runEnd];#stores the end time
       startNum[length(startNum)+1]<-NA;#stores the start num
       endNum[length(endNum)+1]<-round(timeNum[runEnd],digits=timeDigits);#stores the end num
       #we don't want to store this in this case.
       #AvgDOs[length(AvgDOs)+1]<-avgDO; #stores the avg DO for each event
       duration[length(duration)+1]<-NA;
       eventDO[length(eventDO)+1]<-NA;
       type[length(type)+1]<-ThresholdLabels[j];
       #we don't want to store this in this case.
       #durs[length(durs)+1]<-(timeNum[runEnd]-timeNum[runStart])*24*60; #stores event durations
       }
       #Every other case
       else{dataReport[length(dataReport)+1]<-(paste("Found",ThresholdLabels[j],"Event #",NumberOfRuns," Avg DO:",
       avgDOround," Start/End:",timeString[runStart],"-",timeString[runEnd]));
       #for the output reports
       startTime[length(startTime)+1]<-timeString[runStart];#stores the start time
       endTime[length(endTime)+1]<-timeString[runEnd];#stores the end time
       startNum[length(startNum)+1]<-round(timeNum[runStart],digits=timeDigits);#stores the start num
       endNum[length(endNum)+1]<-round(timeNum[runEnd],digits=timeDigits);#stores the end num
       AvgDOs[length(AvgDOs)+1]<-avgDO; #stores the avg DO for each event
       duration[length(duration)+1]<-(timeNum[runEnd]-timeNum[runStart])*24*60;
       eventDO[length(eventDO)+1]<-avgDO;
       type[length(type)+1]<-ThresholdLabels[j];
       durs[length(durs)+1]<-(timeNum[runEnd]-timeNum[runStart])*24*60; #stores event durations
       }
       
     };
      #regardless of whether we ran analysis or not, set these varaibles to 0 so the loop
      #knows we're not currently in a run anymore.
      runStart<-0;
      runEnd<-0;
    };
  };
};
runCounts[j]<-NumberOfRuns; #stores the total events for this threshold detected.
dataReport[length(dataReport)+1]<-"";#adds a blank line to the data file to separate blocks.

avgAvgDO[j]<-mean(AvgDOs); #updates the average average DO for the second report file.
avgDur[j]=mean(durs);#updates the average average duration for the second report file.
};



#adds the summary lines to the data report.
dataReport[length(dataReport)+1]<-"Events Detected:";
for(k in 1:length(Thresholds)){
  dataReport[length(dataReport)+1]<-paste(ThresholdLabels[k],"Events:",runCounts[k]);
};

#compiles the reports for saving
report1<-data.frame(
type=type,
durationMinutes=duration,
startTime=startTime,
endTime=endTime,
startNum=startNum,
endNum=endNum,
eventDO=eventDO);

report1$type<-as.character(report1$type);
report1$startTime<-as.character(report1$startTime);
report1$endTime<-as.character(report1$endTime);

report2<-data.frame(
type=ThresholdLabels,
threshold=Thresholds,
numberFound=runCounts,
avgDurMinutes=avgDur,
avgAvgDO=avgAvgDO);

report2$type<-as.character(report2$type);

#saves out the report files.
output1<-gsub(".csv","_Events.csv",filename);
write.csv(report1, file=output1,row.names=FALSE);

output2<-gsub(".csv","_Summary.csv",filename);
write.csv(report2, file=output2,row.names=FALSE);

#prints the report.
print(dataReport);