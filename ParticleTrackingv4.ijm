//Motility analysis tools
//Designed to automatically track swimming microorganisms and return useful statistics about their behaviour
//Expects an input image to be a time series of images of microorganisms swimming

//Particle detection settings
var subtractBackground=true;	//Whether or not to subtract the average image over all frames from each individual frame prior to analysis
var gaussianBlurRadius=4;	//Radius of gaussian blur in pixels to use to simplify the image prior to particle detection
var blackBackground=true;	//Whether or not the background (non-cell) part of the image is black
var noiseTolerance=100;		//Noise tolerance in pixel intensity units for particle detection

//Particle track analysis settings
var timePerFrame=2;		//Time per frame, in 100s of miliseconds
var imageScale=0.65;		//Image scale in microns per pixel [NB. old default 0.1538, factor of 4.226 too small]
var maximumSearchRange=15;	//Distance in pixels over which to search to connect points into tracks

var minimumTrackLength=50;	//Minimum length of track to analyse, in 100s of miliseconds
var smoothingFactor=0;		//Rolling window size (in 100ms frames) for averaging out jitter in particle motion
var evaluationInterval=20;	//Window size (in 100ms frames) over which to calculate motion statistics

var tumblePersistenceThreshold=0.75; //Persistence to classify as tumbling behaviour
var tumbleVelocityThreshold=15;    //Velocity um/s to classify as tumbling behaviour

//Data arrays
var particleLocationsNumber=0;
var particleLocationsFrameOffset=0;
var particleLocationsX=newArray(0);
var particleLocationsY=newArray(0);
var particleTracksNumber=0;
var particleTracksX=newArray(0);
var particleTracksY=newArray(0);

//Image properties
var width=0;
var height=0;
var depth=0;
var sourceID=0;

//Other general variables
var displayedTrack=0;
var cropRadius=92;

//gaussianBlurRadius=0;
//detectParticles(1, nSlices()+1);
//buildParticleTracks(0, false);
//subtractDrift(0);

macro "Particle Detection Settings [q]" {
	Dialog.create("Particle Detection Settings");
	Dialog.addCheckbox("Subtract background", subtractBackground);
	Dialog.addCheckbox("Black background", blackBackground);
	Dialog.addNumber("Gaussian blur radius (px):", gaussianBlurRadius);
	Dialog.addNumber("Noise tolerance:", noiseTolerance);
	Dialog.show();
	subtractBackground=Dialog.getCheckbox();
	blackBackground=Dialog.getCheckbox();
	gaussianBlurRadius=maxOf(0, Dialog.getNumber());
	noiseTolerance=maxOf(1, round(Dialog.getNumber()));
}

macro "Preview Particle Detection [w]" {
	d=getSliceNumber();
	detectParticles(d, d+1);
}

macro "---" {
}

macro "Particle Track Analysis Settings [e]" {
	Dialog.create("Particle Track Analysis Settings");
	Dialog.addNumber("Time per frame (100s of ms):", timePerFrame);
	Dialog.addNumber("Image scale (microns per pixel):", imageScale);
	Dialog.addNumber("Maximum search range:", maximumSearchRange);
	Dialog.addNumber("Minimum track length (100s of ms):", minimumTrackLength);
	Dialog.addNumber("Smoothing factor (100s of ms):", smoothingFactor);
	Dialog.addNumber("Evaluation interval (100s of ms):", evaluationInterval);
	Dialog.show();
	timePerFrame=maxOf(1, round(Dialog.getNumber()));
	imageScale=Dialog.getNumber();
	maximumSearchRange=maxOf(1, Dialog.getNumber());
	minimumTrackLength=maxOf(1, round(Dialog.getNumber()));
	smoothingFactor=maxOf(0, round(Dialog.getNumber()));
	evaluationInterval=maxOf(1, round(Dialog.getNumber()));
}

macro "Recalculate Particle Track Statistics [r]" {
	recalculateStats(0, false);
}

macro "Print Population-wide Statistics [y]" {
	printPopulationStats(0, false);
}

var driftSubtractionIndex=0;

macro "Subtract Drift From Particle Tracks" {
	Dialog.create("Subtract Drift From Particle Tracks");
	options=newArray("Translation", "Translation & Scale");
	Dialog.addChoice("Drift subtraction method:", options, options[driftSubtractionIndex]);
	Dialog.show();
	choice=Dialog.getChoice();
	for (i=0; i<lengthOf(options); i++) {
		if (choice==options[i]) {
			driftSubtractionIndex=i;
		}
	}
	subtractDrift(driftSubtractionIndex);
	recalculateStats(0, false);
}

var intervalStart=1;
var intervalEnd=100;
var intervalStep=1;

macro "Calculate Population Persistence Statistics [t]" {
	Dialog.create("Calculate Population Persistence Statistics");
	Dialog.addNumber("Evaluation interval start (100s of ms):", intervalStart);
	Dialog.addNumber("Evaluation interval end (100s of ms):", intervalEnd);
	Dialog.addNumber("Evaluation interval step (100s of ms):", intervalStep);
	Dialog.show();
	intervalStart=Dialog.getNumber();
	intervalEnd=Dialog.getNumber();
	intervalStep=Dialog.getNumber();
	grabPersistencePlot(intervalStart, intervalEnd, intervalStep);
}	

macro "---" {
}


macro "Particle Motility Analysis [a]" {
	detectParticles(1, nSlices()+1);
	buildParticleTracks(0, false);
}

macro "---" {
}

macro "Particle Track Plotter [z]" {
	id=plotParticleTracks(false);
}

macro "Particle Track Plotter (centred tracks)" {
	id=plotParticleTracks(true);
}

macro "Display Track [x]" {
	selectImage(sourceID);
	Dialog.create("Display Track");
	Dialog.addCheckbox("Overwrite existing analysis", subtractBackground);
	Dialog.show();
	displayedTrack=minOf(maxOf(round(Dialog.getNumber()), 0), particleTracksNumber);
	start=extractTrack(displayedTrack);
}

macro "Extract Cell [c]" {
	selectImage(sourceID);
	Dialog.create("Display Track");
	Dialog.addNumber("TrackID", displayedTrack);
	Dialog.addNumber("Crop radius (px)", cropRadius);
	Dialog.show();
	displayedTrack=minOf(maxOf(round(Dialog.getNumber()), 0), particleTracksNumber);
	cropRadius=minOf(maxOf(round(Dialog.getNumber()), 1), particleTracksNumber);
	start=extractTrack(displayedTrack);
	getSelectionBounds(tx, ty, tw, th);
	if (tx!=0 || ty!=0 || tw!=getWidth() || th!=getHeight) {
		getSelectionCoordinates(x, y);
		newImage("trackID_"+displayedTrack, bitDepth()+"-bit black", cropRadius*2, cropRadius*2, lengthOf(x));
		setBatchMode(true);
		tgt=getImageID();
		for (i=0; i<lengthOf(x); i++) {
			selectImage(sourceID);
			setSlice(start+i+1);
			makeRectangle(x[i]-cropRadius, y[i]-cropRadius, cropRadius*2, cropRadius*2);
			setPasteMode("Copy");
			run("Copy");
			selectImage(tgt);
			setSlice(i+1);
			run("Paste");
		}
	}
}

macro "---" {
}

var batchOverwrite=false;
var batchRecursive=true;

macro "Batch Analyse Folder" {
	Dialog.create("Batch Analyse Folder");
	Dialog.addCheckbox("Overwrite existing analysis", batchOverwrite);
	Dialog.addCheckbox("Analyse subdirectories", batchRecursive);
	Dialog.addMessage("Click OK then select the directory to analyse in the next window");
	Dialog.show();
	batchOverwrite=Dialog.getCheckbox();
	batchRecursive=Dialog.getCheckbox();
	path=getDirectory("");
	analysePath(path, batchOverwrite, batchRecursive);
}

function detectParticles(startFrame, endFrame) {
	//Record source image id and dimensions
	sourceID=getImageID();
	width=getWidth();
	height=getHeight();
	depth=nSlices();
	//Exit if a stack is not found as the active image
	if (depth==1) {
		exit("Error: Stack required");
	}
	//Start analysis by setting batch mode and clearing any selection
	setBatchMode(true);
	run("Select None");
	//If requested, determine the background image
	if (subtractBackground==true) {
		if (blackBackground==true) {
			run("Z Project...", "start="+1+" stop="+depth+" projection=[Min Intensity]");
		} else {
			run("Z Project...", "start="+1+" stop="+depth+" projection=[Average Intensity]");
		}
		bg=getImageID();
	}
	//Setup the particle location array data
	particleLocationsFrameOffset=startFrame-1;
	particleLocationsX=newArray(0);
	particleLocationsY=newArray(0);
	particleLocationsNumber=0;
	//For each frame to analyse
	for (i=0; i<endFrame-startFrame; i++) {
		showProgress((i+1)/(endFrame-startFrame));
		showStatus("Finding particles...");
		//Select the source image
		selectImage(sourceID);
		setSlice(i+startFrame);
		//Duplicate the current slice for processing and particle analysis
		run("Duplicate...", "title=tmp");
		tmp=getImageID();
		//Set to 32bit if non-black background to allow negative signal after background subtraction
		if (blackBackground!=true) {
			run("32-bit");
		}
		//If requested, subtract the bacground
		if (subtractBackground==true) {
			selectImage(bg);
			run("Select All");
			setPasteMode("Subtract");
			run("Copy");
			selectImage(tmp);
			run("Select All");
			run("Paste");
			run("Select None");
			setPasteMode("Copy");
		}
		//If set, perform a mean-fold of the signal for a non-black image background
		if (blackBackground!=true) {
			selectImage(tmp);
			run("Select All");
			getRawStatistics(area, mean);
			if (mean<0) {
				run("Macro...", "code=v=abs(v+"+abs(mean)+")");
			} else {
				run("Macro...", "code=v=abs(v-"+mean+")");
			}
		}
		//Upscale for sub-pixel particle detection
		us=1; //Upscale factor
		if (us!=1) {
			run("Size...", "width="+getWidth()*us+" height="+getHeight()*us+" constrain average interpolation=Bilinear");
		}
		//If non-zero, perform a gaussian blur of the image
		if (gaussianBlurRadius!=0) {
			selectImage(tmp);
			run("Gaussian Blur...", "sigma="+gaussianBlurRadius*us);
		}
		//Find particle locations with the user specified noise tolerance, and correct for the upscale
		run("Find Maxima...", "noise="+noiseTolerance+" output=[Point Selection]");
		getSelectionCoordinates(x, y);
		if (us!=1) {
			for (j=0; j<lengthOf(x); j++) {
				x[j]=x[j]/us;
				y[j]=y[j]/us;
			}
		}
		//If necessary, pad the length of the selection coordinates to the number of particles recorded in the data arrays
		//Otherwise, if necessary, pad the length of the data arrays to accomadate the selection coordinates
		if (particleLocationsNumber>lengthOf(x)) {
			//Pad with invalid locations (ie. negative positions)
			paddingData=newArray(particleLocationsNumber-lengthOf(x));
			Array.fill(paddingData, -1);
			x=Array.concat(x, paddingData);
			y=Array.concat(y, paddingData);
		} else if (particleLocationsNumber<lengthOf(x)) {
			//Generate new particle locations data arrays to hold data
			particleLocationsXNew=newArray(0);
			particleLocationsYNew=newArray(0);
			//Pad with invalid locations (ie. negative positions)
			paddingData=newArray(lengthOf(x)-particleLocationsNumber);
			Array.fill(paddingData, -1);
			for (l=0; l<=i; l++) {
				particleLocationsXTemp=Array.slice(particleLocationsX, particleLocationsNumber*(l), particleLocationsNumber*(l+1));
				particleLocationsYTemp=Array.slice(particleLocationsY, particleLocationsNumber*(l), particleLocationsNumber*(l+1));
				particleLocationsXTemp=Array.concat(particleLocationsXTemp, paddingData);
				particleLocationsYTemp=Array.concat(particleLocationsYTemp, paddingData);
				particleLocationsXNew=Array.concat(particleLocationsXNew, particleLocationsXTemp);
				particleLocationsYNew=Array.concat(particleLocationsYNew, particleLocationsYTemp);
			}
			particleLocationsX=particleLocationsXNew;
			particleLocationsY=particleLocationsYNew;
			particleLocationsNumber=lengthOf(x);
		}
		//Append the newly detected particles for this time point to the array
		particleLocationsX=Array.concat(particleLocationsX, x);
		particleLocationsY=Array.concat(particleLocationsY, y);
		//Close the temporary image and display the particles detected on the source image
		selectImage(tmp);
		close();
		selectImage(sourceID);
		setSlice(i+startFrame);
		makeSelection("Points", x, y);
	}
	//Tidy up open images
	if (subtractBackground==true) {
		selectImage(bg);
		close();
	}
	selectImage(sourceID);
}

function buildParticleTracks(file, useFile) {
	//Reset the data arrays and variables for recording the track data
	particleTracksX=newArray(0);
	particleTracksY=newArray(0);
	particleTracksNumber=0;
	//Loop through all slices
	for (i=0; i<depth-1; i++) {
		showProgress((i+1)/depth);
		//Loop through all particles on each slice
		for (j=0; j<particleLocationsNumber; j++) {
			showStatus("Building particle tracks... Slice: "+i+" Track: "+particleTracksNumber);
			//If the particle has a valid location (i.e. is not negative)...
			if (particleLocationsX[j+i*particleLocationsNumber]!=-1) {
				//Start a track at that particle
				//Setup arrays to record the track
				trackX=newArray(depth);
				trackY=newArray(depth);
				//Populate arrays with invalid locations (ie. negative) prior to analysis
				Array.fill(trackX, -1);
				Array.fill(trackY, -1);
				//Record the start location
				trackX[i]=particleLocationsX[j+i*particleLocationsNumber];
				trackY[i]=particleLocationsY[j+i*particleLocationsNumber];
				//Set a flag to indicate that the track should be continued
				cont=true;
				//Loop through all slices from that point onwards
				for (l=i; l<depth-1; l++) {
					//But only continue if the flag allows it
					if (cont==true) {
						//Setup a distance measure and best match index for finding the nearest particle on the next slice
						d1=pow(maximumSearchRange, 2);
						index=-1;
						//Loop through all particles on the next slice
						for (k=0; k<particleLocationsNumber; k++) {
							//Grab the current location and the test particle location
							x1=trackX[l];
							y1=trackY[l];
							x2=particleLocationsX[k+(l+1)*particleLocationsNumber];
							y2=particleLocationsY[k+(l+1)*particleLocationsNumber];
							//If the test particle has a valid location
							if (x2!=-1) {
								//Determine the distance from the current location to the test particle
								if (l-i==0) {
									//Simple distance squared if the first time step
									d2=pow(x1-x2, 2)+pow(y1-y2, 2);
								} else {
									//Otherwise distance squared from the predicted particle location based on the previous two time points
									x3=trackX[l-1];
									y3=trackY[l-1];
									x4=x1+x1-x3;
									y4=y1+y1-y3;
									d2=pow(x4-x2, 2)+pow(y4-y2, 2);
								}
								//If the distance is better than the min distance and any previous match then record the distance and match index
								if (d2<d1) {
									d1=d2;
									index=k;
								}
							}
						}
						//If a particle within the minimum connection distance has been found then record the track data
						//Otherwise flag the track to finish
						if (index!=-1) {
							trackX[l+1]=particleLocationsX[index+(l+1)*particleLocationsNumber];
							trackY[l+1]=particleLocationsY[index+(l+1)*particleLocationsNumber];
							//Set the particle data for that point to -1 to prevent any further tracks connecting to it too
							particleLocationsX[index+(l+1)*particleLocationsNumber]=-1;
							particleLocationsY[index+(l+1)*particleLocationsNumber]=-1;
						} else {
							cont=false;
							//Set l manually to break the loop
							l=depth;
						}
					}
				}
				//Append the new track data to the array and display the track, if the length criterion is met
				//Measure track length
				length=0;
				for (l=0; l<depth; l++) {
					if (trackX[l]!=-1) {
						length++;
					}
				}
				//Append data
				if (length*timePerFrame>=minimumTrackLength) {
					particleTracksX=Array.concat(particleTracksX, trackX);
					particleTracksY=Array.concat(particleTracksY, trackY);
					emptyArray=newArray(lengthOf(trackX));
					start=extractTrack(particleTracksNumber);
					//if (particleTracksNumber==0) {
					//	grabTrackStatistics(file, useFile, particleTracksNumber, true);
					//} else {
					//	grabTrackStatistics(file, useFile, particleTracksNumber, false);
					//}
					particleTracksNumber++;
				}
			}
		}
	}
	
	recalculateStats(file, useFile);
	
	showStatus(particleTracksNumber+ " particle tracks measured");
}

var directionBins=72;
var directionMax=PI
var directionNWeight=newArray(directionBins);
var directionVWeight=newArray(directionBins);

function recalculateStats(file, useFile) {
	selectImage(sourceID);
	//Arrays for non-track specific statistics
	directionNWeight=newArray(directionBins);
	directionVWeight=newArray(directionBins);
	//Do the analysis
	for (i=0; i<particleTracksNumber; i++) {
		showStatus("Calculating particle track statistics");
		showProgress((i+1)/particleTracksNumber);
		start=extractTrack(i);
		if (i==0) {
			grabTrackStatistics(file, useFile, i, true);
		} else {
			grabTrackStatistics(file, useFile, i, false);
		}
	}
}

function subtractDrift(mode) {
	selectImage(sourceID);
	//Concatenate all particle x/y locations and dx/dy steps
	px=newArray();
	py=newArray();
	dx=newArray();
	dy=newArray();
	for (i=0; i<particleTracksNumber; i++) {
		start=extractTrack(i);
		if (selectionType()!=-1) {
			getSelectionCoordinates(x, y);
			if (lengthOf(x)>0) {
				cdx=newArray(lengthOf(x)-1);
				cdy=newArray(lengthOf(y)-1);
				for (j=0; j<lengthOf(x)-1; j++) {
					cdx[j]=x[j]-x[j+1];
					cdy[j]=y[j]-y[j+1];
				}
				x=Array.slice(x, 0, lengthOf(cdx));
				y=Array.slice(y, 0, lengthOf(cdy));
				px=Array.concat(px, x);
				py=Array.concat(py, y);
				dx=Array.concat(dx, cdx);
				dy=Array.concat(dy, cdy);
			}
		}
	}
	if (lengthOf(px)>10) {
		//Fit to a straight line to determine mean translation and scaling drifts
		Fit.doFit("Straight Line", px, dx);
		xDrift=Fit.p(0);
		xScale=Fit.p(1);
		Fit.doFit("Straight Line", py, dy);
		yDrift=Fit.p(0);
		yScale=Fit.p(1);
		//Correct particle locations
		for (i=0; i<particleTracksNumber; i++) {
			start=extractTrack(i);
			if (selectionType()!=-1) {
				getSelectionCoordinates(x, y);
				if (lengthOf(x)>0) {
					for (j=0; j<lengthOf(x); j++) {
						if (mode==0) {
							particleTracksX[j+start+i*depth]=x[j]+xDrift*j;
							particleTracksY[j+start+i*depth]=y[j]+yDrift*j;
						} else if (mode==1) {
							particleTracksX[j+start+i*depth]=x[j]+xDrift*j+xScale*x[j]*j;
							particleTracksY[j+start+i*depth]=y[j]+yDrift*j+yScale*y[j]*j;
						}
					}
				}
			}
		}
	}
}

function grabTrackStatistics(file, useFile, trackID, header) {
	//Setup the start of the file/output
	if (header==true) {
		if (useFile!=true) {
			//Reset the output log
			print("\\Clear");
			print(
				"trackID", "trackLength(s)", "totalDisplacement(um)", "meanVelocity(um/s)", "totalDistance(um)",
				"meanSpeed(um/s)", "minSpeed(um/s)", "maxSpeed(um/s)", "stdevSpeed(um/s)",
				"meanPersistence", "minPersistence", "maxPersistence", "stdevPersistence",
				"tumblingPropensityV", "tumblingPropensityP",
				"positionX(px)", "positionY(px)", "displacementX(um)", "displacementY(um)"
			);
		} else {
			print(file,
				"trackID trackLength(s) totalDisplacement(um) meanVelocity(um/s) totalDistance(um) meanSpeed(um/s) minSpeed(um/s) maxSpeed(um/s) stdevSpeed(um/s) meanPersistence minPersistence maxPersistence stdevPersistence tumblingPropensityV tumblingPropensityP positionX(px) positionY(px) displacementX(um) displacementY(um)"
			);
		}
	}

	//Make a selection of the current track
	start=extractTrack(trackID);
	getSelectionBounds(tx, ty, tw, th);
	if (tx!=0 || ty!=0 || tw!=getWidth() || th!=getHeight) {
		getSelectionCoordinates(x, y);

		//Interpolate the track to 1 point per 100ms
		xScaled=newArray((lengthOf(x)-1)*timePerFrame+1);
		yScaled=newArray((lengthOf(y)-1)*timePerFrame+1);
		for (i=0; i<lengthOf(x)-1; i++) {
			for (j=0; j<timePerFrame; j++) {
				xScaled[j+i*timePerFrame]=x[i]+(x[i+1]-x[i])*j/timePerFrame;
				yScaled[j+i*timePerFrame]=y[i]+(y[i+1]-y[i])*j/timePerFrame;
			}
		}
		xScaled[lengthOf(xScaled)-1]=x[lengthOf(x)-1];
		yScaled[lengthOf(yScaled)-1]=y[lengthOf(y)-1];
		x=xScaled;
		y=yScaled;
		makeSelection("Polyline", x, y);

		//Smooth the track according to the smoothing factor
		if (smoothingFactor!=0) {
			xSmooth=newArray(lengthOf(xScaled)-smoothingFactor);
			ySmooth=newArray(lengthOf(yScaled)-smoothingFactor);
			for (i=0; i<lengthOf(xScaled)-smoothingFactor; i++) {
				sumX=0;
				sumY=0;
				for (j=0; j<smoothingFactor; j++) {
					sumX+=xScaled[i+j];
					sumY+=yScaled[i+j];
				}
				xSmooth[i]=sumX/smoothingFactor;
				ySmooth[i]=sumY/smoothingFactor;
			}
			x=xSmooth;
			y=ySmooth;
			makeSelection("Polyline", x, y);
		}

		if (lengthOf(x)>evaluationInterval*2) {
			//Return track statistics according to the evaluation interval
			//Total displacement and average velocity
			totaldX=(x[0]-x[lengthOf(x)-1]);
			totaldY=(y[0]-y[lengthOf(y)-1]);
			Dsum=pow(totaldX*totaldX+totaldY*totaldY, 0.5);				//D		//Displacement - total
			Vmean=Dsum/lengthOf(x);							//dD V		//Velocity - average

			Array.getStatistics(x, Xmin, Xmax, Xmean, Xstdev);			//X			//X position - minimum, maximum, mean, standard deviation
			Array.getStatistics(y, Ymin, Ymax, Ymean, Ystdev);			//Y			//Y position - minimum, maximum, mean, standard deviation
			
			//Variables for tumbling
			tumP=0;
			tumV=0;
			
			//Total distance travelled, speed and velocity
			dX=newArray(lengthOf(x)-evaluationInterval);
			dY=newArray(lengthOf(y)-evaluationInterval);
			dD=newArray(lengthOf(x)-evaluationInterval);
			S=newArray(lengthOf(x)-evaluationInterval);
			V=newArray(lengthOf(x)-evaluationInterval);
			for (i=0; i<lengthOf(x)-evaluationInterval; i++) {
				dX[i]=(x[i+evaluationInterval]-x[i]);
				dY[i]=(y[i+evaluationInterval]-y[i]);
				for (j=0; j<evaluationInterval; j++) {
					cdX=(x[i+1]-x[i]);
					cdY=(y[i+1]-y[i]);
					dD[i]=dD[i]+pow(cdX*cdX+cdY*cdY, 0.5);
				}
				S[i]=dD[i];
				V[i]=pow(dX[i]*dX[i]+dY[i]*dY[i], 0.5);
				//Record tumbling counts
				if (V[i]<tumbleVelocityThreshold) {
					tumV++;
				}
				//Record global direction stats
				angle=atan2(dY[i], dX[i]);
				while (angle<0) {
					angle+=PI*2;
				}
				while (angle>=PI*2) {
					angle-=PI*2;
				}
				bin=floor(directionBins*angle/(PI*2));
				if (dX[i]!=0 && dY[i]!=0) {
					directionNWeight[bin]=directionNWeight[bin]+1;
				}
				directionVWeight[bin]=directionVWeight[bin]+V[i];
			}
			Array.getStatistics(S, Smin, Smax, Smean, Sstdev);			//dD S		//Speed - minimum, maximum, mean, standard deviation
			Array.getStatistics(V, Vmin, Vmax, Vmean2, Vstdev);			//V			//Velocity - minimum, maximum, mean (alternative calculation method), standard deviation
			Lsum=Smean*(lengthOf(x)-evaluationInterval)/evaluationInterval;		//D		//Distance - total

			//Directional persistence & persistence tumble detection
			P=newArray(lengthOf(dX)-evaluationInterval);
			Pcount=0;
			for (i=0; i<lengthOf(dX)-evaluationInterval; i++) {
				Pnum=dX[i]*dX[i+evaluationInterval]+dY[i]*dY[i+evaluationInterval];
				Pden=V[i]*V[i+evaluationInterval];
				if (Pnum!=0) {
					P[Pcount]=(dX[i]*dX[i+evaluationInterval]+dY[i]*dY[i+evaluationInterval])/(V[i]*V[i+evaluationInterval]);
					//Record persistence thresholded tumbling
					if (P[Pcount]<tumblePersistenceThreshold) {
						tumP++;
					}
					Pcount++;
				}
			}
			P=Array.trim(P, Pcount);
			Array.getStatistics(P, Pmin, Pmax, Pmean, Pstdev);			//P		//Directional persistence - minimum, maximum, mean, standard deviation

			//Change in speed and velocity
			dV=newArray(lengthOf(dX)-evaluationInterval);
			dS=newArray(lengthOf(dX)-evaluationInterval);
			for (i=0; i<lengthOf(dX)-evaluationInterval; i++) {
				dV[i]=V[i]-V[i+evaluationInterval];
				dS[i]=S[i]-S[i+evaluationInterval];
			}
			Array.getStatistics(dV, dVmin, dVmax, dVmean, dVstdev);			//dV		//Change in velocity - minimum, maximum, mean, standard deviation
			Array.getStatistics(dS, dSmin, dSmax, dSmean, dSstdev);			//dS		//Change in speed - minimum, maximum, mean, standard deviation

			timeScale=1/(0.1*evaluationInterval);

			//    track    length (s)
			if (useFile!=true) {
				print(
					trackID, lengthOf(x)*0.1,
					Dsum*imageScale, Vmean*imageScale/0.1, Lsum*imageScale,
					Smean*imageScale*timeScale, Smin*imageScale*timeScale, Smax*imageScale*timeScale, Sstdev*imageScale*timeScale,
					Pmean, Pmin, Pmax, Pstdev,
					tumV/(lengthOf(x)-evaluationInterval), tumP/Pcount,
					Xmean, Ymean, totaldX*imageScale, totaldY*imageScale
				);
			} else {
				print(file,
					""+trackID+" "+lengthOf(x)*0.1+" "+
					Dsum*imageScale+" "+Vmean*imageScale/0.1+" "+Lsum*imageScale+" "+
					Smean*imageScale*timeScale+" "+Smin*imageScale*timeScale+" "+Smax*imageScale*timeScale+" "+Sstdev*imageScale*timeScale+" "+
					Pmean+" "+Pmin+" "+Pmax+" "+Pstdev+" "+
					tumV/(lengthOf(x)-evaluationInterval)+" "+tumP/Pcount+" "+
					Xmean+" "+Ymean+" "+totaldX*imageScale+" "+totaldY*imageScale
				);
			}
		}
	}
}

function extractTrack(trackID) {
	//Setup arrays to record the track data
	count=0;
	start=0;
	xs=newArray(depth);
	ys=newArray(depth);
	//Loop through all slices
	for (l=0; l<depth; l++) {
		//Only record the track data if it has a valid location
		if (particleTracksX[l+trackID*depth]!=-1) {
			if (count==0) {
				start=l;
			}
			xs[count]=particleTracksX[l+trackID*depth];
			ys[count]=particleTracksY[l+trackID*depth];
			count++;
		}
	}
	//Trim the track data to positions with valid locations
	xs=Array.trim(xs, count-1);
	ys=Array.trim(ys, count-1);
	//If the track is non-zero length then make a selection displaying it
	if (lengthOf(xs)!=0) {
		makeSelection("Polyline", xs, ys);
	}
	return start;
}

function grabPersistencePlot(intervalStart, intervalEnd, intervalStep) {
	print("\\Clear");
	print("Time", "DirAutocorrel", "DirAutocorrelStdev", "DirAutocorrelStderr", "SpAutocorrel", "SpAutocorrelStdev", "SpAutocorrelStderr", "N");
	timeScale=1/(0.1*evaluationInterval);
	for (autocorelInterval=intervalStart; autocorelInterval<=intervalEnd; autocorelInterval+=intervalStep) {
		showStatus("Calcuating population directional persistence. Interval: "+evaluationInterval/10);
		allPC=newArray();
		allSC=newArray();
		for (trackID=0; trackID<particleTracksNumber; trackID++) {
			start=extractTrack(trackID);
			getSelectionBounds(tx, ty, tw, th);
			if (tx!=0 || ty!=0 || tw!=getWidth() || th!=getHeight) {
				getSelectionCoordinates(x, y);
				if (lengthOf(x)>evaluationInterval+autocorelInterval) {
					dX=newArray(lengthOf(x)-evaluationInterval);
					dY=newArray(lengthOf(y)-evaluationInterval);
					S=newArray(lengthOf(x)-evaluationInterval);
					for (i=0; i<lengthOf(x)-evaluationInterval; i++) {
						dX[i]=(x[i+evaluationInterval]-x[i]);
						dY[i]=(y[i+evaluationInterval]-y[i]);
						S[i]=pow(dX[i]*dX[i]+dY[i]*dY[i], 0.5);
					}

					PC=newArray(lengthOf(dX)-autocorelInterval);
					SC=newArray(lengthOf(dX)-autocorelInterval);
					PCcount=0;
					SCcount=0;
					for (i=0; i<lengthOf(dX)-autocorelInterval; i++) {
						PC[PCcount]=(dX[i]*dX[i+autocorelInterval]+dY[i]*dY[i+autocorelInterval])/(S[i]*S[i+autocorelInterval]);
						if (isNaN(PC[PCcount])==false) {
							PCcount++;
						}
						SC[SCcount]=pow(S[i], 2)/(S[i]*S[i+autocorelInterval]);
						if (isNaN(SC[SCcount])==false && S[i]!=0 && S[i+autocorelInterval]!=0) {
							SCcount++;
						}
					}
					PC=Array.trim(PC, PCcount);
					SC=Array.trim(SC, SCcount);

					Array.getStatistics(PC, PCmin, PCmax, PCmean, PCstdev);
					if (isNaN(PCmean)==false) {
						allPC=Array.concat(allPC, PCmean);
					}
					Array.getStatistics(SC, SCmin, SCmax, SCmean, SCstdev);
					if (isNaN(SCmean)==false) {
						allSC=Array.concat(allSC, SCmean);
					}
				}
			}
		}
		Array.getStatistics(allPC, allPCmin, allPCmax, allPCmean, allPCstdev);
		Array.getStatistics(allSC, allSCmin, allSCmax, allSCmean, allSCstdev);
		print(autocorelInterval/10, allPCmean, allPCstdev, allPCstdev/pow(lengthOf(allPC), 0.5), allSCmean, allSCstdev, allSCstdev/pow(lengthOf(allSC), 0.5), lengthOf(allPC));
	}
}

function plotParticleTracks(centered) {
	newImage("Tracks", "32-bit black", width, height, 1);
	setColor(-1);
	fillRect(0, 0, width, height);
	for (i=0; i<particleTracksNumber; i++) {
		showProgress((i+1)/particleTracksNumber);
		showStatus("Plotting tracks...");
		start=extractTrack(i);
		if (selectionType()!=-1) {
			getSelectionBounds(tx, ty, tw, th);
			if (tx!=0 || ty!=0 || tw!=getWidth() || th!=getHeight) {
				getSelectionCoordinates(x, y);
				setColor(i);
				if (centered==false) {
					ox=0;
					oy=0;
				} else {
					ox=-x[0]+getWidth()/2;
					oy=-y[0]+getHeight()/2;
				}
				for (j=0; j<lengthOf(x)-1; j++) {
					drawLine(ox+x[j], oy+y[j], ox+x[j+1], oy+y[j+1]);
				}
			}
		}
	}
	run("Select None");
	setMinAndMax(0, particleTracksNumber);
	run("Spectrum");
	getLut(r, g, b);
	r[0]=0;
	g[0]=0;
	b[0]=0;
	setLut(r, g, b);
	return getImageID();
}

function printPopulationStats(file, useFile) {
	delim=",";
	str="";
	str+="Angle data\r\n";
	for (i=0; i<directionBins; i++) {
		angle=directionMax*i/directionBins;
		str+=""+angle+delim+directionNWeight[i]+delim+directionVWeight[i]+"\r\n";
	}
	if (useFile==true) {
		print(file, str);
	} else {
		print(str);
	}
}

function analysePath(path, overwrite, recursive) {
	files=getFileList(path);
	setBatchMode(true);
	for (i=0; i<lengthOf(files); i++) {
		if (File.isDirectory(path+files[i])==true && recursive==true) {
			analysePath(path+files[i], overwrite, recursive);
		} else {
			if ((endsWith(files[i], ".tif")==true || endsWith(files[i], ".czi")) && endsWith(files[i], "_tracks.tif")==false) {
				if (overwrite==true || File.exists(path+files[i]+"_trackstats.txt")==false) {
					if (endsWith(files[i], ".tif")==true) {
						print("Tif: "+path+files[i]);
						type="tif";
						name=substring(files[i], 0, lengthOf(files[i])-lengthOf(".tif"));
						series=1;
					} else if (endsWith(files[i], ".czi")==true) {
						print("Czi: "+path+files[i]);
						type="czi";
						name=substring(files[i], 0, lengthOf(files[i])-lengthOf(".czi"));
						run("Bio-Formats Macro Extensions");
						Ext.setId(path+files[i]);
						Ext.getSeriesCount(series);
					}
					//Setup output directories
					if (series==1) {
						outpath=path;
					} else {
						if (File.exists(path+name)==false) {
							File.makeDirectory(path+name);
						}
						outpath=path+name+File.separator();
					}
					//For each series
					for (k=0; k<series; k++) {
						if (type=="tif") {
							open(path+files[i]);
							title=files[i];
						} else if (type=="czi" && series==1) {
							run("Bio-Formats Windowless Importer", "open=["+path+files[i]+"] use_virtual_stack");
							vir=getImageID();
							run("Duplicate...", "duplicate");
							dup=getImageID();
							selectImage(vir);
							close();
							selectImage(dup);
							title=files[i];
						} else if (type=="czi") {
							Ext.setSeries(k);
							Ext.openImagePlus(path+files[i]);
							title=getTitle();
							print("Series "+k+": "+title);
						}
						if (overwrite==true || File.exists(outpath+title+"_trackstats.txt")==false) {
							//Checks for blank slices that would break the analysis
							run("Select None");
							for (j=0; j<nSlices(); j++) {
								setSlice(j+1);
								getRawStatistics(area, mean, min, max);
								if (min==max) {
									while (nSlices()>j) {
										setSlice(j+1);
										run("Delete Slice");
									}
								}
							}
							//Only continue if there is more than one slice
							if (nSlices()>1) {
								file=File.open(outpath+title+"_trackstats.txt");
								detectParticles(1, nSlices()+1);
								buildParticleTracks(file, true);
								plotParticleTracks(false);
								saveAs(outpath+title+"_tracks.tif");
								plotParticleTracks(true);
								saveAs(outpath+title+"_trackscentered.tif");
								close();
								File.close(file);
								file=File.open(outpath+title+"_dirstats.txt");
								printPopulationStats(file, true);
								File.close(file);
								file=File.open(outpath+title+"_trackstats-driftcorrected.txt");
								subtractDrift(1);
								recalculateStats(file, true);
								plotParticleTracks(false);
								saveAs(outpath+title+"_tracks-driftcorrected.tif");
								plotParticleTracks(true);
								saveAs(outpath+title+"_trackscentered-driftcorrected.tif");
								close();
								File.close(file);
								file=File.open(outpath+title+"_dirstats-driftcorrected.txt");
								printPopulationStats(file, true);
								File.close(file);
							}
						}
						close();
						//setBatchMode(false);
					}
				}
			}
		}
	}
}
