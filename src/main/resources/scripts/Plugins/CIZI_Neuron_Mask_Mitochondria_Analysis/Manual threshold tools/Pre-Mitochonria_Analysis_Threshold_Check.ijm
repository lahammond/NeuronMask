// Pre-Masking Threshold Checking Tool
// Author: 	Luke Hammond (lh2881@columbia.edu)
// Cellular Imaging | Zuckerman Institute, Columbia University
// Date:	22nd February 2018
//	
// This macro allows opens an image and creates a MIP for checking the threshold
// save the threshold value in a csv file, in image order, seperated by ","

// Initialization
requires("1.41h");
starttime = getTime()
run("Options...", "iterations=8 count=1 black do=Nothing");
run("Set Measurements...", "fit redirect=None decimal=3");

Dialog.create("Pre-Object Analysis Threshold MIP creator");
SampleSelect = newArray("Axons", "Primary Branches");
Dialog.addChoice("Regions to create MIPs from:", SampleSelect, "Primary Branches");
Dialog.addNumber("Number of masked channels:", 2);
//Dialog.addNumber("Median filter radius (px)", 1);
Dialog.addNumber("Rolling ball radius (px), set to 0 for none", 0);
Dialog.addCheckbox("Use iterative outlier filter (alternative to bg subtract)?",true);
Dialog.show();

Region = Dialog.getChoice();
Channels = Dialog.getNumber();
//FSZ = Dialog.getNumber();
BGSubSZ = Dialog.getNumber();
ROF = Dialog.getCheckbox();

input = getDirectory("Input directory");

if (Region == "Primary Branches") {
	regionfolder = input + "Extracted_Regions/Primary_Branches/";
} 
if (Region == "Axons") {
	regionfolder = input + "Extracted_Regions/Axons/";
}


run("Clear Results"); 
print("\\Clear");
print("Pre-Object Analysis Threshold MIP Creator running:");
setBatchMode(true);
File.mkdir(input + "Region_Threshold_List");		
File.mkdir(input + "Region_Threshold_MIPs");

files = getFileList(regionfolder);
if (files.length == 0){
	print("No images found!");
	//print("***Please ensure Raw data is in a subfolder called 'Raw_data'.***");
}
Array.sort( files );			


//iterate over all files
for(i=0; i<files.length; i++) {				

	// get the name of the current file
	image = files[i];							
	print("\\Update1: processing file " + (i+1) +" of " + files.length);
	looptimeS = getTime();
	// if the current file is a an image then:
	if (indexOf(image, ".nd2") != -1 || indexOf(image, ".tif") != -1) {			
		print("\\Update2: Importing Z-stack...");
		run("Bio-Formats Importer", "open=[" + regionfolder + image + "] autoscale color_mode=Default specify_range concatenate_series open_all_series view=Hyperstack stack_order=XYCZT c_begin=2 c_end="+(1+Channels)+" c_step=1");
		ImageBitDepth = bitDepth();
		if (ImageBitDepth == 8) {
			run("16-bit");
		}
		rename(image);
		//rename as Tif - is their an easier way than this?
		source_Title= getTitle();
		Tif_Title = tiftitle(source_Title);
		getDimensions(width, height, channels, slices, frames);
		rename("image"); 
		if (channels == 1) {
			rename("C1-image");
		}

		//if (FSZ > 0) {
		//	run("Median...", "radius="+FSZ+" stack");
		//}
		
		//run("Median...", "radius=1 stack");
		
		if (BGSubSZ > 0){
			run("Subtract Background...", "rolling="+ BGSubSZ +" stack");
		}

		if (ROF == true){
			if (channels == 1) {
				run("Duplicate...", "title=C1raw duplicate");
				RemoveOutliersFilter("C1-image");
				run("Merge Channels...", "c1=C1raw c2=C1-image create");
				rename("image"); 
				
			}
			if (channels ==2) {
				run("Split Channels");
				selectWindow("C1-image");
				run("Duplicate...", "title=C1raw duplicate");
				selectWindow("C2-image");
				run("Duplicate...", "title=C2raw duplicate");
				RemoveOutliersFilter("C1-image");
				RemoveOutliersFilter("C2-image");
				run("Merge Channels...", "c1=C1raw c2=C1-image c3=C2raw c4=C2-image create");
				rename("image"); 
			}
		}		
			
		
		run("Z Project...", "projection=[Max Intensity]");
		
		save(input + "Region_Threshold_MIPs/" + Tif_Title);
		close();
		selectWindow("image");
		close();
	}
}
run("New... ", "name=[Thresholds] type=Table");
selectWindow("Thresholds");
saveAs("text", input +"Region_Threshold_List/Region_ThresholdsCh1.csv");
if (Channels == 2) {
	saveAs("text", input +"Region_Threshold_List/Region_ThresholdsCh2.csv");
}
print("\\Update3:Finished!");




function RemoveOutliersFilter(imagename) {
	selectWindow(imagename);
	getDimensions(w2, h2, c2, slices, f2);
	rename("ROFimage");
	run("Duplicate...", "title=bgstack duplicate");
	run("Z Project...", "projection=[Max Intensity]");
	for(i=0; i<10; i++) {	
		run("Remove Outliers...", "radius=5 threshold=0 which=Bright");
		//run("Morphological Filters", "operation=Dilation element=Disk radius=3");
	}
	run("Morphological Filters", "operation=Dilation element=Disk radius=10");
	rename("morph");
	selectWindow("MAX_bgstack");
	close();
	selectWindow("morph");
	rename("MAX_bgstack");
	
	for (i=0; i<slices; i++) {
		selectWindow("MAX_bgstack");
		run("Select All");
		run("Copy");
		selectWindow("bgstack");
		setSlice(i+1);
		run("Paste");
	}
	selectWindow("MAX_bgstack");
	close();
	
	imageCalculator("Subtract create stack", "ROFimage","bgstack");
	selectWindow("bgstack");
	close();
	selectWindow("ROFimage");
	close();
	selectWindow("Result of ROFimage");
	rename(imagename);
}		


		
function tiftitle(imagename){
	nl=lengthOf(imagename);
	nl2=nl-4;
	Sub_Title=substring(imagename,0,nl2);
	Tif_Title=Sub_Title+".tif";
	return Tif_Title;
}

function cleanupROI() {
	CountROImain=roiManager("count"); 
		if (CountROImain == 1) {
			roiManager("delete");
			CountROImain=0;
		} else if (CountROImain > 1) {
			ROIarrayMain=newArray(CountROImain); 
			for(n=0; n<CountROImain;n++) { 
	       		ROIarrayMain[n] = n; 
				} 
			roiManager("Select", ROIarrayMain);
			roiManager("Combine");
			roiManager("Delete");
			ROIarrayMain=newArray(0);
			CountROImain=0;
		}		
}

function clearslice() {
	run("Select All");
	run("Clear", "slice");
	run("Select None");
}

function clearROI0() {
	roiManager("Select", 0);
	run("Clear", "slice");
}
function deleteROIs() {
	roiManager("Deselect");
	roiManager("Delete");
	run("Select None");
}

function clearOutsideROI0() {
	roiManager("Select", 0);
	run("Clear Outside", "slice");
}

function clearROIsubarray() {
	CountROIsub=roiManager("count"); 
	ROIarraysub=newArray(CountROIsub); 
	for(m=0; m<CountROIsub;m++) { 
	   ROIarraysub[m] = m; 
	} 
	roiManager("Select", ROIarraysub);
	roiManager("Combine");
	run("Clear", "slice");
}

function clearOutsideROIsubarry() {
	CountROIsub=roiManager("count"); 
	ROIarraysub=newArray(CountROIsub); 
	for(m=0; m<CountROIsub;m++) { 
		ROIarraysub[m] = m; 
	} 
	roiManager("Select", ROIarraysub);
	roiManager("Combine");
	run("Clear Outside", "slice");
}
		