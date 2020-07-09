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

Dialog.create("Pre-Masking Threshold MIP creator");
Dialog.addNumber("Select neuron mask channel:", 1);
Dialog.addNumber("If using Median Filter, select radius (px)", 2);
Dialog.addNumber("Rolling ball radius (px), set to 0 for none", 20);
Dialog.show();

NeCh = Dialog.getNumber();
FSZ = Dialog.getNumber();
NeuronBGSubSZ = Dialog.getNumber();

input = getDirectory("Input directory");
raw = input + "Raw_data/";


run("Clear Results"); 
print("\\Clear");
print("Pre-Masking Threshold MIP Creator running:");
setBatchMode(true);
File.mkdir(input + "Threshold_List");		
File.mkdir(input + "Threshold_MIPs");

files = getFileList(raw);
if (files.length == 0){
	print("No images found!");
	print("***Please ensure Raw data is in a subfolder called 'Raw_data'.***");
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
		run("Bio-Formats Importer", "open=[" + raw + image + "] autoscale color_mode=Default specify_range concatenate_series open_all_series view=Hyperstack stack_order=XYCZT c_begin="+NeCh+" c_end="+NeCh+" c_step=1");
		ImageBitDepth = bitDepth();
		if (ImageBitDepth == 8) {
			run("16-bit");
		}
		rename(image);
		//rename as Tif - is their an easier way than this?
		source_Title= getTitle();
		Tif_Title = tiftitle(source_Title);
		rename("image"); 

		if (FSZ > 0) {
			run("Median...", "radius="+FSZ+" stack");
		}

		
		if (NeuronBGSubSZ > 0){
			run("Subtract Background...", "rolling="+ NeuronBGSubSZ +" stack");
		}
		run("Z Project...", "projection=[Max Intensity]");
		
		save(input + "Threshold_MIPs/" + Tif_Title);
		close();
		selectWindow("image");
		close();
	}
}
run("New... ", "name=[Thresholds] type=Table");
selectWindow("Thresholds");
saveAs("text", input +"Threshold_List/Thresholds.csv");
print("\\Update3:Finished!");
		
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
		