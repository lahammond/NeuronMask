// 3D Neuron Masking Tool
// Author: 	Luke Hammond (lh2881@columbia.edu)
// Cellular Imaging | Zuckerman Institute, Columbia University
// Date:	5th September 2017
//	
// This macro allows for masking over intracellular labelling using a neuron membrane mask.
// An automatic threhsold is set (Li) and this binary image is used to detect the neuron. 
// Subsequent ROIs are used to clear staining outside of the neuron.
// The output is a file showing the thresholded neuron and the cleaned up labels.
// 
// Potential to include an option for manual adjustment of threshold.

// Update 1/2018 v18: Included concatenate series for ND2 multiseries z-stacks. import multiseries, concatenate, hyperstack to stack, deinterleave based on # channels

//** add adjustment size detection for neuron and for inbetween branches!
//** fix up mask conversion step so it can work with both 16 and 8 bit

// Initialization
requires("1.41h");
starttime = getTime()
run("Options...", "iterations=8 count=1 black do=Nothing");
run("Set Measurements...", "fit redirect=None decimal=3");


// Options Window
Dialog.create("Neuron Masking Tool");
Dialog.addNumber("Select neuron mask channel:", 1);
Dialog.addCheckbox("Enhance neuron mask channel for segmentation? Use if processing non deconvolved images.",true);
SampleSelect = newArray("None", "Threshold List Provided", "Median", "Non-Local Means Denoising");
Dialog.addChoice("Filter before applying threshold:", SampleSelect, "Non-Local Means Denoising");
Dialog.addNumber("NL Means Denoising Sigma (40-150 decon, 150-300 Unsharp, test on sample region fist)", 15);
Dialog.addNumber("NL Means Denoising Smoothing Factor", 1);
Dialog.addNumber("If using Median Filter, select radius (px)", 1);
Dialog.addNumber("Rolling ball radius (px), set to 0 for none", 10);
Dialog.addCheckbox("Keep processed neuron channel?",true);

SampleSelect = newArray("Min(>25)", "Huang", "Li", "Moments", "Triangle", "Otsu");
Dialog.addChoice("Automatic Threshold:", SampleSelect, "Li");
Dialog.addNumber("Select minimum object size for neuron mask detection (4 for 0.3um, 15 for 0.1um):", 15);

Dialog.addMessage("-----------------------------------------------------------------------")
Dialog.addNumber("Select channel for masking:", 2);

Dialog.addNumber("Select second channel for masking (0 = none):", 0);


Dialog.addMessage("-----------------------------------------------------------------------")
Dialog.addCheckbox("Remove staining inbetween branches?",true);
Dialog.addNumber("Detection area for cleaning between branches", 20000);

Dialog.show();
NeCh = Dialog.getNumber();
Unsharp = Dialog.getCheckbox();
FType = Dialog.getChoice();
Sigma = Dialog.getNumber();
SmoothFactor = Dialog.getNumber();
FSZ = Dialog.getNumber();
NeuronBGSubSZ = Dialog.getNumber();
KeepNeu = Dialog.getCheckbox();
minSZ = Dialog.getNumber();
TType = Dialog.getChoice();
AnCh = Dialog.getNumber();
//BgSZ = Dialog.getNumber();
//BGSubONb = Dialog.getCheckbox();
AnCh2 = Dialog.getNumber();
//BgSZ2 = Dialog.getNumber();
//BGSubON2b = Dialog.getCheckbox();
InBTW = Dialog.getCheckbox();
InBTWSZ = Dialog.getNumber();

// Preparation
input = getDirectory("Input directory");
raw = input + "Raw_data/";


run("Clear Results"); 
print("\\Clear");
print("3D neuron mask running:");
setBatchMode(true);

files = getFileList(raw);

if (files.length == 0){
	print("No images found!");
	print("***Please ensure Raw data is in a subfolder called 'Raw_data'.***");
}

Array.sort( files );					
File.mkdir(input + "Processed");
if (KeepNeu == true); {
	File.mkdir(input + "Processed_Neuron_Channel");
}
File.mkdir(input + "MIPs");		
File.mkdir(input + "MIPs/Validation_MIPs");	
File.mkdir(input + "MIPs/Branches");
//File.mkdir(input + "MIPs/Terminal_Branches");
//File.mkdir(input + "MIPs/Axons");	
File.mkdir(input + "MIPs/Branches/ROIs");
//File.mkdir(input + "MIPs/Terminal_Branches/ROIs");
//File.mkdir(input + "MIPs/Axons/ROIs");

//iterate over all files
for(i=0; i<files.length; i++) {				

	// get the name of the current file
	image = files[i];							
	print("\\Update1: processing file " + (i+1) +" of " + files.length);
	looptimeS = getTime();
	// if the current file is a an image then:
	if (indexOf(image, ".nd2") != -1 || indexOf(image, ".tif") != -1) {			
		print("\\Update2: Importing Z-stack...");
				
		run("Bio-Formats Importer", "open=[" + raw + image + "] autoscale color_mode=Default concatenate_series open_all_series view=Hyperstack stack_order=XYCZT");
		rename(image);
		ImageBitDepth = bitDepth();
		if (ImageBitDepth == 8) {
			run("16-bit");
		}
		print("\\Update2: Z-stack imported. Detecting neuron and masking channel/s...");
		//rename as Tif - is their an easier way than this?
		source_Title= getTitle();
		Tif_Title = tiftitle(source_Title);
		source_ID = getImageID(); 
		
		getDimensions(width, height, channels, slices, frames);
		getPixelSize(unit, pW, pH);
		run("Properties...", "unit=pixel pixel_width=1.0000 pixel_height=1.0000");
		//middle = parseInt(d2s(slices/2, 1));

		//added steps for concatentation of multiseries TIF files 
		//NOTE Zseries >2GB don't open completely in FIJI if ND2. Must export as TIF from NIS first
		
		if (frames > 1) {
			setBatchMode(false);
			rename("image");
			run("Hyperstack to Stack");
			run("Deinterleave", "how=" + channels + "");
			selectWindow("image #" + NeCh + "");
			rename("C" + NeCh + "-" + image );
			selectWindow("image #" + AnCh + "");
			rename("C" + AnCh + "-" + image );
			setBatchMode(true);
		}
		else {
			run("Split Channels");
		}
		
		selectWindow("C" + NeCh + "-" + image );
		source_ID = getImageID();	

		// THRESHOLD
		run("Duplicate...", "title=Neuron_Raw duplicate");
		selectWindow("C" + NeCh + "-" + image ); 
		//setSlice(middle);
		if (NeCh==AnCh) {
			run("Duplicate...", "title=Neuron_Mask duplicate"); 
			selectWindow("Neuron_Mask");
		}
		
		//run("Gray Scale Attribute Filtering 3D", "operation=Opening attribute=Volume min=50 connectivity=6");
		//run("Remove Outliers...", "radius=5 threshold=0 which=Dark stack");
		
		//run("Despeckle", "stack");
		//run("Enhance Contrast...", "saturated=0.3 process_all use");


		// unsharp to enhance neuron mask if no deconvolution used.
		if (Unsharp==true) {
			run("Unsharp Mask...", "radius=2 mask=0.70 stack");
		}
		
		if (FType=="Median") {
			run("Median...", "radius=FSZ stack");
		}
		
		if (FType=="Non-Local Means Denoising"){
			getDimensions(width, height, channels, slices, frames);
			for (l=1; l<slices+1; l++) {
			setSlice(l);
			run("Non-local Means Denoising", "sigma="+ Sigma +" smoothing_factor="+ SmoothFactor +" slice");
			}
		}

		run("Mean 3D...", "x=1 y=1 z=1");
		
		if (NeuronBGSubSZ > 0){
			run("Subtract Background...", "rolling="+ NeuronBGSubSZ +" stack");
		}
		
		//if (Unsharp==true) {
		//	run("Subtract...", "value=100 stack");
		//}
		
		rename("Neuron_Mask");		
		if (KeepNeu == true) {
			save(input + "Processed_Neuron_Channel/" + Tif_Title);
		}

		
		
		//title = "WaitForUser";
		//msg = "If necessary adjustment the threshold, then click \"OK\".";
		//waitForUser(title, msg);
		//new clean up
		
		run("Remove Outliers...", "radius=2 threshold=0 which=Bright stack");

		if (FType=="Threshold List Provided"){
			open(input + "Threshold_List/Thresholds.csv");
			Threshold = getResult("C"+(i+1)+"");
			setThreshold(Threshold, 65535);
			setOption("BlackBackground", true);
			run("Convert to Mask", "method=Default background=Dark black");
			print("\\Update2: Setting threshold of Image "+(i+1)+ "at "+ Threshold+".");
		
		} else {
			if (TType == "Min"){
				setThreshold(25, 65535);
				setOption("BlackBackground", true);
				run("Convert to Mask", "method=Default background=Dark black");
			}else {
			run("Auto Threshold", "method=" + TType + " white stack use_stack_histogram");
			}
		}

		
		//setAutoThreshold(TType+" dark stack");
		//run("Convert to Mask", "method="+TType+" background=Dark black");
		//Fill holes in processes
		run("Remove Outliers...", "radius=2 threshold=0 which=Bright stack");
		
		run("Options...", "iterations=2 count=1 black do=Nothing");
		run("Dilate", "stack");
		run("Erode", "stack");
		
		//run("Make Binary", "method=Li background=Dark");
		run("Duplicate...", "title=inverted duplicate"); 
		run("Invert", "stack");

		// Peform BG Sub before clearing
		selectWindow("C" + AnCh + "-" + image );
		rename("AnCh_Clean");
		run("Duplicate...", "title=AnChRaw duplicate"); 
		if (AnCh2 > 0) {
			selectWindow("C" + AnCh2 + "-" + image );
			rename("AnCh2_Clean");
			run("Duplicate...", "title=AnChRaw2 duplicate"); 
		}
		
		// Removed BG SUBTRACT STEPS - as we want raw data for later colocalisation analysis
		//if (BgSZ > 0) {
		//	 selectWindow("AnCh_Clean");
		//	 run("Subtract Background...", "rolling=BgSZ stack");
		//	}
		//if (AnCh2 > 0 && BgSZ2 > 0) {
		//	 selectWindow("AnCh2_Clean");
		//	 run("Subtract Background...", "rolling=BgSZ2 stack");
		//	}

		//DETECT NEURON and clear outside
		for (j=1; j<slices+1; j++) {
			//selectWindow("C" + NeCh + "-" + image );
			selectWindow("Neuron_Mask");
			//run("Select None");
			setSlice(j);
			// cleanup ROIs
			cleanupROI();
			// Analyse particles
			run("Analyze Particles...", "size="+ minSZ +"-Infinity add slice");
			selectWindow("AnCh_Clean");
			setSlice(j);
						
			// Perform Clearing	of branches
			CountROIsub=0;
			CountROIsub=roiManager("count"); 
			if (CountROIsub==0) {
				// clear channel of interest
				selectWindow("AnCh_Clean");
				clearslice();
				// clear second channel of interest
				if (AnCh2 > 0) {
					selectWindow("AnCh2_Clean");
					setSlice(j);
					clearslice();
				}
				selectWindow("Neuron_Mask");
				clearslice();
			} else if (CountROIsub==1) {
				// clear channel of interest
				selectWindow("AnCh_Clean");
				clearOutsideROI0();
				// clear second channel of interest
				if (AnCh2 > 0) {
					selectWindow("AnCh2_Clean");
					setSlice(j);
					clearOutsideROI0();
				}
				// clear neuron channel
				selectWindow("Neuron_Mask");
				clearOutsideROI0();
				// clean up ROIs
				deleteROIs();
				CountROIsub=0;
			} else if (CountROIsub>1) {
				// clear channel of interest
				selectWindow("AnCh_Clean");
				clearOutsideROIsubarry();
				// clear second channel of interest
				if (AnCh2 > 0) {
					selectWindow("AnCh2_Clean");
					setSlice(j);
					clearOutsideROIsubarry();
				}
				// clear neuron channel
				selectWindow("Neuron_Mask");
				clearOutsideROIsubarry();
				// clean up ROIs
				deleteROIs();
				ROIarraysub=newArray(0);
				CountROIsub=0;
				} 
		}
		//Clear the places between branches
		if (InBTW==true) {		
			for (k=1; k<slices+1; k++) {
				selectWindow("inverted");
				setSlice(k);
				run("Analyze Particles...", "size=0-"+ InBTWSZ +" add slice");
				selectWindow("AnCh_Clean");
				setSlice(k);
				run("Select None");
				CountROIsub=roiManager("count"); 
				if (CountROIsub==1) {
					clearROI0();
					if (AnCh2 > 0) {
						selectWindow("AnCh2_Clean");
						setSlice(k);
						clearROI0();
					}
					deleteROIs();
					CountROIsub=0;
				} else if (CountROIsub>1) {
					clearROIsubarray();
					if (AnCh2 > 0) {
						selectWindow("AnCh2_Clean");
						setSlice(k);
						clearROIsubarray();
					}
					deleteROIs();
					ROIarraysub=newArray(0);
					CountROIsub=0;
				}
			}
		}
		print("\\Update2: Channel masking complete... almost ready.");

//** need to modify so that 8bit is 8bit mask and 16bit converts to 16bit mask		
		// Convert 8bit mask to 16bit

		selectWindow("Neuron_Mask");
		run("16-bit");
		setMinAndMax(0, 255);
		run("Apply LUT", "stack");
		selectWindow("inverted");
		close();

		if (AnCh2==0) {
			run("Merge Channels...", "c1=Neuron_Mask c2=AnCh_Clean create keep");
			//run("Merge Channels...", "c1=AnCh_Clean c2=Neuron_Mask c4=Neuron_Raw create");
			Stack.setChannel(1);
			run("Grays");
			Stack.setChannel(2);
			run("Green");
			
		} else {	
			//run("Merge Channels...", "c1=AnChRaw C2=AnCh_Clean c3=AnChRaw2 C4=AnCh2_Clean c5=Neuron_Mask c6=Neuron_Raw create");
			run("Merge Channels...", "c1=Neuron_Mask c2=AnCh_Clean c3=AnCh2_Clean create keep");
			Stack.setChannel(1);
			run("Grays");
			Stack.setChannel(2);
			run("Green");
			Stack.setChannel(3);
			run("Red");
		}
		
		run("Properties...", "unit="+unit+" pixel_width="+pW+" pixel_height="+pH+"");
					
		//Save Merged Processed Data
		save(input + "Processed/" + Tif_Title);
		close();

		// Create MIP for validation and Distance Mapping
		selectWindow("Neuron_Mask");
		run("Z Project...", "projection=[Max Intensity]");
		run("8-bit");
		//Save Neuron MIP as single channel
		save(input + "MIPs/Branches/" + Tif_Title);
		//save(input + "MIPs/Terminal_Branches/" + Tif_Title); Why did I take this out because I still need it, putting it back in.
		//save(input + "MIPs/Terminal_Branches/" + Tif_Title);
		//save(input + "MIPs/Axons/" + Tif_Title);
		// then produce validation MIP
		run("16-bit");
		setMinAndMax(0, 255);
		run("Apply LUT");
		selectWindow("AnCh_Clean");
		run("Z Project...", "projection=[Max Intensity]");
		selectWindow("AnChRaw");
		run("Z Project...", "projection=[Max Intensity]");
		selectWindow("Neuron_Raw");
		run("Z Project...", "projection=[Max Intensity]");
		if (AnCh2==0) {
			run("Merge Channels...", "c1=MAX_Neuron_Mask c2=MAX_Neuron_Raw c3=MAX_AnCh_Clean c4=MAX_AnChRaw  create");
			Stack.setChannel(1);
			run("Green");
			Stack.setChannel(2);
			run("Magenta");
			//run("Enhance Contrast...", "saturated=0.3 equalize");
			//resetMinAndMax();
			Stack.setChannel(3);
			run("Green");
			Stack.setChannel(4);
			run("Red");
			//run("Enhance Contrast...", "saturated=0.3 equalize");
			//resetMinAndMax();		
		} else {
			selectWindow("AnCh2_Clean");
			run("Z Project...", "projection=[Max Intensity]");
			selectWindow("AnChRaw2");
			run("Z Project...", "projection=[Max Intensity]");
			run("Merge Channels...", "c1=MAX_Neuron_Mask c2=MAX_Neuron_Raw c3=MAX_AnCh_Clean c4=MAX_AnChRaw c5=MAX_AnCh2_Clean c6=MAX_AnChRaw2  create");
			Stack.setChannel(1);
			run("Green");
			Stack.setChannel(2);
			run("Magenta");
			//run("Enhance Contrast...", "saturated=0.3 equalize");
			//resetMinAndMax();
			Stack.setChannel(3);
			run("Green");
			//run("Enhance Contrast...", "saturated=0.3 equalize");
			//resetMinAndMax();
			Stack.setChannel(4);
			run("Red");
			//run("Enhance Contrast...", "saturated=0.3 equalize");
			//resetMinAndMax();
			Stack.setChannel(5);
			run("Green");
			Stack.setChannel(6);
			run("Red");
		}
		//Save MIPs
		save(input + "MIPs/Validation_MIPs/" + Tif_Title);
						
		close("\\Others");		
		close();
	}
looptimeE=getTime();
loopdif = (looptimeE-looptimeS)/1000;
print("\\Update3: Total processing time for this image = "+ loopdif + " seconds.");
}
// Get time
endtime = getTime();
dif = ((endtime-starttime)/1000)/60;
print("\\Update4: Total processing time = "+ dif + " minutes for "+ i +" files.");
print("\\Update5: Total processing time = "+ (dif/i) + " minutes per file.");
print("\\Update6: ---");
print("\\Update7: Settings:");
print("\\Update8: Filter: "+ FType +" Sigma: " + Sigma + " Smoothing Factor: "+ SmoothFactor +".");
print("\\Update9: BGSub Neuron: "+ NeuronBGSubSZ +" BGSub compartments: " + NeuronBGSubSZ + " Min size for neuron mask:"+ minSZ +" Inbetween Branch Size: "+ InBTWSZ +".");

selectWindow("Log");
saveAs("txt", input+"/3D_Neuron_Mask_Log.txt");

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
