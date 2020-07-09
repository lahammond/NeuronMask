//Neuron Branch Object Analysis tool
//
// Author: 	Luke Hammond
// Cellular Imaging | Zuckerman Institute, Columbia University
// Version: 0.5
// Date:	10th November 2017
//
// For endosome anlaysis in neuron branches
// Expects 3 channel image (distance, cleaned and raw data for 1 channel masked neuron) or
// 5 channel image (distance, cleaned Ch1, raw Ch1, cleaned Ch2, raw Ch2)
//
// Update V0.3 (16th Nov) - Instead of global treshold, subtract minimum background value then
//							run a local otsu threshold to isolate objects
// Update v.05 (22nd Feb) - MAX IP for validating objects rather than 3d Stack - (consider creating a max ip of coloured objects!)
// Update v.06 (19th Jun 18) - updated for branch analysis of mitochondria


// Initialization
requires("1.41h");
starttime = getTime();
run("Options...", "iterations=3 count=1 black do=Nothing");
run("Set Measurements...", "fit redirect=None decimal=3");
run("Colors...", "foreground=white background=black selection=yellow");
run("Clear Results");
cleanupROI();

// menu
Dialog.create("3D Neuron Branch Object Analysis");

Dialog.addCheckbox("Detect objects using provided threshold list/s?",false);
Dialog.addNumber("Number of masked channels:", 1);
Dialog.addMessage("-")
Dialog.addMessage("If not using a provided threhsold list please use approxmate cytosol levels below.")
Dialog.addNumber("Channel 1 approximate cytosol background intensity value:", 100);
Dialog.addNumber("Channel 2 approximate cytosol background intensity value:", 100);
Dialog.addMessage("-")
Dialog.addNumber("Rolling ball radius (px)", 0);
Dialog.addCheckbox("Use iterative outlier filter (alternative to bg subtract)?",true);



//Dialog.addNumber("Number of channels to analyze", 1);
// add text stating that the Neuron_Regions directory should be selected as the input directory


Dialog.show();
ThList = Dialog.getCheckbox();
NoOfChannels = Dialog.getNumber();
Thr = Dialog.getNumber();
Thr2 = Dialog.getNumber();
BGSZ = Dialog.getNumber();
ROF = Dialog.getCheckbox();

// Preparation
input = getDirectory("Input directory:");
print("\\Clear");
print("3D Neuron Branch Object Analysis running:");
setBatchMode(true);
Axon = false;
PB = false;
TB = false;

// SECOND PASS ANALYSIS OPTION - DETECT INTEGRIN, THEN FILTER OUT OBJECTS LESS THAN 2SD BELOW MEAN
SecondPassC1 = false;
SecondPassC2 = false;




if (File.exists(input + "Extracted_Regions/Branches/"))
	BranchFolder = true;
//if (File.exists(input + "Extracted_Regions/Primary_Branches/"))
//	PB = true;
//if (File.exists(input + "Extracted_Regions/Terminal_Branches/"))
//	TB = true;

if (File.exists(input + "Analyzed/") == 0) {
	File.mkdir(input + "Analyzed");
}


if (BranchFolder == true){
	as = getTime();
	if (File.exists(input + "Analyzed/Branches") == 0) {
		File.mkdir(input + "Analyzed/Branches");
	}
	File.mkdir(input + "Analyzed/Branches/Counts");
	File.mkdir(input + "Analyzed/Branches/Coloc");
	Branches = getFileList(input + "Extracted_Regions/Branches");
	for(i=0; i<Branches.length; i++) {	
		Branch = Branches[i];
		cleanupROI();
		if (ThList == true){
			open(input + "Region_Threshold_List/Region_ThresholdsCh1.csv");
			ThresholdCh1 = getResult("C"+(i+1)+"");
			if (NoOfChannels == 2) {
				open(input + "Region_Threshold_List/Region_ThresholdsCh2.csv");
				ThresholdCh2 = getResult("C"+(i+1)+"");
			}
		closewindow("Results");
		} else {
			ThresholdCh1 = Thr;
			ThresholdCh2 = Thr2;	
		}
		print("\\Update1: Processing cell " + (i+1) + " of " + (Branches.length) + ".");
		branchextractor(input + "Extracted_Regions/Branches/",input + "Analyzed/Branches/", Branch, ThresholdCh1, ThresholdCh2);	
					
	}
	ae = getTime();
	atime = (ae-as)/1000;
	print("\\Update1: Processing time = ", (atime/60), "minutes. ", (atime/(i-1)), "seconds per cell");
} 
/*
if (PB == true){
	ps = getTime();
	if (File.exists(input + "Analyzed/Primary_Branches") == 0) {
		File.mkdir(input + "Analyzed/Primary_Branches");
	}
	File.mkdir(input + "Analyzed/Primary_Branches/Counts");
	File.mkdir(input + "Analyzed/Primary_Branches/Coloc");
	branches = getFileList(input + "Extracted_Regions/Primary_Branches");
	for(i=0; i<branches.length; i++) {	
		Branch = branches[i];
		cleanupROI();
		if (ThList == true){
			open(input + "Region_Threshold_List/Region_ThresholdsCh1.csv");
			ThresholdCh1 = getResult("C"+(i+1)+"");
			if (NoOfChannels == 2) {
				open(input + "Region_Threshold_List/Region_ThresholdsCh2.csv");
				ThresholdCh2 = getResult("C"+(i+1)+"");
			}
		closewindow("Results");
		} else {
			ThresholdCh1 = Thr;
			ThresholdCh2 = Thr2;
		}
		print("\\Update2: Primary branches detected, processing branch " + (i+1) + " of " + (branches.length) + ".");
		branchextractor(input + "Extracted_Regions/Primary_Branches/",input + "Analyzed/Primary_Branches/", Branch, ThresholdCh1, ThresholdCh2);		
	}
	pe = getTime();
	ptime = (pe-ps)/1000;
	print("\\Update2: Primary branch processing time = ", (ptime/60), "minutes. ", (ptime/(i-1)), "seconds per axon");
} else {
	print("\\Update2: No primary branches to process. Looking for terminal branches...");
}

if (TB == true){
	ts = getTime();
	if (File.exists(input + "Analyzed/Terminal_Branches") == 0) {
		File.mkdir(input + "Analyzed/Terminal_Branches");
	}
	File.mkdir(input + "Analyzed/Terminal_Branches/Counts");
	File.mkdir(input + "Analyzed/Terminal_Branches/Coloc");
	branches = getFileList(input + "Extracted_Regions/Terminal_Branches");
	for(i=0; i<branches.length; i++) {	
		Branch = branches[i];
		cleanupROI();
		if (ThList == true){
			open(input + "Region_Threshold_List/Region_ThresholdsCh1.csv");
			ThresholdCh1 = getResult("C"+(i+1)+"");
			if (NoOfChannels == 2) {
				open(input + "Region_Threshold_List/Region_ThresholdsCh2.csv");
				ThresholdCh2 = getResult("C"+(i+1)+"");
			}
		} else {
			ThresholdCh1 = Thr;
			ThresholdCh2 = Thr2;
		}
		print("\\Update3: Terminal branches detected, processing branch " + (i+1) + " of " + (branches.length) + ".");
		branchextractor(input + "Extracted_Regions/Terminal_Branches/",input + "Analyzed/Terminal_Branches/", Branch, ThresholdCh1, ThresholdCh2);		
	}
	te = getTime();
	ttime = (te-ts)/1000;
	print("\\Update3: Terminal branch processing time = ", (ttime/60), "minutes. ", (ttime/(i-1)), "seconds per axon");	
} else {
	print("\\Update3: No terminal branches to process.");
}
*/
endtime = getTime();
dif = (endtime-starttime)/1000;
print("\\Update4: Total processing time = ", (dif/60), "minutes.");




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
function tiftitlespace(imagename){
	nl=lengthOf(imagename);
	nl2=nl-5;
	Sub_Title=substring(imagename,0,nl2);
	Tif_Title=Sub_Title+".tif";
	return Tif_Title;
}

function branchextractor(directory, analysisdir, Branch, ThresholdCh1, ThresholdCh2){
	run("Bio-Formats Importer", "open=[" + directory + Branch +"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

	Output_title = subtitle(Branch);
	
	rename("Branch");
	getDimensions(width, height, channels, slices, frames);
	getPixelSize(MainUnit, MainW, MainH);
	
	run("Split Channels");
	selectWindow("C2-Branch");
	run("Duplicate...", "title=C2-BranchRaw duplicate");
	selectWindow("C2-Branch");
	
	// Image processing before endosome analysis. If not using decon then several steps required.
	if (ThList == false) {
		if (BGSZ > 0){
			run("Subtract Background...", "rolling="+ BGSZ +" stack");
		}
		if (ROF == true){
			RemoveOutliersFilter("C2-Branch");
			rename("C2-Branch");
			print("Setting threshold for channel 1 to: "+ThresholdCh1);
			setThreshold(ThresholdCh1, 65535);
			run("Convert to Mask", "method=Default background=Dark black");
		} else {
			print("Setting the minimum intensity for channel 1 to: "+ThresholdCh1);
			run("Min...", "value="+ ThresholdCh1 +" stack");
			setMinAndMax(0, 4095);
			call("ij.ImagePlus.setDefault16bitRange", 12);
			run("8-bit");
			// otsu started to behave strangly on FIJI update 16/1/18 - replaced with mean
			run("Auto Local Threshold", "method=Phansalkar radius=10 parameter_1=0 parameter_2=0 white stack");
		}
	} else {
		//run("Median...", "radius=1");
		if (BGSZ > 0){
			run("Subtract Background...", "rolling="+ BGSZ +" stack");
		}
		if (ROF == true){
			RemoveOutliersFilter("C2-Branch");
			rename("C2-Branch");
			print("Setting threshold for channel 1 to: "+ThresholdCh1);
			setThreshold(ThresholdCh1, 65535);
			run("Convert to Mask", "method=Default background=Dark black");
		} else {
			print("Setting the minimum intensity for channel 1 to: "+ThresholdCh1);
			run("Min...", "value="+ ThresholdCh1 +" stack");
			setMinAndMax(0, 4095);
			call("ij.ImagePlus.setDefault16bitRange", 12);
			run("8-bit");
			// otsu started to behave strangly on FIJI update 16/1/18 - replaced with mean
			run("Auto Local Threshold", "method=Phansalkar radius=10 parameter_1=0 parameter_2=0 white stack");
		}
	}

	if (SecondPassC1 == true) {
		run("Duplicate...", "title=C2-BranchSP duplicate");
		selectWindow("C2-Branch");
	}

	//settings to analyse raw intensities of channel 1 (C3)
	run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels" +
		" integrated_density mean_gray_value std_dev_gray_value minimum_gray_value" +
		" centroid centre_of_mass show_masked_image_(redirection_requiered)" +
		" dots_size=5 font_size=10" +
		" store_results_within_a_table_named_after_the_image_(macro_friendly)" + 
		" redirect_to=C2-BranchRaw");
	selectWindow("C2-Branch");
	getStatistics(area, mean, min, max, std, histogram);
	if (max < 128) threshold = max;
  	else threshold = 128;
	run("3D Objects Counter", "threshold="+threshold+" slice=63 min.=2 max.=33343088 objects statistics summary");
	
	selectWindow("Statistics for C2-Branch redirect to C2-BranchRaw");
	saveAs("Results", analysisdir +"/Counts/"+ Output_title + "_Ch1_Objects.csv");
	closewindow(Branch + "_Ch1_Objects.csv");


	
	// Read values from this table
	open(analysisdir +"/Counts/"+ Output_title + "_Ch1_Objects.csv");
	
	NumObjects = (nResults);
	ObjectVols = newArray(NumObjects);
	ObjectSurfs = newArray(NumObjects);
	ObjectVoxels = newArray(NumObjects);
	ObjectIntDen = newArray(NumObjects);
	ObjectMeans = newArray(NumObjects);
	for(j=0; j<NumObjects; j++) {
		ObjectVols[j] = getResult("Volume (micron^3)", j);
		ObjectSurfs[j] = getResult("Surface (micron^2)", j);
		ObjectVoxels[j] = getResult("Nb of obj. voxels", j);
		ObjectIntDen[j] = getResult("IntDen", j);
		ObjectMeans[j] = getResult("Mean", j);
				
	}
	

	//For validation
	selectWindow("Masked image for C2-Branch redirect to C2-BranchRaw");
	rename("C2-Mask");
	run("16-bit");
	selectWindow("Objects map of C2-Branch redirect to C2-BranchRaw");
	rename("ObjectMap");
	run("16-bit");
	run("Merge Channels...", "c1=C2-Mask c2=C2-BranchRaw c3=ObjectMap create keep ignore");
	Stack.setChannel(3);
	run("3-3-2 RGB");
	rename("MergeCh1");
	run("Z Project...", "projection=[Max Intensity]");
	saveAs("Tiff", analysisdir + Output_title + "_Ch1_Objects_Overlay.tif");
	close();
	closewindow("MergeCh1");
	closewindow("C2-Mask");
	closewindow("C2-BranchRaw");

	// switched to object map as C2-Mask was giving different objects out!
	selectWindow("ObjectMap");
	rename("C2-Mask");
	

	// settings to analyse distance (in pixels) to soma of channel 1 (C1)
	run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels mean_gray_value dots_size=5 font_size=10" +
		" store_results_within_a_table_named_after_the_image_(macro_friendly)" +
		" redirect_to=C1-Branch");
	selectWindow("C2-Branch");

	getStatistics(area, mean, min, max, std, histogram);
	if (max < 128) threshold = max;
  	else threshold = 128;
	
	run("3D Objects Counter", "threshold="+threshold+" slice=63 min.=2 max.=33343088 statistics summary");
	
	selectWindow("Statistics for C2-Branch redirect to C1-Branch");
	saveAs("Results", analysisdir +"/Counts/"+ Output_title + "_Ch1_Objects.csv");
	closewindow(Output_title + "_Ch1_Objects.csv");
	closewindow("C2-Branch");

	// Read values from this table
	open(analysisdir +"/Counts/"+ Output_title + "_Ch1_Objects.csv");
	
	ObjectDistancePx = newArray(NumObjects);
	ObjectDistanceUM = newArray(NumObjects);
	for(j=0; j<NumObjects; j++) {
		ObjectDistancePx[j] = getResult("Mean", j);
		ObjectDistanceUM[j] = (ObjectDistancePx[j]*MainW);
				
	}
	
	//Create a annotated table --
	title1 = "Annotated_Objects"; 
	title2 = "["+title1+"]"; 
	f=title2; 
	run("New... ", "name="+title2+" type=Table"); 
	print(f,"\\Headings:Volume (micron^3)\tSurface Area (micron^2)\tTotal Voxels\tIntegrated Density\tMean Intensity\tDistance from soma (um)\tDistance from soma (px)"); 
	for(j=0; j<NumObjects; j++){ 
		print(f,ObjectVols[j]+"\t"+ObjectSurfs[j]+"\t"+ObjectVoxels[j]+"\t"+ObjectIntDen[j]+"\t"+ObjectMeans[j]+"\t"+ObjectDistanceUM[j]+"\t"+ObjectDistancePx[j]);
	} 
	selectWindow(title1);

	run("Text...", "save=["+ analysisdir +"/Counts/"+ Output_title + "_Ch1_Objects.csv]");
	closewindow(title1);


	
		

	if (channels == 3){
		//settings to analyse raw intensities of channel 1 (C3)
		selectWindow("C3-Branch");
		run("Duplicate...", "title=C3-BranchRaw duplicate");
		selectWindow("C3-Branch");
		if (ThList == false) {
			//run("Median...", "radius=1");
			if (BGSZ > 0){
				run("Subtract Background...", "rolling="+ BGSZ +" stack");
			}
			if (ROF == true){
				RemoveOutliersFilter("C3-Branch");
				rename("C3-Branch");
				print("Setting threshold for channel 2 to: "+ThresholdCh2);
				setThreshold(ThresholdCh2, 65535);
				run("Convert to Mask", "method=Default background=Dark black");
			} else {
				print("Setting the minimum intensity for channel 1 to: "+ThresholdCh2);
				run("Min...", "value="+ ThresholdCh2 +" stack");
				setMinAndMax(0, 4095);
				call("ij.ImagePlus.setDefault16bitRange", 12);
				run("8-bit");
				// otsu started to behave strangly on FIJI update 16/1/18 - replaced with mean
				run("Auto Local Threshold", "method=Phansalkar radius=10 parameter_1=0 parameter_2=0 white stack");
			}
		} else {
			//run("Median...", "radius=1");
			if (BGSZ > 0){
				run("Subtract Background...", "rolling="+ BGSZ +" stack");
			}
			if (ROF == true){
				RemoveOutliersFilter("C3-Branch");
				rename("C3-Branch");
				print("Setting threshold for channel 2 to: "+ThresholdCh2);
				setThreshold(ThresholdCh2, 65535);
				run("Convert to Mask", "method=Default background=Dark black");
			} else {
				print("Setting the minimum intensity for channel 1 to: "+ThresholdCh2);
				run("Min...", "value="+ ThresholdCh2 +" stack");
				setMinAndMax(0, 4095);
				call("ij.ImagePlus.setDefault16bitRange", 12);
				run("8-bit");
				// otsu started to behave strangly on FIJI update 16/1/18 - replaced with mean
				run("Auto Local Threshold", "method=Phansalkar radius=10 parameter_1=0 parameter_2=0 white stack");
			}
		}	

		if (SecondPassC1 == true) {
			run("Duplicate...", "title=C2-BranchSP duplicate");
			selectWindow("C2-Branch");
		}


		
				
		run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels" +
			" integrated_density mean_gray_value std_dev_gray_value minimum_gray_value" +
			" centroid centre_of_mass show_masked_image_(redirection_requiered)" +
			" dots_size=5 font_size=10" +
			" store_results_within_a_table_named_after_the_image_(macro_friendly)" + 
			" redirect_to=C3-BranchRaw");
		
		getStatistics(area, mean, min, max, std, histogram);
		if (max < 128) threshold = max;
  		else threshold = 128;
		
		run("3D Objects Counter", "threshold="+threshold+" slice=63 min.=2 max.=33343088 objects statistics summary");
//		
		selectWindow("Statistics for C3-Branch redirect to C3-BranchRaw");
		saveAs("Results", analysisdir + "/Counts/"+ Branch + "_Ch2_Objects.csv");
		closewindow(Branch + "_Ch2_Objects.csv");
		//For validation
		selectWindow("Masked image for C3-Branch redirect to C3-BranchRaw");
		rename("C3-Mask");
		run("16-bit");
		selectWindow("Objects map of C3-Branch redirect to C3-BranchRaw");
		rename("ObjectMap");
		run("16-bit");
		// why C2 mask ok with out this step? convert c3 to 16bit check
		//run("16-bit");
		run("Merge Channels...", "c1=C3-Mask c2=C3-BranchRaw c3=ObjectMap create keep ignore");
		Stack.setChannel(3);
		run("3-3-2 RGB");
		
		rename("MergeCh2");
		run("Z Project...", "projection=[Max Intensity]");
		saveAs("Tiff", analysisdir + Branch + "_Ch2_Objects_Overlay.tif");
		close();
		selectWindow("MergeCh2");
		close();
		selectWindow("C3-Mask");
		close();
		closewindow("C3-BranchRaw");

		//Switched to objectmap as C3-Mask giving strange results
		selectWindow("ObjectMap");
		rename("C3-Mask");



		
		// crashed here?		to do with threshold set in 3d Objects counter.
		// settings to analyse distance (in pixels) to soma of channel 1 (C1)
		run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels mean_gray_value dots_size=5 font_size=10" +
			" store_results_within_a_table_named_after_the_image_(macro_friendly)" +
			" redirect_to=C1-Branch");
		selectWindow("C3-Branch");
		
		getStatistics(area, mean, min, max, std, histogram);
		if (max < 128) threshold = max;
  		else threshold = 128;
		
		run("3D Objects Counter", "threshold="+threshold+" slice=63 min.=2 max.=33343088 statistics summary");
		
		selectWindow("Statistics for C3-Branch redirect to C1-Branch");
		saveAs("Results", analysisdir +"/Counts/"+ Branch + "_Ch2_Distance.csv");
		closewindow(Branch + "_Ch2_Distance.csv");
		closewindow("C3-Branch");
				
		//For Colocalisation area and intensity analysis
		selectWindow("C2-Mask");
		setThreshold(1, 65535);
		setOption("BlackBackground", true);
		run("Convert to Mask", "method=Default background=Dark calculate black");
		selectWindow("C3-Mask");
		setThreshold(1, 65535);
		setOption("BlackBackground", true);
		run("Convert to Mask", "method=Default background=Dark calculate black");
		imageCalculator("AND create stack", "C2-Mask","C3-Mask");
		selectWindow("Result of C2-Mask");
		//measuring distance of colocalised regions!
		//run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels mean_gray_value dots_size=5 font_size=10 store_results_within_a_table_named_after_the_image_(macro_friendly) redirect_to=C1-Branch");
		//run("3D Objects Counter", "threshold=100 slice=63 min.=2 max.=33343088 statistics summary");
		
		// new version - corrects for when no objects in image and max threshold does not exist!
		getStatistics(area, mean, min, max, std, histogram);
		if (max < 128) threshold = max;
  		else threshold = 128;
		run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels mean_gray_value dots_size=5 font_size=10 store_results_within_a_table_named_after_the_image_(macro_friendly) redirect_to=C1-Branch");
		run("3D Objects Counter", "threshold="+threshold+" slice=63 min.=2 max.=33343088 statistics summary");
		
		selectWindow("Statistics for Result of C2-Mask redirect to C1-Branch");
		saveAs("Results", analysisdir +"/Coloc/"+ Branch + "_coloc_area_dist.csv");
		closewindow(Branch + "_coloc_area_dist.csv");
		//check this line
		filesopenlist = getList("image.titles");
		if (lengthOf(filesopenlist) > 3) {
			run("Merge Channels...", "c1=C2-Mask c2=C3-Mask c3=[Result of C2-Mask] create");
			selectWindow("Composite");
			run("Z Project...", "projection=[Max Intensity]");
			saveAs("Tiff", analysisdir + Branch + "_ColocOverlay.tif");
			closewindow("MAX_Composite");
			closewindow("Composite");
			closewindow("C1-Branch");	
		}
		
	} else {
		closewindow("C1-Branch");
		closewindow("C2-Mask");
	}
}


function closewindow(windowname) {
	if (isOpen(windowname)) { 
      		 selectWindow(windowname); 
       		run("Close"); 
  		} 
}
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

function subtitle(imagename){
	nl=lengthOf(imagename);
	nl2=nl-4;
	Sub_Title=substring(imagename,0,nl2);
	
	return Sub_Title;
}