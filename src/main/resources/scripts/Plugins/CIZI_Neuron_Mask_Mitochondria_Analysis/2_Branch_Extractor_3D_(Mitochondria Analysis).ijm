// 3d Region Extractor
 
// Author: 	Luke Hammond
// Cellular Imaging | Zuckerman Institute, Columbia University
// Version: 0.6
// Date:	27th October 2017
//	
// This macro was written to allow ROIs saved in a file to be used to extract
// selected regions from corresponding 3D neuron images - masked by neuron masking tool. 
// It also copies in a geodesic distance map image so objects can be measured relative to their
// distance to the soma.

// Input directory should contain 3 subfolders: 
//		/Neurons (contains 3D neuron images)
//		/ROIs (contains saved ROI lists)
//		/Distance (contains distance map images - 16bit)

//*Future update: can be modified to create a distance map from ROI 0. then use this distance map to create a
//distance channel for each extracted branch. distances are in pixels for distance map images.

// Update v3: multiple branch ROIs are merged together to create one extracted branch image.
// Update v4 12/19/2017: Added 3D distance map, which can be used to measure volume at distance from soma
//						The volume of the neuron in 10um intervals are printed to a volume results table 
// Update v5 1/31/2018: Ignores subfolders when processing images using new function " ImageFilesOnlyArray" - call by " arrayname = ImageFilesOnlyArray(arrayname)
// Update v6 3/2/2018: Updated folder structure


// Initialization
requires("1.41h");
starttime = getTime();
run("Options...", "iterations=3 count=1 black do=Nothing");
run("Set Measurements...", "fit redirect=None decimal=3");
run("Colors...", "foreground=white background=black selection=yellow");
run("Clear Results");
cleanupROI();

// menu
Dialog.create("3D Neuron Region Extraction Tool");

Dialog.addNumber("Number of masked channels", 1);
//Dialog.addCheckbox("Extract Axons?",true);
//Dialog.addCheckbox("Extract Primary Branches?",true);
//Dialog.addCheckbox("Extract Terminal Branches?",true);

Dialog.show();
NumOfChan = Dialog.getNumber();
//AxOn = Dialog.getCheckbox();
//PBOn = Dialog.getCheckbox();
//TBOn = Dialog.getCheckbox();

// Preparation
input = getDirectory("Input directory");
print("\\Clear");
print("3D branch extractor running:");
File.mkdir(input + "Extracted_Regions");
File.mkdir(input + "Distance_Maps");
File.mkdir(input + "Analyzed");
setBatchMode(true);

files = getFileList(input + "Processed");

files = Array.sort( files );
	
//iterate over all files
for(i=0; i<files.length; i++) {	
	//Process Axons
	//update has -1 as folder also contains log from previous run - ?
	print("\\Update1:Processing neuron " + (i+1) + " of " + (files.length) + ".");
	branchesrois = getFileList(input + "MIPs/Branches/ROIs");
	branchesrois = Array.sort( branchesrois );

	if (lengthOf(branchesrois) == 0 ) {
		exit("No axon ROIs to process. Please create ROIs around cell bodies and run again.");
	} else {
		
		
		File.mkdir(input + "Extracted_Regions/Branches");
		File.mkdir(input + "Analyzed/Branches");
		File.mkdir(input + "Analyzed/Branches/Volume");
		File.mkdir(input + "Distance_Maps/Branches");

		
		//Create Distance Map
		Branches =  getFileList(input + "MIPs/Branches");
		Branches = ImageFilesOnlyArray(Branches);
		Branches = Array.sort( Branches );
		MIP = Branches[i];
		image = files[i];
		ROI = branchesrois[i];
		cleanupROI();
		//Assumes C1 is neuron mask and only imports this.
		Output_title = subtitle(image);
		run("Bio-Formats Importer", "open=[" + input +"MIPs/Branches/"+ MIP + "] autoscale color_mode=Default specify_range view=Hyperstack stack_order=XYCZT c_begin=1 c_end=1 c_step=1");
		rename("Neuron");
		run("Make Binary");
		run("Duplicate...", "title=Marker");
		roiManager("Open", input +"MIPs/Branches/ROIs/"+ ROI + "");
		roiManager("Select", 0);
		run("Clear Outside");
		run("Select None");
		run("Geodesic Distance Map", "marker=Marker mask=Neuron distances=[Chessknight (5,7,11)] output=[16 bits] normalize");
		rename("Distance");
		//Tif_Title = tiftitlespace(MIP);
		saveAs("Tiff", input + "Distance_Maps/Branches/"+ MIP);
		rename("Distance");
		selectWindow("Marker");
		close();
		selectWindow("Neuron");
		close();
		
		// Open Neuron Image "Processed" and create the distance channel
		run("Bio-Formats Importer", "open=["+ input +"Processed/" + image + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
		rename("Processed");
		getDimensions(width, height, ChNum, slices, frames);
		run("Split Channels");
		selectWindow("C1-Processed");
		rename("DistMask");
		run("Duplicate...", "title=DistanceStack duplicate");
		selectWindow("DistMask");
		setMinAndMax(0, 255);
		run("Apply LUT", "stack");
		run("Invert", "stack");
		
		
		selectWindow("DistanceStack");
		for (z=0; z<slices; z++){
			setSlice(z+1);
			selectWindow("Distance");
			run("Select All");
			run("Copy");
			selectWindow("DistanceStack");
			run("Paste");
		}
		selectWindow("Distance");
		close();
		imageCalculator("Transparent-zero stack", "DistanceStack","DistMask");
		//added in close diskmask
		//selectWindow("DistMask");
		//close();
				
		if (NumOfChan == 1){
			run("Merge Channels...", "c1=DistanceStack c2=C2-Processed create");
		} else {
			run("Merge Channels...", "c1=DistanceStack c2=C2-Processed c3=C3-Processed create");
			Stack.setChannel(3);
			run("Red");
		}
		selectWindow("Composite");
		rename("Processed");

		//Extract and create regions
		//CountROImain=roiManager("count");
		//ROIarrayMain=newArray(CountROImain);
		
		//REMOVED AS AXON SHOULD ONLY BE ONE ROI
		//for(r=1; r<ROIarrayMain.length; r++) {
		//	print("\\Update2: extracting axon " + (i+1) +"/" + ROIarrayMain.length +".");
		//	selectWindow("Processed");
		//	roiManager("Select", r);
		//	run("Duplicate...", "title=["+ image +"_region_"+ r +"] duplicate");
		//	run("Clear Outside", "stack");
		//	//run("Bio-Formats Exporter", "save=["+ input +"Neuron_Regions/"+ neuron +"_region_"+ i +".ome.tif] export compression=Uncompressed");
		//	saveAs("Tiff", input +"Neuron_Regions/Axons/"+ image +"_region_"+ r +".tif");
		//}

		selectWindow("Processed");
		//roiManager("Select", 1);
		//run("Duplicate...", "title=["+ image +"_Branch] duplicate");
		//run("Clear Outside", "stack");
		saveAs("Tiff", input +"Extracted_Regions/Branches/"+ Output_title);
		roiManager("Deselect");
		roiManager("Delete");		
		run("Select None");
		ROIarrayMain=newArray(0);
		CountROImain=0;
		rename("Region");
		run("Split Channels");
		closewindow("C2-Region");
		if (NumOfChan == 2){
			closewindow("C3-Region");	
		}
		selectWindow("C1-Region");
		

		tablename = "[" + image + "Volume table]";
		run("New... ", "name="+tablename+" type=Table");
		f = tablename;
		print(f, "\\Headings:Distance from soma (um)\tTotalArea (um2)\tVolume (um3)");
	
		for(m=1; m<1001; m++) {	
			n = round(((m+99)/10));
			
			sum = sum3DArea("C1-Region", m, (m+99));
			m = m+100;
			print(f,  n+ "\t" + sum + "\t" + (sum*0.2) );
		}
		
		closewindow("C1-Region");
		selectWindow(image + "Volume table");
		//selectWindow("Summary of ");
		saveAs("text", input +"Analyzed/Branches/Volume/"+ Output_title + "_Volume.csv");
		closewindow(image + "Volume table");
	}

/*
	//Process Dendrites (instead of processing each ROI we need to combine to one image)
	dendritesrois = getFileList(input + "MIPs/Primary_Branches/ROIs");
	Array.sort( dendritesrois );
	if (lengthOf(dendritesrois) == 0 || PBOn == false){
		print("\\Update3: No Primary Branch ROIs to process.");
	} else {
		print("\\Update3: Processing Primary Branches.");
		File.mkdir(input + "Extracted_Regions/Primary_Branches");
		File.mkdir(input + "Analyzed/Primary_Branches");
		File.mkdir(input + "Analyzed/Primary_Branches/Volume");
		File.mkdir(input + "Distance_Maps/Primary_Branches");

		//Create Distance Map
		dendrites =  getFileList(input + "MIPs/Primary_Branches");
		dendrites = ImageFilesOnlyArray(dendrites);
		Array.sort( dendrites );
		
		MIP = dendrites[i];
		image = files[i];
		ROI = dendritesrois[i];
		cleanupROI();
		//Assumes C1 is neuron mask and only imports this.
		run("Bio-Formats Importer", "open=[" + input +"MIPs/Primary_Branches/"+ MIP + "] autoscale color_mode=Default specify_range view=Hyperstack stack_order=XYCZT c_begin=1 c_end=1 c_step=1");
		rename("Neuron");
		run("Make Binary");
		run("Duplicate...", "title=Marker");
		roiManager("Open", input +"MIPs/Primary_Branches/ROIs/"+ ROI + "");
		roiManager("Select", 0);
		run("Clear Outside");
		run("Select None");
		run("Geodesic Distance Map", "marker=Marker mask=Neuron distances=[Chessknight (5,7,11)] output=[16 bits] normalize");
		rename("Distance");
		//Tif_Title = tiftitlespace(MIP);
		saveAs("Tiff", input + "Distance_Maps/Primary_Branches/"+ MIP);
		rename("Distance");
		selectWindow("Marker");
		close();
		selectWindow("Neuron");
		close();
		
		// Open Neuron Image "Processed" and create the distance channel
		filesopenlist = getList("image.titles");
		if (lengthOf(filesopenlist) > 1) {
			print("\\Update4: Dataset already imported, pasting new distance map and proceeding to remove ROIs.");
			selectWindow("Processed");
			run("Select None");
			run("Split Channels");
			selectWindow("C1-Processed");
			rename("DistanceStack");
			
		}
		else {
			run("Bio-Formats Importer", "open=["+ input +"Processed/" + image + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
			rename("Processed");
			getDimensions(width, height, ChNum, slices, frames);
			run("Split Channels");
			selectWindow("C1-Processed");
			rename("DistMask");
			run("Duplicate...", "title=DistanceStack duplicate");
			selectWindow("DistMask");
			setMinAndMax(0, 255);
			run("Apply LUT", "stack");
			run("Invert", "stack");
			
		}
			
		
		selectWindow("DistanceStack");
		for (z=0; z<slices; z++){
			setSlice(z+1);
			selectWindow("Distance");
			run("Select All");
			run("Copy");
			selectWindow("DistanceStack");
			run("Paste");
		}
		selectWindow("Distance");
		close();
		imageCalculator("Transparent-zero stack", "DistanceStack","DistMask");
		//added in close diskmask
		//selectWindow("DistMask");
		//close();
		

		if (NumOfChan == 1){
			run("Merge Channels...", "c1=DistanceStack c2=C2-Processed create");
		} else {
			run("Merge Channels...", "c1=DistanceStack c2=C2-Processed c3=C3-Processed create");
			Stack.setChannel(3);
			run("Red");
		}
		selectWindow("Composite");
		rename("Processed");

		//Extract and create regions
		CountROImain=roiManager("count");
		ROIarrayMain=newArray(CountROImain);
		if (CountROImain==2) {
			print("\\Update2: extracting primary dendrite.");
			selectWindow("Processed");
			roiManager("Select", 1);
			run("Duplicate...", "title=["+ image +"_region] duplicate");
			run("Clear Outside", "stack");
			//run("Bio-Formats Exporter", "save=["+ input +"Neuron_Regions/"+ neuron +"_region_"+ i +".ome.tif] export compression=Uncompressed");
			saveAs("Tiff", input +"Extracted_Regions/Primary_Branches/"+ image +"_region.tif");
		} else if (CountROImain>2) {
			roiManager("Select", 0);
			roiManager("Delete");
			ROIarraysub=newArray((CountROImain-1)); 
			for(m=0; m<(CountROImain-1);m++) { 
        		ROIarraysub[m] = m; 
				} 
			roiManager("Select", ROIarraysub);
			roiManager("Combine");
			run("Duplicate...", "title=["+ image +"_region] duplicate");
			run("Clear Outside", "stack");
			//run("Bio-Formats Exporter", "save=["+ input +"Neuron_Regions/"+ neuron +"_region_"+ i +".ome.tif] export compression=Uncompressed");
			saveAs("Tiff", input +"Extracted_Regions/Primary_Branches/"+ image +"_region.tif");
						
		}
		rename("Region");
		run("Split Channels");
		closewindow("C2-Region");
		if (NumOfChan == 2){
			closewindow("C3-Region");	
		}
		selectWindow("C1-Region");
		

		tablename = "[" + image + "Volume table]";
		run("New... ", "name="+tablename+" type=Table");
		f = tablename;
		print(f, "\\Headings:Distance from soma (um)\tTotalArea (um2)\tVolume (um3)");
	
		for(m=1; m<1001; m++) {	
			n = round(((m+99)/10));
			
			sum = sum3DArea("C1-Region", m, (m+99));
			m = m+99;
			print(f,  n+ "\t" + sum + "\t" + (sum*0.2) );
		}
		n = "remaining";
		sum = sum3DArea("C1-Region", 1001, 32000);
		print(f,  n+ "\t" + sum + "\t" + (sum*0.2) );
		
		
		
		closewindow("C1-Region");
		
		
		//saveAs(f, input +"Neuron_Regions/Primary_Branches/Analyzed/"+ image + "_Primary_Branch_Volume.csv");
		//selectWindow("Summary of ");	
		selectWindow(image + "Volume table");
		saveAs("text", input +"Analyzed/Primary_Branches/Volume/"+ image + "_Primary_Branch_Volume.csv");
		closewindow(image + "Volume table");
		
		roiManager("Deselect");
		roiManager("Delete");
		ROIarrayMain=newArray(0);
		ROIarraysub=newArray(0);
		CountROImain=0;
		
	}

	//Process Terminal Branches
	terminalsrois = getFileList(input + "MIPs/Terminal_Branches/ROIs");
	Array.sort( terminalsrois );
	
	if (lengthOf(terminalsrois) == 0 || TBOn == false) {
		print("\\Update4: No Terminal Branch ROIs to process.");
	} else {
		print("\\Update4: Processing Terminal Branches.");
		File.mkdir(input + "Extracted_Regions/Terminal_Branches");
		File.mkdir(input + "Analyzed/Terminal_Branches");
		File.mkdir(input + "Analyzed/Terminal_Branches/Volume");
		File.mkdir(input + "Distance_Maps/Terminal_Branches");

		//Create Distance Map
		terminals =  getFileList(input + "MIPs/Terminal_Branches");
		terminals = ImageFilesOnlyArray(terminals);
		Array.sort( terminals );
		
		MIP = terminals[i];
		image = files[i];
		ROI = terminalsrois[i];
		cleanupROI();
		//Assumes C1 is neuron mask and only imports this.
		run("Bio-Formats Importer", "open=[" + input +"MIPs/Terminal_Branches/"+ MIP + "] autoscale color_mode=Default specify_range view=Hyperstack stack_order=XYCZT c_begin=1 c_end=1 c_step=1");
		rename("Neuron");
		run("Make Binary");
		run("Duplicate...", "title=Marker");
		roiManager("Open", input +"MIPs/Terminal_Branches/ROIs/"+ ROI + "");
		roiManager("Select", 0);
		run("Clear Outside");
		run("Select None");
		run("Geodesic Distance Map", "marker=Marker mask=Neuron distances=[Chessknight (5,7,11)] output=[16 bits] normalize");
		rename("Distance");
		//Tif_Title = tiftitlespace(MIP);
		saveAs("Tiff", input + "Distance_Maps/Terminal_Branches/"+ MIP);
		rename("Distance");
		selectWindow("Marker");
		close();
		selectWindow("Neuron");
		close();
		
		// Open Neuron Image "Processed" and create the distance channel
		filesopenlist = getList("image.titles");
		if (lengthOf(filesopenlist) > 1) {
			print("\\Update5: Dataset already imported, pasting new distance map and proceeding to remove ROIs.");
			selectWindow("Processed");
			run("Select None");
			run("Split Channels");
			selectWindow("C1-Processed");
			rename("DistanceStack");
			
		}
		else {
			run("Bio-Formats Importer", "open=["+ input +"Processed/" + image + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
			rename("Processed");
			getDimensions(width, height, ChNum, slices, frames);
			run("Split Channels");
			selectWindow("C1-Processed");
			rename("DistMask");
			run("Duplicate...", "title=DistanceStack duplicate");
			selectWindow("DistMask");
			setMinAndMax(0, 255);
			run("Apply LUT", "stack");
			run("Invert", "stack");
			
		}
			
		
		selectWindow("DistanceStack");
		for (z=0; z<slices; z++){
			setSlice(z+1);
			selectWindow("Distance");
			run("Select All");
			run("Copy");
			selectWindow("DistanceStack");
			run("Paste");
		}
		selectWindow("Distance");
		close();
		imageCalculator("Transparent-zero stack", "DistanceStack","DistMask");
		//added in close diskmask
		//selectWindow("DistMask");
		//close();

		
		if (NumOfChan == 1){
			run("Merge Channels...", "c1=DistanceStack c2=C2-Processed create");
		} else {
			run("Merge Channels...", "c1=DistanceStack c2=C2-Processed c3=C3-Processed create");
			Stack.setChannel(3);
			run("Red");
		}
		selectWindow("Composite");
						
		rename("Processed");

		//Extract and create regions
		CountROImain=roiManager("count");
		ROIarrayMain=newArray(CountROImain);
		if (CountROImain==2) {
			print("\\Update2: extracting terminals.");
			selectWindow("Processed");
			roiManager("Select", 1);
			run("Crop");
			//run("Duplicate...", "title=["+ image +"_region] duplicate");
			run("Clear Outside", "stack");
			run("Select None");
			//run("Bio-Formats Exporter", "save=["+ input +"Neuron_Regions/"+ neuron +"_region_"+ i +".ome.tif] export compression=Uncompressed");
			saveAs("Tiff", input +"Extracted_Regions/Terminal_Branches/"+ image +"_region.tif");
		} else if (CountROImain>2) {
			roiManager("Select", 0);
			roiManager("Delete");
			ROIarraysub=newArray((CountROImain-1)); 
			for(m=0; m<(CountROImain-1);m++) { 
        		ROIarraysub[m] = m; 
				} 
			roiManager("Select", ROIarraysub);
			roiManager("Combine");
			run("Crop");
			//run("Duplicate...", "title=["+ image +"_region] duplicate");
			run("Clear Outside", "stack");
			//run("Bio-Formats Exporter", "save=["+ input +"Neuron_Regions/"+ neuron +"_region_"+ i +".ome.tif] export compression=Uncompressed");
			saveAs("Tiff", input +"Extracted_Regions/Terminal_Branches/"+ image +"_region.tif");			
		}
		
		rename("Region");
		run("Split Channels");
		closewindow("C2-Region");
		if (NumOfChan == 2){
			closewindow("C3-Region");	
		}
		selectWindow("C1-Region");
		tablename = "[" + image + "Volume table]";
		run("New... ", "name="+tablename+" type=Table");
		print(tablename, "\\Headings:Distance from soma (um)\tTotalArea (um2)\tVolume (um3)");
		n = "all";
		sum = sum3DArea("C1-Region", 1, 32000);
		
		print(tablename,  n+ "\t" + sum + "\t" + (sum*0.2));

		//Need to add total number of terminal branches?
			
		closewindow("C1-Region");
		//selectWindow("Summary of ");
		selectWindow(image + "Volume table");
		saveAs("text", input +"Analyzed/Terminal_Branches/Volume/"+ image + "_Terminal_Volume.csv");	
		//saveAs("Results", directory +"/Analyzed/"+ Branch + "coloc_area_dist.csv");		
		closewindow(image + "Volume table");
		
		roiManager("Deselect");
		roiManager("Delete");
		ROIarrayMain=newArray(0);
		ROIarraysub=newArray(0);
		CountROImain=0;
		
	}
*/	

closewindow("Processed");
closewindow("DistMask");
	
}

endtime = getTime();
dif = (endtime-starttime)/1000;
print("\\Update6:Finished!");
print("\\Update7: Processing time =", (dif/60), "minutes total. " + (dif/i) + "seconds per neuron.");

selectWindow("Log");
saveAs("txt", input+"/3D_Region_Extractor_Log.txt");


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

function sum3DArea(imagename, min, max) {
	run("Set Measurements...", "area limit redirect=None decimal=3");
	selectWindow(imagename);
	setThreshold(min, max);
	run("Analyze Particles...", "clear summarize stack");
	Sum = 0; 
	lines1 = split(getInfo(), "\n"); 
	headings = split(lines1[0], "\t"); 
	for (j=1; j<lines1.length; j++) { 
	        values = split(lines1[j], "\t"); 
	        a = values[2]; 
	        sum = sum + a; 
	}
	closewindow("Summary of "+ imagename);
	return sum;
}


function sumcolumn(col) {
	Sum = 0; 
	//selectWindow("Results"); 
	lines1 = split(getInfo(), "\n"); 
	headings = split(lines1[0], "\t"); 
	for (j=1; j<lines1.length; j++) { 
	        values = split(lines1[j], "\t"); 
	        a = values[col]; 
	        sum = sum + a; 
	}
	return sum;
}

function closewindow(windowname) {
	if (isOpen(windowname)) { 
      		 selectWindow(windowname); 
       		run("Close"); 
  		} 
}

function ImageFilesOnlyArray (arr) {
	//pass array from getFileList through this e.g. NEWARRAY = ImageFilesOnlyArray(NEWARRAY);
	setOption("ExpandableArrays", true);
	f=0;
	files = newArray;
	for (i = 0; i < arr.length; i++) {
		if(endsWith(arr[i], ".tif") || endsWith(arr[i], ".nd2") ) {   //if it's a tiff image add it to the new array
			files[f] = arr[i];
			f = f+1;
		}
	}
	arr = files;
	Array.sort(arr);
	return arr;
}

function subtitle(imagename){
	nl=lengthOf(imagename);
	nl2=nl-4;
	Sub_Title=substring(imagename,0,nl2);
	
	return Sub_Title;
}
