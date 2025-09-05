/*
 Macro to convert nbit images into 8bit tiff images.
 */

setBatchMode(true);
arg_in=split(getArgument," ");

input = arg_in[0]; //path to image
output = arg_in[1]; //path to save images
index = arg_in[2]; //path to image series
name = arg_in[3] // name of the file to save
processFile(input, output, index, name);

function processFile(input, output, index, name) {
	
	path = input;
	print("Processing: " + path);
	
	//full_name = File.getName(path);
	//name = split(full_name, ".");
	//name = name[0];

	run("Bio-Formats Importer", "open=" + path +" color_mode=Default open_files display_metadata rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_list="+index);
    run("Remove Outliers...", "radius=1 threshold=0 which=Bright");
	
	run("8-bit");
	Stack.setChannel(1);
	run("Grays");
	Stack.setChannel(2);
	run("Red");
	Stack.setChannel(3);
	run("Yellow");
	
	saveAs("Tiff", output + File.separator + name);
	close();
}
