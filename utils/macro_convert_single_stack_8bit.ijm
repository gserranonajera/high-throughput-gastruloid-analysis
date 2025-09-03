/*
 Macro to convert nbit images into 8bit tiff images.
 */

setBatchMode(true);
arg_in=split(getArgument," ");

input = arg_in[0]; //path to image
output = arg_in[1]; //path to save images

processFile(input, output);

function processFile(input, output) {
	
	path = input;
	print("Processing: " + path);
	
	full_name = File.getName(path);
	name = split(full_name, ".");
	name = name[0];

    run("Bio-Formats Importer", "open=" + path +" autoscale color_mode=Default concatenate_series open_all_series display_metadata rois_import=[ROI manager] view=Hyperstack stack_order=Default");	
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
