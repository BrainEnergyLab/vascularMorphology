import ij.IJ; //Import this so that we can use imageJ to open the image

//This is basically declaring a function
public ArrayList listFilesAndFilesSubDirectories(String directoryName, String substring) {

	//Here we declare our folder as a file variable for our directory and create a file array of the files within folder
	File folder = new File(directoryName);
	File [] listOfFiles = folder.listFiles();
	ArrayList fileLocations = new ArrayList();
	

	//Loop through the files in the file list
	for (File file:listOfFiles) {
		//If the file we're checking is a file and not a directory
		if (file.isFile()) {
			//And if it contains the substring we're interested in within its full path
			if(file.getAbsolutePath().contains(substring)) {
				//We print the full name and open the image
				System.out.println(file.getAbsolutePath());
				fileLocations.add(file.getAbsolutePath());
				//imp=IJ.openImage(file.getAbsolutePath);
				//imp.show();
			}
		//If the file we're checking is a directory, then we run the whole thing on that directory
		} else if (file.isDirectory()) {
			listFilesAndFilesSubDirectories(file.getAbsolutePath(), substring);           
		}
	}

	return fileLocations;
}

//Example inputs for the function (called a method in java?)
String dir = "E:/Dropbox (Brain Energy Lab)/Everything/2P data/Devin/";
String toFind = "WideFOV";

System.out.println("Start");
listFilesAndFilesSubDirectories(dir, toFind);
System.out.println("End");
