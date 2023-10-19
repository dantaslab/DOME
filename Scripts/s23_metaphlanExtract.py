#!/usr/bin/env python3
""" Extracts metaphlan data at every taxonmic level

Walks through every .txt file in the root of a given directory, and for each 
text file, creates files by taxonomic level in a subfolder named "Extracted",
writing to the newly created file only the lines from the original file that
match the correct number of delimiters.

Author: Bejan Mahmud, bmahmud@wustl.edu (via Alaric D'Souza)

Example usage: python3 metaphlanExtract.py d05_metaphlan/

Args:
    path: Path of directory with unprocessed .txt files requiring extraction

"""

import glob, os, re, sys

def main(path):
    # Go to the specified directory
    # Note that default directory is location of script file
    if(os.path.exists(path)):
        os.chdir(path)
    else:
        sys.stderr.write("Error: Given path not valid")
        sys.exit(1)

    # Create a new folder to contain the output files
    extractedPath = os.path.join(os.getcwd(), "Extracted")
    if not os.path.exists(extractedPath):
        #make extracted path
        os.makedirs(extractedPath)

        #initialize dictionary of taxonmic levels
        taxdict={
            "phyla":2,
            "class":3,
            "order":4,
            "family":5,
            "genus":6,
            "species":7}

        #make taxonomic folders
        for taxa in taxdict.keys():
            os.makedirs(os.path.join(extractedPath,taxa))

    # Iterate through all the text files of the given directory    
    for file in glob.glob('*.txt'):
        #read file to array
        with open(file,'r') as f:
            f.readline()
            myfile=[line.strip().split("\t") for line in f.readlines()]

            newfile = []
            newline = ""

            #Skip the lines taht start with "#"
            #Some lines also report "additional species", starting with "k__***" after the abundance value. Such additional species items will be removed so that all of the lines have abundance values as the last item
            for item in myfile:
                newline = item
                if newline[0].startswith("#"):
                    continue
                elif newline[-1].startswith("k"):
                    newline = newline[:-1]
                newline[0] = newline[0].split("|")
                newfile.append(newline)

        #write out taxa specific files
        for taxa,lenline in taxdict.items():
            with open(os.path.join(extractedPath,taxa,file.replace(".txt","_"+taxa+".txt")),"w+") as outfile:
                #write header line
                outfile.write(taxa+"\tabundance\n")
                #write relative abundance data
                outfile.write("\n".join(["\t".join([line[0][-1],line[-1]]) for line in newfile if len(line[0])==lenline]))

if __name__ == "__main__":
    print(sys.argv)
    if len(sys.argv) == 1:
        main(os.path.dirname(os.path.realpath(__file__)))
    elif len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        sys.stderr.write("Error: Invalid number of parameters\n")
        sys.exit(1)