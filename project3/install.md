 # Installation and Execution Guide 
 
 ## Prerequisites
 Before proceeding, ensure the following tools are installed on your system:
 • GCC (GNU Compiler Collection)
 • Make utility
 
 ## Steps to Install and Run the Code
 
 1. Download the Repository
 Download or clone the Project3 repository to your local system:
 git clone <repository url >
 Navigate into the project3 directory: cd project3
 
 2. Understanding the Directory Structure
 The repository has the following structure:
 • Makefile: Contains commands for make, make run, and make clean.
 • src/: Contains the source files (dynamics.c, utils.c, error.c).
 • data/: Contains input files and will hold the generated output files.
 
 3. Build the Project
 Compile the source code using the Makefile:
 make
 This will generate the necessary executable files.

 4. Run the Code
 Execute the program using the following command:
 make run
 The program will process the input files located in the data/ directory and
 generate the output files in the same folder.
 
5. Clean the Project
 To remove any compiled files or intermediate build files, use the following command:
 make clean
 
 # Output
 The output files will be generated in the data/ directory after successful execution.
 For any issues or troubleshooting, refer to the project documentation or contact the repository maintainer. 
