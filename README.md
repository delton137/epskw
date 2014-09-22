
Code to calculate the nonlocal dielectric function eps(k,omega) and static structure factor for water. Reads .xtc or .xyz files 

COMPILATION INSTRUCTIONS: To compile you need the xtc library libxdrf.a. The library file is included, but chances are it will not work on your system and you will need to compile your own version. 

To obtain this library, download the source code from ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.tar.gz , edit the makefile as necessary and make. Copy libxdrf.a into the /epskw folder. Run "make" to compile


Output files are outputted in a form readible by xmgrace. To view the figures use: 
xmgrace -nxy "output_filename"

2014 Dan Elton 

