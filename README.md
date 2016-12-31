# epskw
A code to calculate the nonlocal dielectric function eps(k,omega) and static structure factor for trajectory files containing pure water. Reads .xtc or .xyz files 

## COMPILATION INSTRUCTIONS
To compile you need the xtc library libxdrf.a. A precompiled version is included in this repo,  but chances are it will not work on your system and you will need to compile your own version. To obtain this library, download the source code from [ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.tar.gz](ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.tar.gz), edit the makefile as necessary and make. Copy libxdrf.a into the /epskw folder. Run "make" to compile

## Running the program 
Edit the fixed format input file epskw.inp to specify your coordinate file, output title, and options for calculation.  
to run, use: 
`.epskw.x < epskw.inp`

## Plotting the results  
Output files are outputted in a form readible by xmgrace. To view the figures use:  

`xmgrace -nxy *output_filename*`


## License 
Copyright 2014-2016 Daniel C. Elton

This software is licensed under The MIT License (MIT)
Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software
without restriction, including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software, and to permit 
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

