# epskw
This code analyzes molecular dynamics trajectory files. It currently can read .xyz or .xtc files and works with water, methanol, and acetonitrile. It can easily be extended to accomidate any simulation of a pure polar liquid. 

The code can compute the following quantities:

* static longitudinal and transverse nonlocal susceptibiliies, chi_L(k,0), chi_T (k,0)
* longitudinal and transverse polarization correlation functions phiL(k, t), phiT(k, t) 
* nonlocal susceptibiliies, chi(k,0), by feeding the output from this code into the kw.m code
* distance decompositions chi_L(k,0,R), chi_T (k,0,R), phiL(k, t, R), phiT(k, t, R), for the smallest k-value accessible in the system.  
* static structure factors S_L(k,0) and S_T(k,0) and (experimental feature, may not be normalized correctly)

## Compilation instructions
To compile you need the xtc library libxdrf.a. A precompiled version is included in this repo,  but chances are it will not work on your system and you will need to compile your own version. To obtain this library, download the source code from [ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.tar.gz](ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.tar.gz), edit the makefile as necessary and make. Copy libxdrf.a into the /epskw folder. Run "make" to compile

## Setting the input options 
The included file epskw.inp gives an example of an input file, which is fixed format. To use a small set of k values (10-20), set *use small k set ?* to *.t.*. To use a dense set of k-points, select a maximum number of points to use using the diagonals and edges of the box. 

The following forcefield models are supported: 
*TIP4P, TIP4P2005, TIP4P2005f, SPCE, TIP3P, 
*TTM3F (requires extra file containing charges of all atoms in each timestep in the order OHHOHHOHH.., no spaces between timesteps)
*methanol (GAFF) 
*methanolH1 (H1+3 parameters)
*acetonitrile (GAFF)
*generic

To add your own forcefield model, edit main_stuff.f90 

## Running the program 
to run, use: 
`./epskw.x < epskw.inp`

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

