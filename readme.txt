Instructions for compiling and running 2pi.model.f code
Details regarding the model can be found in our paper in the Journal of Computational Neuroscience
https://doi.org/10.1007/s10827-022-00815-x

The following files are needed:
2pi.in = input file for changing some model parameters
   	 this can be modified without recompiling the code
2pi.model.f = main program
const.f = module file that contains constants
ode.par = parameter file used by main program for ode solver

Follow these steps to compile the code
1) compile const.f module
 	gfortran -fdefault-real-8 -c const.f 

2) compile 2pi.model.f code
	mpif90 -fdefault-real-8 -c 2pi.model.f

3) link module and code
	mpif90 const.o 2pi.model.o -o 2pi.exe

4) code can be executed using MPI. Max processors is from 2pi.in (Jeend - Jestart)/Jestep + 1
	mpirun -np xx 2pi.exe
  (xx = number of processors)

The code produces the following output files:
1) frequency.out.int.Jex.x.Jiy.y.gabaz.z
	contains the firing frequency (in Hz) of the interneuron as a function of 
	pyramidal cell injected current (Je), interneuron injected current (Ji),
 	gaba maximum conductance (ggaba), and time
	Files are labeled as x.x = Je value, y.y = Ji value, z.z=ggaba value
2) frequency.out.pyr.Jex.x.Jiy.y.gabaz.z
	same as (1), but for the pyramidal cell
3) if writeyn=1, membrane.pot.out.Jex.x.Jiy.y.gabaz.z
	contains the pyramidal cell and interneuron membrane potential (mV)
  	as a function of time - these are large files
