! this module contains the constants for the 2PI Model code
!  '2pi.model.f'
! It must be compiled first, before the main code is compiled.
	
! compile this module with
! gfortran -fdefault-real-8 -c const.f

	module const
	implicit none
	save

! declare variables
	real gNa,gp,gAHP,gCa,gk,C,phi,epsi,tauCa
	real RTF,temp,gasconst
	real beta,gam,tau
	real F,S,vol,tauinv
	real gkL,gClL,gNaL
	real epsK,Kbath,pi,rhopump,Nasat,Ksat
	real rhoKcc,rhoNKcc,taugaba
	real gglut1,gglut2,tauglut,ECa

	real C_int,phi_int,epsi_int,gCa_int,tauCa_int
	real gp_int,gNa_int,gAHP_int,gCl_int,gK_int
	real EL_int,gL_int,Eglut
	real gkL_int,gNaL_int

! set constants for pyramidal cell and synapses
	parameter(pi=acos(-1.))
	parameter(gk = 80.,gNa = 100.,gp = 1.,gAHP = 1.5)	!mS/cm^2 potassium, sodium, persistent sodium, calcium dependent potassium max conductances	
	parameter(C = 1.,phi = 1.,epsi = 0.002,tauCa = 80.)	!membrane capacitance microF/cm^2, Conversion b/w Ca current density and ion concentration change (mmol/C cm), time const Ca (ms)
	parameter(F = 96485.3)		!faraday constant A.s/mol
	parameter(temp = 37.+273.)	!temperature Kelvin
	parameter(gasconst = 8.3145)	!J/mol.K
	
	parameter(RTF = 1000.*gasconst*temp/F) 		!R*T/F in mV
	parameter(beta = 4.)				!ratio of intracellular to extracellular space
	parameter(vol = 1.4368e-9)			!cm^3	volume of cell
	parameter(S = 4.*pi*(3.*vol/(4.*pi))**(2./3.))	!surface of cell
	parameter(gam = S/(F*vol))			!M/C.cm conversion between current density and rate of ion concentration change
	parameter(tau = 1000.)				!conversion ms to s
	parameter(tauinv = 1./tau)			
	parameter(gkL = 0.05)				!mS/cm^2 max potassium leak conductance
	parameter(gClL = 0.015)				!mC/cm^2 max chloride leak conductance
	parameter(gNaL = 0.0015)		 	!mS/cm^2 max sodium leak conductance
	parameter(rhopump = 0.25)			!mM/s Na/K pump strength
	parameter(Nasat = 22.)			 	!mM pyramidal cell intracellular sodium concentration for pump at half capacity
	parameter(Ksat = 3.5)				!mM extracellular potassium concentration for pump at half capacity
	parameter(rhoKcc = 0.3)				!mM/s pyramidal cell max KCC2 cotransporter strength
	parameter(rhoNKcc = 0.1)			!mM/s pyramidal cell max NKCC1 cotransporter strength
	parameter(epsK = 0.4)				!1/s rate factor potassium loss to environment
	parameter(Kbath = 3.5)			 	!mM environment potassium concentration
	parameter(taugaba = 9.)				!ms decay constant for GABA synapse
	parameter(gglut1 = 0.1,gglut2=gglut1)		!mS/cm^2 max glutamate conductance
	parameter(tauglut = 3.)				!ms decay constant for excitatory synapse on interneuron
	parameter(gCa = 1.)				!mS/cm^2 max Calcium conductance
	parameter(ECa = 120.)				!mV calcium equilibrium potential
	parameter(Eglut = 0.)				!mV glutamate equilibrium potential

! set constants for interneuron 
	parameter(C_int = 1.)				!microF/cm^2 membrane capacitance
	parameter(phi_int = 5.)				!temperature factor for gating variables
	parameter(gNa_int = 35.)			!mS/cm^2 max sodium conductance
	parameter(gK_int = 9.)				!mS/cm^2 max potassium conductance
	parameter(EL_int = -65.)			!mV leak equilibrium potential 
	parameter(gL_int = 0.1)				!mS/cm^2 max leak conductance
	parameter(gkL_int = 0.08276)			!mS/cm^2 max potassium leak conductance
	parameter(gNaL_int = 0.0172)		 	!mS/cm^2 max sodium leak conductance

	
	end module const
