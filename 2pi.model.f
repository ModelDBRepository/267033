! This program is based on the model by A. L. Harris and W. Stein
! See details in BioRXIV 
!	https://biorxiv.org/cgi/content/short/2021.04.25.441350v1
! It is a 2 cell model with connected pyramidal neuron and interneuron
!  through gaba and glutamate synapses
! The pyramidal cell and interneuron reversal potentials vary
!  with ion concentrations
! Ion concentrations are updated in real-time

! This code needs the following files: ode.par, 2pi.in, const.f

!Follow these steps to compile the code
!1) compile const.f module
! 	gfortran -fdefault-real-8 -c const.f 
!2) compile 2pi.model.f code
!	mpif90 -fdefault-real-8 -c 2pi.model.f
!3) link module and code
!	mpif90 const.o 2pi.model.o -o 2pi.exe
!4) code can be executed using MPI fortran. Max processors is (Jeend - Jestart)/Jestep + 1
!	mpirun -np xx 2pi.exe
!  (xx = number of processors)

!The code produces the following output files:
!1) frequency.out.int.Jex.x.Jiy.y.gabaz.z
!	contains the firing frequency (in Hz) of the interneuron as a function of 
!	pyramidal cell injected current (Je), interneuron injected current (Ji),
! 	gaba maximum conductance (ggaba), and time
!	Files are labeled as x.x = Je value, y.y = Ji value, z.z=ggaba value
!2) frequency.out.pyr.Jex.x.Jiy.y.gabaz.z
!	same as (1), but for the interneuron
!3) if writeyn=1, membrane.pot.out.Jex.x.Jiy.y.gabaz.z
!	contains the pyramidal cell and interneuron membrane potential (mV)
!  	as a function of time - these are large files



	program twocell
	implicit none

	include 'ode.par'

!MPI	
	include 'mpif.h'
	integer num_proc,my_rank,ierr
	integer status(MPI_STATUS_SIZE)
	integer begNJe,finNJe,aveNJe,extra
!end MPI


! Declare variables
	real ystart(nvar)
	real ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),xp,yp
	real dxsav,hdid,hnext,htry,x,dydx,y,yscal,errmax,htemp,xnew
	real yerr,ytemp,SAFETY,PGROW,PSHRINK,ERRCON,dnm(6)
	integer Nt,Nx,Ncl,kmax,Nstep
	integer kount,nok,nbad,n
	character*20 real_to_char
	real ti,tf,eps,h1,dxsav_den
	external derivs
	real timei,timef,runtime,gamim,ggaba
	real Jestart,Jeend,Jestep,Je
	integer cntJe,NJe,cntJi_int,NJi_int,i,cntggaba,Nggaba
	character*3 dum,dumJi,dumgaba
	real Ji_intstart,Ji_intend,Ji_intstep,Ji_int
	real dum1,dum2,dum3,dum4
	real ggabastart,ggabaend,ggabastep
	integer writeyn
	real V0,n0,h0,Ca0,K_i0,K_o0
	real Na_i0,Na_o0,Cl_i0,Cl_o0,lils0
	real lilse0,lilsi0,V0_int,n0_int,h0_int
	real K_i0_int,Na_i0_int


	call cpu_time (timei)

!MPI
	call MPI_INIT(ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,num_proc,ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

!end MPI


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read input file								
										
	call readinput(ti,tf,ystart,Nstep,gamim,ggabastart
     &	,ggabaend,ggabastep,Jestart,Jeend,Jestep,Ji_intstart,
     &	Ji_intend,Ji_intstep,writeyn,n0,h0,Ca0,K_i0,K_o0,
     &	Na_i0,Na_o0,Cl_i0,Cl_o0,lils0,lilse0,lilsi0,V0_int,
     &	n0_int,h0_int,V0,K_i0_int,Na_i0_int)				
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! setup loop over pyramidal cell injected current Je
	NJe = nint((Jeend - Jestart)/Jestep) + 1

	aveNJe = NJe/num_proc
	extra = mod(NJe,num_proc)

	begNJe = my_rank*aveNJe+1
	if(num_proc.eq.1)then
	  finNJe = (my_rank+1)*aveNJe
	else
	  finNJe = (my_rank+1)*aveNJe+(my_rank/(num_proc-1))*extra
	endif

	
! loop over Je values - this loop is parallel	
	do cntJe = begNJe,finNJe
	  Je = Jestart + (cntJe-1)*Jestep

	  write(dum,"(F3.1)")Je

! loop over interneuron injecte current values Ji 
	  NJi_int = nint((Ji_intend - Ji_intstart)/Ji_intstep) + 1
	  do cntJi_int = 1,NJi_int
	    Ji_int = Ji_intstart + (cntJi_int -1)*Ji_intstep

	    write(dumji,"(F3.1)")Ji_int

! loop over gaba conductance ggaba
	    Nggaba = nint((ggabaend - ggabastart)/ggabastep) + 1
	    do cntggaba = 1,Nggaba
	      ggaba = ggabastart + (cntggaba - 1)*ggabastep

	      write(dumgaba,"(F3.1)")ggaba

! open/create membrane potential files
	      if(writeyn.eq.1)then
	        open(1000*cntJi_int+100*cntggaba+cntJe,file=
     &'membrane.pot.out.Je'//dum//'.Ji'//dumji//'.gaba'//dumgaba)
	        write(1000*cntJi_int+100*cntggaba+cntJe,*)           
     &'time (s)     V_pyr (mV)    V_int (mV)'
	      endif

! open/create frequency file
	      write(dum,"(F3.1)")Je
	      open(10000*cntJi_int+1000*cntggaba+10*cntJe,file=
     &'frequency.out.pyr.Je'//dum//'.Ji'//dumji//'.gaba'//dumgaba)
	      write(10000*cntJi_int+1000*cntggaba+10*cntJe,*)
     &		'Je     Ji     ggaba     time (s)      freq (Hz)'

	      open(100000*cntJi_int+10000*cntggaba+100*cntJe,file=
     &'frequency.out.int.Je'//dum//'.Ji'//dumji//'.gaba'//dumgaba)
	      write(100000*cntJi_int+10000*cntggaba+100*cntJe,*)
     &		'Je     Ji     ggaba     time (s)      freq_p (Hz)'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve ODEs with Runge-Kutta 	 	
!  See Numerical recipes												

	      call rkdumb(ystart,ti,tf,nstep,derivs,gamim,ggaba,Je
     &	    ,cntJe,Ji_int,cntJi_int,cntggaba,writeyn)										

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! if writing to membrane potential files
	      if(writeyn.eq.1)then
! write input parameters to output file
		write(1000*cntJi_int+100*cntggaba+cntJe,*)
     &		'#V0 = ',V0
		write(1000*cntJi_int+100*cntggaba+cntJe,*)
     &		'#n0,h0,Ca0 = ',n0,h0,Ca0
		write(1000*cntJi_int+100*cntggaba+cntJe,*)
     &		'#K_i0,K_o0 = ',K_i0,K_o0
		write(1000*cntJi_int+100*cntggaba+cntJe,*)
     &		'#Na_i0,Na_o0 = ',Na_i0,Na_o0
		write(1000*cntJi_int+100*cntggaba+cntJe,*)
     &		'#Cl_i0,Cl_o0 = ',Cl_i0,Cl_o0
		write(1000*cntJi_int+100*cntggaba+cntJe,*)
     &		'#lils0,lilse0,lilsi0 = ',
     &		lils0,lilse0,lilsi0
		write(1000*cntJi_int+100*cntggaba+cntJe,*)
     &		'#V0_int = ',V0_int
		write(1000*cntJi_int+100*cntggaba+cntJe,*)
     &		'#n0_int,h0_int, = ',
     &		n0_int,h0_int
		write(1000*cntJi_int+100*cntggaba+cntJe,*)
     &		'#K_i0_int,Na_i0_int = ',K_i0_int,Na_i0_int
		write(1000*cntJi_int+100*cntggaba+cntJe,*)
     &		'#Nstep = ',Nstep
		write(1000*cntJi_int+100*cntggaba+cntJe,*)
     &		'#gami = ',gamim 

	        close(1000*cntJi_int+100*cntggaba+cntJe)
	      endif

! write input parameters to output file pyramidal
		write(10000*cntJi_int+1000*cntggaba+10*cntJe,*)
     &		'#V0 = ',V0
		write(10000*cntJi_int+1000*cntggaba+10*cntJe,*)
     &		'#n0,h0,Ca0 = ',n0,h0,Ca0
		write(10000*cntJi_int+1000*cntggaba+10*cntJe,*)
     &		'#K_i0,K_o0 = ',K_i0,K_o0
		write(10000*cntJi_int+1000*cntggaba+10*cntJe,*)
     &		'#Na_i0,Na_o0 = ',Na_i0,Na_o0
		write(10000*cntJi_int+1000*cntggaba+10*cntJe,*)
     &		'#Cl_i0,Cl_o0 = ',Cl_i0,Cl_o0
		write(10000*cntJi_int+1000*cntggaba+10*cntJe,*)
     &		'#lils0,lilse0,lilsi0 = ',
     &		lils0,lilse0,lilsi0
		write(10000*cntJi_int+1000*cntggaba+10*cntJe,*)
     &		'#V0_int = ',V0_int
		write(10000*cntJi_int+1000*cntggaba+10*cntJe,*)
     &		'#n0_int,h0_int = ',
     &		n0_int,h0_int
		write(10000*cntJi_int+1000*cntggaba+10*cntJe,*)
     &		'#K_i0_int,Na_i0_int = ',K_i0_int,Na_i0_int
		write(10000*cntJi_int+1000*cntggaba+10*cntJe,*)
     &		'#Nstep = ',Nstep
		write(10000*cntJi_int+1000*cntggaba+10*cntJe,*)
     &		'#gami = ',gamim 

	      close(10000*cntJi_int+1000*cntggaba+10*cntJe)

! write input parameters to output file interneuron
		write(100000*cntJi_int+10000*cntggaba+100*cntJe,*)
     &		'#V0 = ',V0
		write(100000*cntJi_int+10000*cntggaba+100*cntJe,*)
     &		'#n0,h0,Ca0 = ',n0,h0,Ca0
		write(100000*cntJi_int+10000*cntggaba+100*cntJe,*)
     &		'#K_i0,K_o0 = ',K_i0,K_o0
		write(100000*cntJi_int+10000*cntggaba+100*cntJe,*)
     &		'#Na_i0,Na_o0 = ',Na_i0,Na_o0
		write(100000*cntJi_int+10000*cntggaba+100*cntJe,*)
     &		'#Cl_i0,Cl_o0 = ',Cl_i0,Cl_o0
		write(100000*cntJi_int+10000*cntggaba+100*cntJe,*)
     &		'#lils0,lilse0,lilsi0 = ',
     &		lils0,lilse0,lilsi0
		write(100000*cntJi_int+10000*cntggaba+100*cntJe,*)
     &		'#V0_int = ',V0_int
		write(100000*cntJi_int+10000*cntggaba+100*cntJe,*)
     &		'#n0_int,h0_int ',
     &		n0_int,h0_int
		write(100000*cntJi_int+10000*cntggaba+100*cntJe,*)
     &		'#Nstep = ',Nstep
		write(100000*cntJi_int+10000*cntggaba+100*cntJe,*)
     &		'#gami = ',gamim 

	      close(100000*cntJi_int+10000*cntggaba+100*cntJe)

	    enddo	!ggaba loop

	  enddo		!Ji_int loop

	enddo	!Je loop

	call cpu_time (timef)
      	runtime=(timef-timei)
	print*,'total runtime was ',runtime,' seconds'

! MPI
	call MPI_FINALIZE(ierr)
!end MPI

	stop
	end program twocell



! Subroutines (mostly) in order that they are called

! This routine contains the differential equations to be integrated
! modified from Numerical Recipes for the model here
*=====================================================================*
	subroutine derivs(x,y,dydx,gamim,ggaba,Je,Ji_int)
*	These are the derivatives of the original functions, which are
*	ystart(i)
*=====================================================================*

	use const
	implicit none

	include 'ode.par'

	dimension y(nmax),dydx(nmax)	
	real x,y,dydx
	real V,n,h,Ca,m
	real IL,Ik,INa,INap,IAHP
	real an,bn,ah,bh,am,bm
	real ma1,minf
	real EL,Ek,ENa,ECl
	real c_k_i,c_k_o,c_Na_i,c_Na_o,c_Cl_i,c_Cl_o
	real IkL,Ipump
	real IKcc,INKcc,IKi,Isink
	real IClL,Igaba,lils
	real INaL
	real se,si,Iglute,Igluti

	real V_int,n_int,h_int,Ca_int
	real IL_int,Ik_int,INa_int,IAHP_int
	real an_int,bn_int,ah_int,bh_int
	real ma1_int
	real ECa_int,ECl_int
	real m_int,INaL_int,IKL_int,IClL_int,minf_int
	real am_int,bm_int
	real gamim,gami,ggaba,Je,Ji_int
	real Ek_int,ENa_int
	real c_k_i_int,c_Na_i_int,Ipump_int
	real ICa

	gami = gamim*gam

!=======USER SUPPLIED==================================
! these are the coupled DE's for the pyramidal neuron
!  and the concentrations

	V = y(1)
	n = y(2)
	h = y(3)
	Ca = y(4)

! concentrations are labeled with c_xxx_i or c_xxx_o
! the i/o refers to inside or outside, xxx is the ion

	c_k_i = y(5)
	c_k_o = y(6)
	c_Na_i = y(7)
	c_Na_o = y(8)
	c_Cl_i = y(9)
	c_Cl_o = y(10)
	lils = y(11)
	se = y(12)
	si = y(13)


	Ek = RTF*log(c_k_o/c_k_i)
	ENa = RTF*log(c_Na_o/c_Na_i)
	ECl = -RTF*log(c_Cl_o/c_Cl_i)

	dydx(11) = -lils/taugaba
	Igaba = ggaba*lils*(V - ECl)

	dydx(12) = -se/tauglut
	dydx(13) = -si/tauglut

	IkL = gkL*(V-Ek)
	IClL = gClL*(V - ECl)
	INaL = gNaL*(V - ENa)

	IL = IkL + IClL + INaL

	an = 0.032*(V+52.)/(1. - exp(-(V+52.)/5.))
	bn = 0.5*exp(-(V+57.)/40.)

	ah = 0.128*exp(-(V+50.)/18.)
	bh = 4./(1. + exp(-(V+27.)/5.))
	
	am = 0.32*(V+ 54.)/(1. - exp(-(V+ 54.)/4.))
	bm = .28*(V+27)/(exp((V+ 27.)/5.)-1.)

	ma1 = 1./(1. + exp(-(V+25.)/2.5))

	minf = am/(am+bm)
	m = minf

	Ik = gk*n**4.*(V - Ek)
	INa = gNa*m**3.*h*(V - ENa)
	INap = gp*m**3.*(V - ENa)
	IAHP = gAHP*(Ca/(Ca + 1.))*(V - Ek)
	ICa = gCa*ma1*(V - ECa)

	Ipump = (rhopump/gam)/( (1.+exp((Nasat-c_Na_i)/3.))
     &	  *(1.+exp(Ksat - c_K_o)) )
	IKcc = rhokcc*log(c_K_i*c_Cl_i/(c_K_o*c_Cl_o))
	INKcc = rhoNKcc*( log(c_K_i*c_Cl_i/(c_K_o*c_Cl_o))
     &	  + log(c_Na_i*c_Cl_i/(c_Na_o*c_Cl_o)) )
     &	  /(1.+exp(16. - c_K_o))
	Isink = epsK*(c_K_o - Kbath)
	Iglute = gglut1*se*(V - Eglut)


	dydx(5) = -tauinv*(gam*(Ik + IAHP + IKL - 2.*Ipump) 
     &	  + IKcc + INKcc)
	dydx(7) = tauinv*( -gam*(INa + INap + INaL + 3.*Ipump) 
     &	  - INKcc )
	dydx(9) = tauinv*( gam*(Igaba + IClL) - IKcc - 2.*INKcc )
	dydx(10) = -beta*dydx(9)


	dydx(1) = (Je - IL - Ik - INa - INap - Igaba
     &	  - Iglute - IAHP - Ipump)/C
	dydx(2) = phi*(an*(1.-n) - bn*n)
	dydx(3) = phi*(ah*(1.-h) - bh*h)
	dydx(4) = -epsi*gCa*ma1*(V - ECa) - Ca/tauCa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! these are the coupled DE's for the interneuron and concentrations
	V_int = y(14)
	n_int = y(15)
	h_int = y(16)
	c_k_i_int = y(17)
	c_Na_i_int = y(18)

	Ek_int = RTF*log(c_k_o/c_k_i_int)
	ENa_int = RTF*log(c_Na_o/c_Na_i_int)

! USER - if you want to write the ion concentrations as a function of time
!	uncomment the next line - the files are large

!	write(27,*)x/1000.,c_k_i,c_Na_i,c_k_o,c_Na_o,c_k_i_int,c_Na_i_int

	am_int = 0.1*(V_int + 35.)/(1. - exp(-(V_int + 35.)/10.))
	bm_int = 4.*exp(-(V_int + 60.)/18.)
	ah_int = 0.07*exp(-(V_int + 58.)/20.)
	bh_int = 1./(1. + exp(-(V_int + 28.)/10.))
	an_int = 0.01*(V_int + 34.)/(1. - exp(-(V_int + 34.)/10.))
	bn_int = 0.125*(exp(-(V_int + 44.)/80.))

	minf_int = am_int/(am_int + bm_int)
	m_int = minf_int

	IkL_int = gkL_int*(V_int - Ek_int)
	INaL_int = gNaL_int*(V_int - ENa_int)
	IL_int = IkL_int + INaL_int
	
	Ipump_int = (rhopump/gami)/( (1.+exp((Nasat-c_Na_i_int)/3.))
     &	  *(1.+exp(Ksat - c_K_o)) )
	Ik_int = gk_int*n_int**4.*(V_int - Ek_int)
	INa_int = gNa_int*m_int**3.*h_int*(V_int - ENa_int)

	Igluti = gglut2*si*(V_int - Eglut)

	dydx(14) = (Ji_int - IL_int - Ik_int - INa_int 
     &	  - Igluti - Ipump_int)/C_int

	dydx(15) = phi_int*(an_int*(1.-n_int) - bn_int*n_int)
	dydx(16) = phi_int*(ah_int*(1.-h_int) - bh_int*h_int)

	dydx(17) = -tauinv*(gami*(Ik_int + IKL_int - 2.*Ipump_int) )
	dydx(18) = tauinv*(-gami*(INa_int + INaL_int + 3.*Ipump_int))

	dydx(6) = tauinv*(gam*beta*(Ik + IAHP + IkL - 2.*Ipump)
     &	  + beta*(IKcc + INKcc) - Isink) - beta*dydx(17) 		
	dydx(8) = -beta*(dydx(7) + dydx(18))


!=======END USER SUPPLIED===============================


	return
	end



!-----------------------------------------------------------------------
! This function converts a real number into a character, so that it can be
! used for naming files, etc.

	function real_to_char(angle)
!-----------------------------------------------------------------------
		
	real angle
	character*20 cangle
	character*10 frmt
	character*20 real_to_char
	

	if(angle.gt.9999)then
	   print*,'angle is too big for function real_to_char'
	   stop	   
	else if(angle.lt.-9999)then
	   print*,'angle is too small for function real_to_char'
	   stop	
	else if(angle.eq.0)then
	   frmt='(F7.5)'
	else if(abs(angle).lt.10)then
	   frmt='(F7.5)'
	else
	   frmt='(F7.5)'
	endif

	write(cangle,frmt)angle

	real_to_char=cangle

	end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Reads input file
	subroutine readinput(ti,tf,ystart,Nstep,gamim,ggabastart,
     &	ggabaend,ggabastep,Jestart,Jeend,Jestep,Ji_intstart,
     &	Ji_intend,Ji_intstep,writeyn,n0,h0,Ca0,K_i0,K_o0,
     &	Na_i0,Na_o0,Cl_i0,Cl_o0,lils0,lilse0,lilsi0,V0_int,
     &	n0_int,h0_int,V0,K_i0_int,Na_i0_int)

	implicit none
	include 'ode.par'

	real ti,tf,V0,n0,h0,Ca0,K_i0,K_o0
	real Na_i0,Na_o0,Cl_i0,Cl_o0,lils0
	integer Nstep,writeyn
	real ystart(nvar)
	real lilse0,lilsi0,V0_int,n0_int,h0_int
	real gamim,Jestart,Jeend,Jestep
	real Ji_intstart,Ji_intend,Ji_intstep
	real ggabastart,ggabaend,ggabastep
	real K_i0_int,Na_i0_int

	open(2,file='2pi.in')
	
	read(2,*)ti
	read(2,*)tf
	read(2,*)V0
	read(2,*)n0,h0,Ca0
	read(2,*)K_i0,K_o0
	read(2,*)Na_i0,Na_o0
	read(2,*)Cl_i0,Cl_o0
	read(2,*)lils0,lilse0,lilsi0
	read(2,*)V0_int
	read(2,*)n0_int,h0_int
	read(2,*)K_i0_int
	read(2,*)Na_i0_int
	read(2,*)Nstep
	read(2,*)gamim
	read(2,*)ggabastart,ggabaend,ggabastep
	read(2,*)Jestart,Jeend,Jestep
	read(2,*)Ji_intstart,Ji_intend,Ji_intstep
	read(2,*)writeyn

	if(nstep.gt.nstpmx)then
	  print*,'nstep > nstpmax'
	  print*,'increase nstpmax'
	  stop
	endif

! put initial conditions into array ystart
	ystart(1) = V0
	ystart(2) = n0
	ystart(3) = h0
	ystart(4) = Ca0

	ystart(5) = K_i0
	ystart(6) = K_o0
	ystart(7) = Na_i0
	ystart(8) = Na_o0
	ystart(9) = Cl_i0
	ystart(10) = Cl_o0

	ystart(11) = lils0
	ystart(12) = lilse0	
	ystart(13) = lilsi0

	ystart(14) = V0_int
	ystart(15) = n0_int
	ystart(16) = h0_int

	ystart(17) = K_i0_int
	ystart(18) = Na_i0_int

	return
	end

! From Numerical Recipes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs,gamim,ggaba,Je,
     &	  Ji_int)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	INTEGER n,NMAX
	REAL h,x,dydx(n),y(n),yout(n)
	EXTERNAL derivs
	PARAMETER (NMAX=50) 	!Set to the maximum number of functions.
!Given values for the variables y(1:n) and their derivatives dydx(1:n) known at x, use
!the fourth-order Runge-Kutta method to advance the solution over an interval h and return
!the incremented variables as yout(1:n), which need not be a distinct array from y. The
!user supplies the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
	INTEGER i
	REAL h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)
	real gamim,ggaba,Je,Ji_int

	hh=h*0.5
	h6=h/6.
	xh=x+hh

	do i=1,n 	!First step.
	  yt(i)=y(i)+hh*dydx(i)
	enddo 

	call derivs(xh,yt,dyt,gamim,ggaba,Je,Ji_int) 	!Second step.

	do i=1,n
	  yt(i)=y(i)+hh*dyt(i)
	enddo 

	call derivs(xh,yt,dym,gamim,ggaba,Je,Ji_int) 	!Third step.

	do i=1,n
	  yt(i)=y(i)+h*dym(i)
	  dym(i)=dyt(i)+dym(i)
	enddo 

	call derivs(x+h,yt,dyt,gamim,ggaba,Je,Ji_int) 	!Fourth step.

	do i=1,n 	!Accumulate increments with proper weights.
	  yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
	enddo 

	return
	END

! From Numerical Recipes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE rkdumb(vstart,x1,x2,nstep,derivs,gamim,ggaba,Je,
     &	cntJe,Ji_int,cntJi_int,cntggaba,writeyn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	include 'ode.par'

	INTEGER nstep
!Maximum number of functions and
!maximum number of values to be stored.
	REAL x1,x2,vstart(nvar),xx(NSTPMX),y(NMAX,NSTPMX)
	EXTERNAL derivs
!	COMMON /path/ xx,y 	!Storage of results.
! USES rk4
!Starting from initial values vstart(1:nvar) known at x1 use fourth-order Runge-Kutta to
!advance nstep equal increments to x2. The user-supplied subroutine derivs(x,v,dvdx)
!evaluates derivatives. Results are stored in the common block path. Be sure to dimension
!the common block appropriately.
	INTEGER i,k,j,cntJe,cntJi_int,cntggaba,writeyn,j_int
	REAL h,x,dv(NMAX),v(NMAX)
	real dVp(NSTPMX),dVpp,freq,tspike(NSTPMX),dVp_int(NSTPMX)
	real gamim,ggaba,Je,Ji_int,dVpp_int,tspike_int(NSTPMX)
	real freq_int
	integer printme

	do i=1,nvar 	!Load starting values.
	  v(i)=vstart(i)
	  y(i,1)=v(i)
	enddo 

	xx(1)=x1
	x=x1
	h=(x2-x1)/nstep	


! write the functions to the output file

	j = 1
	printme = 1

	do k=1,nstep 	!Take nstep steps.

	  if(v(14) .gt. 0.)then		!v(14) is interneuron V
	    v(11) = 1.
	  endif

	  if(v(1) .gt. 0.)then		!v(1) is pyramidal neuron V
	    v(12) = 1.
	    v(13) = 1.
	  endif

	  call derivs(x,v,dv,gamim,ggaba,Je,Ji_int)
	  call rk4(v,dv,nvar,x,h,v,derivs,gamim,ggaba,Je,Ji_int)

	  if(x+h.eq.x)then
	    print*,'stepsize not significant in rkdumb'
	    stop
	  endif

	  x=x+h
	  xx(k+1)=x 	!Store intermediate steps.

	  do i=1,nvar
	    y(i,k+1)=v(i)
	  enddo 

!write membrane potentials to file
	  if(writeyn.eq.1)then
	    write(1000*cntJi_int+100*cntggaba+cntJe,'(6(f6.2,3x))')
     &	x/1000.,y(1,k),y(14,k)		!time,V_pyramid,V_int
	  endif

!calculate frequency as a function of time and write to file
! 1st and 2nd derivs of V for pyramidal neuron
	  if(k.gt.2)then
	    dVp(k) = (y(1,k) - y(1,k-1))/h
	    dVpp = (y(1,k) - 2.*y(1,k-1) + y(1,k-2))/h**2.

	    dVp_int(k) = (y(14,k) - y(14,k-1))/h
	    dVpp_int = (y(14,k) - 2.*y(14,k-1) + y(14,k-2))/h**2.
	  endif

	  if(k.gt.3)then
	    if(dVp(k).lt.0. .and. dVp(k-1).gt.0. .and. dVpp.lt.0. 	!1st deriv changes from + to -
     &	    .and. y(1,k).gt.0.)then					!and 2nd deriv<0 and V_pyramid>0
	      tspike(j) = x/1000.
	
	      if(j.gt.1)then
		freq = 1./(tspike(j) - tspike(j-1))
		write(10000*cntJi_int+1000*cntggaba+10*cntJe,'(5(f6.2,3x))')
     &	Je,Ji_int,ggaba,tspike(j-1),freq

	      endif		!j.gt.1

	      j = j+1
	    endif		!if a spike happens

! interneuron frequency
	    if(dVp_int(k).lt.0. .and. dVp_int(k-1).gt.0.  	!1st deriv changes from + to -
     &	    .and. dVpp_int.lt.0. .and. y(14,k).gt.0.)then	!and 2nd deriv<0 and V_int>0
	      tspike_int(j_int) = x/1000.
	
	      if(j_int.gt.1)then
		freq_int = 1./(tspike_int(j_int) - tspike_int(j_int-1))
		write(100000*cntJi_int+10000*cntggaba+100*cntJe,'(5(f6.2,3x))')
     &	Je,Ji_int,ggaba,tspike_int(j_int-1),freq_int

	      endif	!j_int.gt.1

	      j_int = j_int+1
	    endif	!if a spike happens

	  endif		!k.gt.3	      

	enddo 


	return
	END


