plot_QxQyVSdeltap : macro ={
 	call, file='ptc/PTC.macro';

	create, table=PTCchroma, column=dp,Qx0,Qy0;

	dpmax=0.8e-2; !10e-3; !0.6e-2; !10e-3; !0.7e-2;
 	nmax=10;
 	n=-nmax;
	dp:=dpmax*n/nmax;
 	!i could have set up the dpmin instead. If I don't, the way the loop is set, it will go from -dpmax to +dpmax
	while (n <= nmax)!I could do while dp<blah then...
 		{
         		value, dp,Qx0,Qy0;
          		exec, PTCtunes;
         		fill, table=PTCchroma;

         		n=n+1;
 		};

	
	!write, table=PTCchroma;

	write, table=PTCchroma, file="PTC-QxQydp.ptc";

	plot, table=PTCchroma, haxis=dp,hmin=-1.0e-2,hmax=1.0e-2, vaxis=qx0,qy0, vmin=0.08,vmax=0.24 title='PTC - nonlinear chromaticity - dpmax0.8e-2',noversion=true, file=nonlinear_chroma_sps_multipoles;


};
