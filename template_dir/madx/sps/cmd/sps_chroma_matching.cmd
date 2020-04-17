match_chroma(QPH_val, QPV_val) : macro ={
	
	value QPH_val;
	value QPV_val;
	! match chromaticity
	match, sequence=sps;
	vary, name=QPH;
	vary, name=QPV;
	global, dq1=QPH_val, dq2=QPV_val;
	jacobian, calls=1000, tolerance=1e-20;
	endmatch;
	twiss;

	

};
