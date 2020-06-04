 option, -echo;

/***************************************************
 * SPECIFY THE SPS BEAM CONFIGURATION, ENERGY and EMITTANCES
 ***************************************************/

 call, file = 'sps/beams/lhc_beam_injection.beamx';

/***************************************************
 * files for SPS
 ***************************************************/

 call, file = 'sps/elements/sps.ele';
 call, file = 'sps/aperture/aperturedb_1.dbx';
 call, file = 'sps/aperture/aperturedb_2.dbx';
 call, file = 'sps/sequence/sps.seq';
 !call, file = 'sps/aperture/aperturedb_3.dbx';
 call, file = 'sps/strength/SPS_Q20.str';
 call, file = 'sps/strength/elements.str';
 option, echo;


/***************************************************
 * use
 ***************************************************/

 use, sequence=SPS;

/***************************************************
 * Build thin model
 ***************************************************/ 
 myvalue(xx,yy,zz): macro = {myval = table(xx,yy,zz);};


 freq = 400;
 CRAVITY.1 : RFMULTIPOLE, VOLT=0, FREQ=freq, TILT=0, KNL:={knl1_cravity1}, PNL:={pnl1_cravity1};
 CRAVITY.2 : RFMULTIPOLE, VOLT=0, FREQ=freq, TILT=0, KNL:={knl1_cravity2}, PNL:={pnl1_cravity2};

 ! install crab cavities
 ! remove markers and instruments
 USE, period=SPS, range=#S/#E;
 select, flag=seqedit, class=instrument;
 select, flag=seqedit, class=marker;
 seqedit, sequence=SPS;
	remove, element=BEGI.10010;		!zero length element
	remove, element=VVFB.21801;		!zero length element
	remove, element=VVFB.21877;		!zero length element
	remove, element=QSPL.31809;		!zero length element
	remove, element=VVFB.61801;		!zero length element
	remove, element=QSPL.61809;		!zero length element
	remove, element=VVFB.61877;		!zero length element
	remove, element=selected;
	install, element=CRAVITY.1, at=6312.7213;
	install, element=CRAVITY.2, at=6313.3213;
	flatten;
 endedit;
 USE, period=SPS, range=#S/#E;


 ! START AT BWS51995 (original s=5243.0323)
 mystart: marker;
 SEQEDIT, sequence=SPS;
	install, element=mystart, at=5243.0323; !equivalent to: at=BWSA.51995->l/2, from=BWSA.51995;
	flatten;
	cycle,start=mystart;
	remove, element=SPS$START;
	remove, element=SPS$END;
	flatten;
 ENDEDIT;
 USE, sequence=SPS;


 ! make thin 
 use, sequence=sps;
 select, flag=makethin, slice=1, thick=false;
 makethin, sequence=sps, style=teapot, makedipedge=true;
 use, sequence=SPS;
 twiss;

 ! match tunes
 match, sequence=sps;
 vary, name=kqd, step=1e-8;
 vary, name=kqf1, step=1e-8;
 global, q1=26.13, q2=26.18;
 jacobian,calls=1000,tolerance=1e-12;
 endmatch;

 twiss;

option, -warn;
save, sequence=sps, beam=true, file=sps_thin.seq;
