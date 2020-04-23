
option,echo,warn;

AssignMultipoles : macro = {

// Dipole Errors
	EOPTION,SEED=98634628;
	Select, flag=error, clear = true;
//	Select,flag=error,class=RBEND, pattern="^MBA";
	Select,flag=error, class=multipole, pattern="^MBA";
	EFCOMP, dkn={0,0,b3a+b3ar*tgauss(3),b4a+b4ar*tgauss(3),b5a+b5ar*tgauss(3),b6a,b7a+b7ar*tgauss(3),0,b9a,0,0,0,0,0,0,0};
	
	Select, flag=error, clear = true;
	!Select,flag=error,class=RBEND, pattern="^MBB";
	Select,flag=error,class=multipole, pattern="^MBB";
	EFCOMP, dkn={0,0,b3b+b3br*tgauss(3),b4b+b4br*tgauss(3),b5b+b5br*tgauss(3),b6b,b7b+b7br*tgauss(3),0,b9b,0,0,0,0,0,0,0};
!        value, b3a+b3ar*tgauss(3);

value, KQD, KQF1, b3a, b3b, b4a, b4b, b4f, b4fa, b4d, kLOF, kLOD, b5a, b5b, b6f, b6d, b7a, b7b;
fill, table=fittedmultipoles;
write, table=fittedmultipoles, file="Fittedmultipoles.txt";
*/
};
