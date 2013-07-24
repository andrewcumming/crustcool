#!/usr/bin/perl


# Bfield
if (0) {
	@Bfield = (1.8e13,5.6e13,5.6e14);
#	@Bfield = (1e13,3e13,1e14,1.8e14,3e14,1e15,1.8e15);
	
	foreach $BB (@Bfield) {

		system "cp init/init.dat.default init.dat";

		system "echo \"Bfield\t".$BB."\n\" >>init.dat";
		system "echo \"angle_mu\t1\n\" >>init.dat";
		system "echo \"Tc\t5e7\n\" >>init.dat";
		system "echo \"precalc\t1\n\" >>init.dat";	
		system "echo \"Edep\t30\n\" >>init.dat";	

		sleep 3;
		system "crustcool";
		sleep 10;

		$command = sprintf("mv gon_out/prof gon_out/prof_B%g_E30_Tc5e7_mu1",$BB);
		system $command;
		print "$command\n";
		sleep 3;
	}
	
}



# Crust heating at different depths
#
if (1) {

#	@rhovec = (3e9,1e10,3e10,1e11,3e11);
	@rhovec = (1e8,3e8,1e9,3e9,1e10,3e10,1e11,3e11);

foreach $rho (@rhovec) {
	system "cp init/init.dat.default init.dat";
	system "echo \"Bfield\t1e14\n\" >>init.dat";
#	system "echo \"angle_mu\t-1\n\" >>init.dat";
	system "echo \"angle_mu\t1\n\" >>init.dat";
	system "echo \"Tc\t5e7\n\" >>init.dat";
	if ($rho > 1e8) {
		system "echo \"precalc\t0\n\" >>init.dat";	
	}

	system "echo \"piecewise\t1\" >>init.dat";
	system "echo \">0\t1e9\" >>init.dat";
	system "echo \">".$rho."\t1e9\" >>init.dat";
	$rho1=1.01*$rho;
	system "echo \">".$rho1."\t-1\" >>init.dat";
	system "echo \">-1\t-1\n\" >>init.dat";

	sleep 3;
	system "crustcool";
	sleep 10;
	# these ones have Tc1e8
#		$command = sprintf("mv gon_out/prof gon_out/prof_new_B1e14_T2e9_rho%g",$rho);
#		$command = sprintf("mv gon_out/prof gon_out/prof_new_B1e14_T2e9_rho%g_mu1",$rho);
#		$command = sprintf("mv gon_out/prof gon_out/prof_new_B1e14_T2e9_Tc5e7_rho%g",$rho);
$command = sprintf("mv gon_out/prof gon_out/prof_new_B1e14_T1e9_Tc5e7_rho%g_mu1_env",$rho);
#$command = sprintf("mv gon_out/prof gon_out/prof_new_B1e14_T2e9_Tc1e8_rho%g_mu1_env",$rho);
	system $command;
	print "$command\n";
#	system "mv gon_out/prof gon_out/prof_B1e14E".$rho."_mu1";
	sleep 3;
}
}



# Makes cooling curves with mu=-1, and then 0.1 to 1.0 in steps of 0.1
# Used to compute net lightcurve from the whole surface
if (0) {

$suffix = "rn";

for ( $i=0; $i<=10; $i++) {

	if ($i==0) {
		$mu = -1.0; 
	} else {
		$mu = 0.1*$i;
	}
	system "cp init/init.dat.steps init/init.dat.$suffix";
	system "echo \"angle_mu\t$mu\n\" >>init/init.dat.$suffix";
	system "echo \"Tc\t1e8\n\" >>init/init.dat.$suffix";
	sleep 3;
	system "crustcool $suffix";
	sleep 10;
	if ($i>0) {
		system "mv gon_out/prof gon_out/prof_".$suffix."_mu$mu";
	} else {
		system "mv gon_out/prof gon_out/prof_$suffix";
	}
	sleep 3
}

system "mv gon_out/prof_".$suffix."_mu1 gon_out/prof_".$suffix."_mu1.0";

}



if (0)
 {
	$ener=10.0;
	system "cp init/init.dat.default init.dat";
	system "echo \"Bfield\t1e14\n\" >>init.dat";
	system "echo \"angle_mu\t1\n\" >>init.dat";
	system "echo \"Edep\t$ener\n\" >>init.dat";
	system "echo \"Tc\t2e8\n\" >>init.dat";
#	if ($ener > 0.3) {
#		system "echo \"precalc\t0\n\" >>init.dat";	
#	}
	sleep 3;
	system "crustcool";
#	sleep 10;
#	system "mv gon_out/prof gon_out/prof_new_B1e14E".$ener."_1e9_mu1";
#	sleep 3;
}





# Different energies deposited in the outer crust
# Used to make a comparison figure with Pons
if (0) {

@Edep = (0.3,1.0,3.0,10.0,30.0,100.0);

foreach $ener (@Edep) {
	system "cp init/init.dat.default init.dat";
	system "echo \"Bfield\t1e14\n\" >>init.dat";
	system "echo \"angle_mu\t1\n\" >>init.dat";
	system "echo \"Edep\t$ener\n\" >>init.dat";
	system "echo \"Tc\t1e8\n\" >>init.dat";
	if ($ener > 0.3) {
		system "echo \"precalc\t0\n\" >>init.dat";	
	}
	sleep 3;
	system "crustcool";
	sleep 10;
	system "mv gon_out/prof gon_out/prof_B1e14E".$ener."_1e9_mu1_env";
	system "mv gon_out/out gon_out/out_B1e14E".$ener."_1e9_mu1_env";
	sleep 3;
}

foreach $ener (@Edep) {
	system "cp init/init.dat.default init.dat";
	system "echo \"Bfield\t1e15\n\" >>init.dat";
	system "echo \"angle_mu\t1\n\" >>init.dat";
	system "echo \"Edep\t$ener\n\" >>init.dat";
	system "echo \"Tc\t1e8\n\" >>init.dat";
	if ($ener > 0.3) {
		system "echo \"precalc\t0\n\" >>init.dat";	
	}
	sleep 3;
	system "crustcool";
	sleep 10;
	system "mv gon_out/prof gon_out/prof_B1e15E".$ener."_1e9_mu1_env";
	system "mv gon_out/out gon_out/out_B1e15E".$ener."_1e9_mu1_env";
	sleep 3;
}


}



if (0) {
foreach $ener (@Edep) {
	system "cp init/init.dat.default init.dat";
	system "echo \"Bfield\t1e14\n\" >>init.dat";
	system "echo \"angle_mu\t-1.0\n\" >>init.dat";
	system "echo \"Edep\t$ener\n\" >>init.dat";
	system "echo \"Tc\t2e8\n\" >>init.dat";
	if ($ener > 0.3) {
		system "echo \"precalc\t0\n\" >>init.dat";	
	}
	sleep 3;
	system "crustcool";
	sleep 10;
	system "mv gon_out/prof gon_out/prof_B1e14E".$ener."_1e9";
	system "mv gon_out/out gon_out/out_B1e14E".$ener."_1e9";
	sleep 3;
}

foreach $ener (@Edep) {
	system "cp init/init.dat.default init.dat";
	system "echo \"Bfield\t1e15\n\" >>init.dat";
	system "echo \"angle_mu\t1\n\" >>init.dat";
	system "echo \"Edep\t$ener\n\" >>init.dat";
	system "echo \"Tc\t2e8\n\" >>init.dat";
	if ($ener > 0.3) {
		system "echo \"precalc\t0\n\" >>init.dat";	
	}
	sleep 3;
	system "crustcool";
	sleep 10;
	system "mv gon_out/prof gon_out/prof_B1e15E".$ener."_1e9_mu1";
	system "mv gon_out/out gon_out/out_B1e15E".$ener."_1e9_mu1";
	sleep 3;
}

foreach $ener (@Edep) {
	system "cp init/init.dat.default init.dat";
	system "echo \"Bfield\t1e15\n\" >>init.dat";
	system "echo \"angle_mu\t-1.0\n\" >>init.dat";
	system "echo \"Edep\t$ener\n\" >>init.dat";
	system "echo \"Tc\t2e8\n\" >>init.dat";
	if ($ener > 0.3) {
		system "echo \"precalc\t0\n\" >>init.dat";	
	}
	sleep 3;
	system "crustcool";
	sleep 10;
	system "mv gon_out/prof gon_out/prof_B1e15E".$ener."_1e9";
	system "mv gon_out/out gon_out/out_B1e15E".$ener."_1e9";
	sleep 3;
}
}

# Make a dummy gon_out/prof so that the plotting routine works!
system "cp gon_out/prof_B1e14E1_1e9 gon_out/prof"
