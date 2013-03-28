#!/usr/bin/perl

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

@Edep = (0.3,1.0,3.0,10.0,30.0,100.0);

foreach $ener (@Edep) {
	system "cp init/init.dat.default init.dat";
	system "echo \"Bfield\t1e14\n\" >>init.dat";
	system "echo \"angle_mu\t1\n\" >>init.dat";
	system "echo \"Edep\t$ener\n\" >>init.dat";
	if ($ener > 0.3) {
		system "echo \"precalc\t0\n\" >>init.dat";	
	}
	sleep 3;
	system "crustcool";
	sleep 10;
	system "mv gon_out/prof gon_out/prof_B1e14E".$ener."_1e9_new";
	sleep 3;
}

# Make a dummy gon_out/prof so that the plotting routine works!
system "cp gon_out/prof_1822_A gon_out/prof"
