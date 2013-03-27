#!/usr/bin/perl

$suffix = "rn2";

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
