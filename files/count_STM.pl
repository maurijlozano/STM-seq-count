#! /usr/bin/env perl
# Count-STM - Programed by Mauricio J Lozano
use strict;
use warnings;

# Open files for barcodes and pools, and loads files line by line into arrays.
open (BARCODES, "<barcodes.txt");
my @barcode_a = <BARCODES>;
close (BARCODES);

open (POOLS, "<poolsk.txt");
my @poolsk_a = <POOLS>;
close (POOLS);

open (POOLS, "<poolsh.txt");
my @poolsh_a = <POOLS>;
close (POOLS);

open (FIRMAS, "<firmask.txt");
my @firmask_a = <FIRMAS>;
close (FIRMAS);

open (FIRMAS, "<firmash.txt");
my @firmash_a = <FIRMAS>;
close (FIRMAS);

# Declare variables for hashes namespace-package.
my %condicion_h;
my %poolsk_h;
my %poolsh_h;
my %firmask_h;
my %firmash_h;
my %resultados_h;

# Transforms barcodes, pools and signature tag into hashes.
# Slits each element in the arrays separated by comma, removes new line characters and stores into hash.
foreach (@barcode_a) {
	my @barcode_a = split(",", $_);
	# convertimos en hash, primero removemos las newline
	chomp(@barcode_a);
	$condicion_h{$barcode_a[0]} = $barcode_a[1];
}

# pooslk
foreach (@poolsk_a) {
	my @poolk_a = split(",", $_);
	chomp(@poolk_a);
	# convertimos en hash
	$poolsk_h{$poolk_a[0]} = $poolk_a[1];
}

# pooslh
foreach (@poolsh_a) {
	my @poolh_a = split(",", $_);
	chomp(@poolh_a);
	# convertimos en hash
	$poolsh_h{$poolh_a[0]} = $poolh_a[1];
}

# firmask -signature tags
foreach (@firmask_a) {
	my @firmak_a = split(",", $_);
	chomp (@firmak_a);
	# convertimos en hash
	$firmask_h{$firmak_a[0]} = $firmak_a[1];
}

# firmash -signature tags
foreach (@firmash_a) {
	my @firmah_a = split(",", $_);
	chomp (@firmah_a);
	# convertimos en hash
	$firmash_h{$firmah_a[0]} = $firmah_a[1];
}

#Opens and loads "path...fasta fileq" line by line and begins analysis
open(STMSEQ, "</media/mauricio/Datos1/STM_seq.fasta");

while (<STMSEQ>) {
	my $seq = $_;
	chomp ($seq);
	# Skips fasta file headers begining with > 
	eval {	
	if ($seq =~ /^>/){		
	}
	else {
		# Searchs barcodes, pools, restriction enzime recognition sequence and mutant-tag and counts.
		# Stores into hash a pair of condition/pool/mutant-tag:count
		foreach my $key_c (keys %condicion_h){
			my $revcom_keyc = revcom($condicion_h{$key_c});
			if ($seq =~ /^(.){0,5}($condicion_h{$key_c})((.){85,90})$/p){
				my $seq1 = $3;				
				if ($key_c =~ /.*-K$/){
				foreach my $key_p (keys %poolsk_h){
					my $revcom_keyp = revcom($poolsk_h{$key_p});	
					if ($seq =~ /^.{85,90}$revcom_keyp(.){0,6}$/){
						foreach my $key_f (keys %firmask_h){
							my $revcom_keyf = revcom($firmask_h{$key_f});
							if ($seq1 =~ /^(.){24,}GGTACC$firmask_h{$key_f}.*/){
								if ($resultados_h{$key_c.'_'.$key_p.'_'.$key_f}){
								$resultados_h{$key_c.'_'.$key_p.'_'.$key_f} += 1;
								}
								else {
								$resultados_h{$key_c.'_'.$key_p.'_'.$key_f} = 1;
								}
							last;
							}
							elsif ($seq1 =~ /^(.){24,}GGTACC$revcom_keyf.*/){
								if ($resultados_h{$key_c.'_'.$key_p.'_'.$key_f}){
								$resultados_h{$key_c.'_'.$key_p.'_'.$key_f} += 1;
								}
								else {
								$resultados_h{$key_c.'_'.$key_p.'_'.$key_f} = 1;
								}
							last;
							}
						}
					last;
					}
				}
				}				
				else {
				foreach my $key_p (keys %poolsh_h){
					my $revcom_keyp = revcom($poolsh_h{$key_p});	
					if ($seq =~ /^.{85,90}$revcom_keyp(.){0,6}$/){
						foreach my $key_f (keys %firmash_h){
							my $revcom_keyf = revcom($firmash_h{$key_f});
							if ($seq1 =~ /^AAGCTT$firmash_h{$key_f}.*/){
								if ($resultados_h{$key_c.'_'.$key_p.'_'.$key_f}){
								$resultados_h{$key_c.'_'.$key_p.'_'.$key_f} += 1;
								}
								else {
								$resultados_h{$key_c.'_'.$key_p.'_'.$key_f} = 1;
								}
							last;
							}
							elsif ($seq1 =~ /^AAGCTT$revcom_keyf.*/){
								if ($resultados_h{$key_c.'_'.$key_p.'_'.$key_f}){
								$resultados_h{$key_c.'_'.$key_p.'_'.$key_f} += 1;
								}
								else {
								$resultados_h{$key_c.'_'.$key_p.'_'.$key_f} = 1;
								}
							last;
							}
						}
					last;
					}

				}				
				}
			last;
			}		
			elsif ($seq =~ /^(.){85,90}$revcom_keyc(.){0,6}$/){
				if ($key_c =~ /.*-K$/){
				foreach my $key_p (keys %poolsk_h){
					if ($seq =~ /^(.){0,5}($poolsk_h{$key_p})((.){85,90})$/p){
						my $seq1 = $3;
						foreach my $key_f (keys %firmask_h){
							my $revcom_keyf = revcom($firmask_h{$key_f});							
							if ($seq1 =~ /^GGTACC$firmask_h{$key_f}.*/){
								if ($resultados_h{$key_c.'_'.$key_p.'_'.$key_f}){
								$resultados_h{$key_c.'_'.$key_p.'_'.$key_f} += 1;
								}
								else {
								$resultados_h{$key_c.'_'.$key_p.'_'.$key_f} = 1;
								}
							last;
							}
							elsif ($seq1 =~ /^GGTACC$revcom_keyf.*/){
								if ($resultados_h{$key_c.'_'.$key_p.'_'.$key_f}){
								$resultados_h{$key_c.'_'.$key_p.'_'.$key_f} += 1;
								}
								else {
								$resultados_h{$key_c.'_'.$key_p.'_'.$key_f} = 1;
								}
							last;
							}
						}
					last;
					}
				}
				}
				else{
				foreach my $key_p (keys %poolsh_h){
					if ($seq =~ /^(.){0,5}($poolsh_h{$key_p})((.){85,90})$/p){
						my $seq1 = $3;
						foreach my $key_f (keys %firmash_h){
							my $revcom_keyf = revcom($firmash_h{$key_f});							
							if ($seq1 =~ /^(.){24,}AAGCTT$firmash_h{$key_f}.*/){
								if ($resultados_h{$key_c.'_'.$key_p.'_'.$key_f}){
								$resultados_h{$key_c.'_'.$key_p.'_'.$key_f} += 1;
								}
								else {
								$resultados_h{$key_c.'_'.$key_p.'_'.$key_f} = 1;
								}
							last;
							}
							elsif ($seq1 =~ /^(.){24,}AAGCTT$revcom_keyf.*/){
								if ($resultados_h{$key_c.'_'.$key_p.'_'.$key_f}){
								$resultados_h{$key_c.'_'.$key_p.'_'.$key_f} += 1;
								}
								else {
								$resultados_h{$key_c.'_'.$key_p.'_'.$key_f} = 1;
								}
							last;
							}
						}
					last;
					}
				}
				}
			last;
			}
		}
	}
	};
	# if an error ocurrs...
	if ($@) {
		#remove # in the next three lines to generate error.log
		#open(ERRORLOG,">>error.log");
		#print ERRORLOG "Hubo un error, identificando la firma. Error: " . $! . "\n";
		#close ERRORLOG;
		}

}
close STMSEQ;

# write to file

open (RESULTADOS,">resultados.txt");

foreach my $key_r (keys %resultados_h){
	print RESULTADOS $key_r.",".$resultados_h{$key_r}."\n";
	}

close RESULTADOS;


sub revcom {
	# A subroutine to reverse complement DNA
	# Get the DNA to be worked on...
	my ($dna) = @_;
	# First we reverse the DNA
	$dna = reverse $dna;
	# Now translate the DNA bases
	$dna =~ tr/ACGTacgt/TGCAtgca/;
	# Return the output
	return $dna;
}
