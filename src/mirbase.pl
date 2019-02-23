#!/usr/bin/perl

# Copyright (C) 2018 - 2019  Sebastian Mackowiak
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


## Downloads mirbase files from www.mirbase.org, version given as argv[0]
## when argv[1] is given then gffs are also downloaded

use strict;
use warnings;

my $cdir=`pwd`;
chomp $cdir;

my $version="14";

my $wget=`wget 2>&1`;
if($wget =~ /help/){
}else{
	die "Please install wget\n";
}

if($ARGV[0]){$version=$ARGV[0];}else{die "No version given as argv1; eg. perl mirbase.pl 21\nwill download mirbase21 files\nif you give argv1 then also all gffs will be downloaded too\n";}

my $gff=0;
if($ARGV[1]){$gff=$ARGV[1];};


if(not -d "$ENV{'HOME'}/mirbase/"){mkdir "$ENV{'HOME'}/mirbase";}

chdir "$ENV{'HOME'}/mirbase/";

$a=`wget -S --spider ftp://mirbase.org/pub/mirbase/$version/ 2>&1`;
my @A=split("\n",$a);
my $nf=0;
foreach my $e(@A){
	if($e =~ /^550/ and $e =~ /No/i){
		$nf=1;
	}
}
if(not $nf){
	print STDERR "mirbase version $version exists\nDownloading files now to $ENV{'HOME'}/mirbase/$version/\n";
	if(not -d $version){mkdir $version;}
	chdir $version;
	my @l=split("\n",$a);
	foreach my $e(@l){
		if($e =~ /(\S+.fa.gz)/){
			if(-f $1){
				print STDERR "File $1 exists, skipping\n";
				next;
			}
			`wget ftp://mirbase.org/pub/mirbase/$version/$1`;
		}

	}

	if($gff){
		if(not -d "gffs"){mkdir "gffs";}
		chdir "gffs";
		$a=`wget -S --spider ftp://mirbase.org/pub/mirbase/$version/genomes/ 2>&1`;
		my @l=split("\n",$a);
		foreach my $e(@l){

			if($e =~ /(\w+.gff3)/){

				if(-f $1){
					print STDERR "File $1 exists, skipping\n";
					next;
				}
				`wget ftp://mirbase.org/pub/mirbase/$version/genomes/$1`;
			}
		}
	}
	chdir $cdir;

}else{
	print STDERR "mirbase version $version does not exist\n";
	$a=`wget -S --spider ftp://mirbase.org/pub/mirbase/ 2>&1`;
	print STDERR "available versions are\n";

	my @l=split("\n",$a);
	my %h;
	foreach my $e(@l){
		if($e =~ /\d+\s+\d+\s+(\d+\.*\d*)/){
			$h{$1}=1;
		}
	}
	for my $k(sort {$a <=> $b} keys %h){
		print "$k\n";
	}
}	

print STDERR "If you want to extract mirnas for a specific species and you have mirdeep2 installed
then you can just run

extract_miRNAs hairpin.fa hsa > hairpin_hsa.fa

	or for more species like

extract_miRNAs hairpin.fa hsa,mmu > hairpin_hsa_mmu.fa

\n";



