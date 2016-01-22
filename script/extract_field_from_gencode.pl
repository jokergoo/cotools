 
use strict;

my $file = $ARGV[0];
my $level = $ARGV[1];  # gene
my $primary_key = $ARGV[2];  # gene_id
my $field = $ARGV[3];  # gene_type

open F, $file or die "cannot open $file\n";

my $a;
my $b;
while(my $line = <F>) {
	if($line =~/^#/) {
		next;
	}

	chomp $line;
	my @tmp = split "\t", $line;
	if($tmp[2] ne $level) {
		next;
	}

	if($line =~/$primary_key "(.*?)"/) {
		$a = $1;
		if($line =~/$field "(.*?)"/) {
			$b = $1;
			print "$a\t$b\n"
		}
	}
}
