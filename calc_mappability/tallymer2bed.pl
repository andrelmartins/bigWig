#!/usr/bin/perl
@chr = @ARGV;
$last = -1;
while(<STDIN>) {
	@SPL   = split(/[\t\s+-]/);
	$chri  = @SPL[0];
	$pos   = @SPL[2];

	if($last != $pos) {
		print @chr[$chri]."\t".$pos."\t".($pos+1)."\n"; 
	}

	$last = $pos;
}
