#! /usr/bin/perl
#
use strict; 
open my $fh, '<', $ARGV[0] or die "$!\n";
while (my $line = <$fh>) {
    chomp($line);
    print reverse($line)."\n";
}
close $fh;
