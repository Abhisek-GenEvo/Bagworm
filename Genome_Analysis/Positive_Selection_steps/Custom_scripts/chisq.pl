use strict;

open(IN, $ARGV[0]);

while(chomp(my $line = <IN>))
{
        my @arr = split(/\t/, $line);
        if ($arr[1] ne '')
        {
          my $l = 2*($arr[1]-$arr[2]);

          system ("printf $arr[0]\"\t\" >> file1 ; /media/Disk1/ABM/paml4.9a/bin/chi2 $arr[3] $l | grep prob >> file1");
        }
}
