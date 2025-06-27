use Data::Dumper

open(i1, $ARGV[0]);
open (o1, ">$ARGV[0]\_roottotipBL.txt");
while(<i1>)
{
    chomp; $name = $_;
    system ("cp raxml_best_tree/$_ .");
    system("Rscript branchlength.R $_");
    my %h;
    open(i2, "$name\_tmp");
    while(chomp(my $l = <i2>))
    {
        @a = split("\t", $l);
        @a1 = split("#", $a[0]);
        $h{$a1[0]}="$a[1]";
    }
    my @txt;
    foreach my $name1 (sort { $h{$b} <=> $h{$a} } keys %h)
    {
        push @txt, $name1;
    }

    if ($txt[0] eq "Eumeta1")
    {
        print o1 "$name\n";
    }

} close(i1); close(o1); close(i2)
