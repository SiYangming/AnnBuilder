#!/usr/bin/perl

# A perl template file for parsing the LocusLink and UniGene source
# file to map GenBank Accession numbers to LocusLink ids.

open(DATA, "</scratch/homes/jzhang/R/library/AnnBuilder/temp/file6ceaf087Rn.data") || die "Can not open";
 
open(OUT, ">/scratch/homes/jzhang/R/library/AnnBuilder/temp/gb2ll.out") || die "Can not open #OUTFILE#";


%outhash = ();

$/ = "//";

while( <DATA> ) {
    $locusid = "NA";
    $gb = "NA";
    @lines = split("\n");   

    foreach $x (@lines) {
        if( $x =~ /^LOCUSLINK/ ) {
	    #@temp = split/\s+|;/, $x;
            @temp = split/LOCUSLINK/, $x;
	    ($locusid = $temp[1]) =~ s/\s+//g;
	    next;
	}

	if( $x =~ /^SEQUENCE\s+ACC=([0-9a-zA-Z]*)[\.]?.*;\s+.*$/x  ) {
             if($locusid ne "NA"){
	         if (exists $outhash{$1}){
	             $outhash{$1} = $outhash{$1} . ";" . $locusid;
                 }else{
                     $outhash{$1} = $locusid;
                 }
             }  
	     next;
   	}
    }
}

close DATA;
open(DATA2, "< /scratch/homes/jzhang/R/library/AnnBuilder/temp/file7c83e458LL_tmpl") || die "Can not open #LLFILE#";

$/ = ">>";

while( <DATA2> ) {
    $locusid = "NA";
    $orgFound = "FALSE";

    @lines = split("\n");

    foreach $x (@lines) {
        if( $x =~ /^LOCUSID/ ) {
            @temp = split/\s+/, $x;
                $locusid = $temp[1];
                next;
        }

        if( $x =~ /^ORGANISM/){
            @temp = split/:\s+/, $x;
            if($temp[1] eq "Homo sapiens"){
                $orgFound = "TRUE";
            }
            next;
        }
        
        if($orgFound eq "TRUE"){
            if( $x =~ /^ACCNUM/ ) {
                @temp = split/\|/, $x;
                @temp = split/\s+/, $temp[0]; 
                if (exists $outHash{$temp[1]}){
                    if($temp[1] ne "none"){
	                $outHash{$temp[1]} = $outHash{$temp[1]} . ";" . 
                                                       $locusid;
                    }
                }else{
                    if($temp[1] ne "none"){
                        $outHash{$temp[1]} = $locusid;
                    }
                }
                next;
            }
        }
    }
}

close DATA2;

foreach $x (keys %outHash){
    @temp = split/;/, $outHash{$x};
    %seen = ();
    $uniq = "NA";
    foreach $item (@temp){
        if($uniq ne "NA"){
            $uniq = $uniq . ";" . $item unless($seen{$item})++;
        }else{
            $uniq = $item;
        }
    }
    print OUT $x, "\t", $uniq, "\n";
}


