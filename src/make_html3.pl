#!/usr/bin/env perl 

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

## author SM
## date 2017-12-07

use strict;
use warnings;
use Getopt::Std;
use Cwd;

use PDF::API2;              ## needed to create a PDF file
use Math::Trig;             ## needed for sinus function
use List::Util qw(min max); ## needed to provide min and max function for lists  
use File::Basename;

my $version="1.0.0";
my $script_id='make_html3.pl';
Usage() if(not $ARGV[0] or $ARGV[0] =~ /-*-he*l*p*/);


## Variable definitions 
my $n_count = 0;
my $k_count = 0;
my $e_count = 0;
my $sig = 0;

my $infile;         
my $pdfs = 0;       ## force pdf creation
my $csv = 0;

my $xposshift = 40;
my $lstruct_multi;
my $sb_obs;

## read in available organisms at the end of the script
my %organisms;
while(<DATA>){
    chomp;
    my $tmp;
    $tmp=$_;
    $_ =~ s/\s+//g;
    $organisms{lc($_)}=$tmp;
}
my $known;

## read in quantify output file
my $id;
my %hash;  
my %seen;


my $created = 1;

my $in;
my $counter=0;

my %struct; ## counts how often a nt is covered reads

my $i;

my $offset = 0;

my $me=0;     ## mature end coordinate
my @desc;

my %mat_pre_arf =();


my $lflank1; ## length(string of left flank)
my $fl1;    ## 
my $lflank2; ## length string of right flank

my $fl2b=0;   ## right flank begin  
my $lloop;   ## string of loop
my $lb=0;     ## starting position of loop

my $lstar;   ## string of star sequence
my $sb=0;     ## starting 
my $lmature; ## string of mature 

my $mb=0;     ## mature begin
my $struct; ## structure string
my $pri_seq;## pri-cursor sequence
my $lenstr=0; 

my $pdf;    ## pdf descriptor
my $page;   ## page descriptor
my $gfx;    ## graphic variable

my $trb;    ## fontvariable
my $text;
my $text2;


my $aligned;                          ## reads reading 
my %hash2;                            ## begin of read in precursor

my %hash2c;                           ## number of reads per read sequence 
my %hash2key;
my %hash2mm;                          ## number of mismatches
my %hash2order;                       ## output order saved
my %hash2seq;
my %hash2sample;
my %star_exp_hit_pos;
my $blat;
my $spacer;   ## length of longest entry
my $spaces;   ## string of spaces to fill up spacer


my %order;                            ## stores begin coordinates of fl1,m,l,s,fl2 sequences
my $multiplier = 3.6;#4.825;               ## minimal distance between two letters

## calculate predefined pdf loci for alignment characters
my %position_hash;
$counter = 0;
for(my $i=0;$i < 200; $i++){
    $position_hash{$counter} = $xposshift+$multiplier+20+$i*$multiplier;
    $counter++;
}


my $yorig = 500; ## 500
my $downy = 50;

my $dline;                            ## line graphic handler
my $first=1;
my $lastx;
my $lasty;

my $final;                            ## final output string of a read
my @pseq;                             ## precursor sequence  
my @rseq;                             ## read sequence

my $totalreads = 0;

my %assign_str;                       ## color assigned to position where letter is drawn
my %assign_str_exp;

my $bpo1=-10;                             ## left nt pos in first bp 
my $bpo2=-10;                             ## right nt pos in first bp 
my $bpo1r=-10;                            ## left nt pos in second bp 
my $bpo2r=-10;                            ## right nt pos in second bp 


my $ffe=0;                              ## first flank end position
my $ff2b=0;                             ## second flank begin position

my @sorted;                           ## array that stores sorted order of fl1,m,l,s,fl2 
my $y=$yorig;                         ## y coordinate


my ($minx,$miny,$maxx,$maxy);        ## min and max x,y coordinates of rna sequence
my @rna;                 ## rna sequence
my @rna_d;
my %xc;                  ## holds x cooridnate of each nt
my %yc;                  ## holds y coordinate of each nt
my $sid="";        

## pdf histogram colors
my $col_star_exp = 'lightskyblue';
my $col_star_obs = 'darkviolet';
$col_star_obs='#e60000';
my $col_mature = 'red';
$col_mature="#0066ff";
my $col_loop = 'orange';

my %hm;
my %hs;
my %hp;

my $mirbase = 0;
my %mature2hairpin;
my %hairpin2mature;  ## some hairpins have more than 1 mature assigned, circumvent this problem
my %hash_q; ## takes up all entries from the quantifier module

my $blast="http://blast.ncbi.nlm.nih.gov/Blast.cgi?QUERY=";
my $blast_query = "&db=nucleotide&QUERY_FROM=&QUERY_TO=&QUERYFILE=&GENETIC_CODE=1&SUBJECTS=&stype=nucleotide&SUBJECTS_FROM=&SUBJECTS_TO=&SUBJECTFILE=&DBTYPE=gc&DATABASE=nr&EQ_MENU=&NUM_ORG=1&EQ_TEXT=&BLAST_PROGRAMS=blastn&PHI_PATTERN=&MAX_NUM_SEQ=100&SHORT_QUERY_ADJUST=on&EXPECT=10&WORD_SIZE=7&MATRIX_NAME=PAM30&MATCH_SCORES=2,-3&GAPCOSTS=5+2&COMPOSITION_BASED_STATISTICS=0&FILTER=L&REPEATS=repeat_9606&FILTER=m&TEMPLATE_LENGTH=0&TEMPLATE_TYPE=0&PSSM=&I_THRESH=&SHOW_OVERVIEW=true&SHOW_LINKOUT=true&GET_SEQUENCE=auauauaauauauauauauuauaa&FORMAT_OBJECT=Alignment&FORMAT_TYPE=HTML&ALIGNMENT_VIEW=Pairwise&MASK_CHAR=2&MASK_COLOR=1&DESCRIPTIONS=100&ALIGNMENTS=100&NEW_VIEW=true&OLD_BLAST=false&NCBI_GI=false&SHOW_CDS_FEATURE=false&NUM_OVERVIEW=100&FORMAT_EQ_TEXT=&FORMAT_ORGANISM=&EXPECT_LOW=&EXPECT_HIGH=&QUERY_INDEX=&CLIENT=web&SERVICE=plain&CMD=request&PAGE=Nucleotides&PROGRAM=blastn&MEGABLAST=&RUN_PSIBLAST=&TWO_HITS=&DEFAULT_PROG=megaBlast&WWW_BLAST_TYPE=&DB_ABBR=&SAVED_PSSM=&SELECTED_PROG_TYPE=blastn&SAVED_SEARCH=true&BLAST_SPEC=&QUERY_BELIEVE_DEFLINE=&DB_DIR_PREFIX=&USER_DATABASE=&USER_WORD_SIZE=&USER_MATCH_SCORES=&USER_FORMAT_DEFAULTS=&NO_COMMON=&NUM_DIFFS=2&NUM_OPTS_DIFFS=1&UNIQ_DEFAULTS_NAME=A_SearchDefaults_1Mn7ZD_2Sq4_1Z58HQ5Jb_23tpbD_167y9p&PAGE_TYPE=BlastSearch&USER_DEFAULT_PROG_TYPE=blastn&USER_DEFAULT_MATCH_SCORES=3.";


########################################################################
##
##
## Reading command line parameters
##
##
########################################################################


## options
my %options=();

getopts("Auf:ot:eq:dab:i:j:lm:M:PW:SzNZGp:",\%options);

## if not options{'P'} we cant use options{'Z'}
$options{'Z'} =0 if(not $options{'P'});

my $trna=0;
$trna=1 if(`head -n1 $options{'q'}|grep -i trna`); ## this we need to get from somewhere

my $snrna=0;
$snrna=1 if(`head -n1 $options{'q'}|grep -i snor`); ## this we need to get from somewhere

my %weighted=();
if($options{'W'} and -s $options{'W'}){
    open IN,$options{'W'} or die "Could not open file $options{'W'}\n";
    while(<IN>){
        if(/(\S+)\s+(\d+)/){
            $weighted{$1}=$2;
        }
    }
    close IN;
}

## everything else given to it corresponds to the samples
my @files_mirnaex=split(",",$options{'M'});
foreach(@files_mirnaex){
    print STDERR "$_ file with miRNA expression values\n";
}


my $time = "na";
if($options{'q'} =~ /expression_analyses\/expression_analyses_(.+)\/miRBase.mrd$/){
    $time=$1;	
}else{
    die "no timestamp matched at options 'q'\n";
}
my $exdir = "expression_analyses/expression_analyses_${time}";

## obtain current working directory
my $cwd = cwd;

## organism parameter
my $org='';
if($options{'t'}){$org=$organisms{$options{'t'}}};



my %confident =();
if($options{b}){
    open IN,"<$options{'b'}" or die "file not found\n";
    while(<IN>){
        next if(/\#/);
        chomp;
        my $r = $_;
        $r =~ s/\|/_/g;
        $confident{$r} = 1;
    }
    close IN;
}

############################################
##
## creating html and pdf files 
##
############################################

PrintQuantifier();

CloseHTML();

system("cp $exdir/expression_${time}.html expression_${time}.html");

print STDERR "###########################
#
#
# You can already open the html overview table 
#
# => expression_${time}.html
#
#
#########################################
";

if(not $options{'d'} or $options{'p'}){
    $mirbase = 1;
    CreateStructurePDFQuantifier(%hash_q);
}
exit;

#############################################
##
## Subroutines
##
#############################################


sub CreateStructurePDFQuantifier{
    my %hash = @_;
    my $filename;
    print STDERR "creating PDF files\n";

    if(not -d "$cwd/pdfs_$time"){
        mkdir "$cwd/pdfs_$time";
    }
    my $tk =scalar keys %hash;
    for(sort { $hash{$b}{"freq_total"} <=> $hash{$a}{"freq_total"} } keys %hash){
        next if($options{'p'} and $options{'p'} ne $_); ## we can make a single pdf given by option -p 
        next if(not $hash{$_}{'pdf'});
        next if(not $hash{$_}{'freq_total'});
        $sid = $_;
        $sid =~ tr/\|/_/;
        %star_exp_hit_pos =();
        $filename = $sid;

        if($mirbase){
            $filename = $sid; #$hairpin2mature{$sid};
        }

        next if ($seen{$filename});
        next if(-f "$cwd/pdfs_$time/$filename.pdf" and not $options{'N'});
        ## reinit variables;
        $i=0;

        $offset = 0;

        $me=0;     ## mature end coordinate
        @desc=();

        $lflank1 = 0; ## length(string of left flank)
        $fl1 = 0;    ## 
        $lflank2 = 0; ## length string of right flank
        $fl2b=-1;   ## right flank begin  
        $lloop = 0;   ## string of loop
        $lstar = 0;   ## string of star sequence
        $sb=$mat_pre_arf{$sid}{'sb'};     ## starting 
        $lmature = 0; ## string of mature 
        $mb= $mat_pre_arf{$sid}{'mb'}{'pos'};     ## mature begin
        $me= $mat_pre_arf{$sid}{'me'}{'pos'}; ## we didnt set any here 
        $lb=length($hash_q{$sid}{"pri_seq"})/2; ## we guess that it maybe somewhere around half the length of the precursor sequence
        $struct = 0; ## structure string
        $pri_seq="";## pri-cursor sequence
        $lenstr=0; 

        $pdf='';    ## pdf descriptor
        $page='';   ## page descriptor
        $gfx='';    ## graphic variable
        $trb='';    ## fontvariable

        %hash2 = ();
        %hash2c = ();
        %hash2mm = ();
        %hash2order = ();
        %order = ();

        $yorig = 500;
        $downy = 50;

        $dline='';                            ## line graphic handler

        $first=1;
        $lastx=0;
        $lasty=0;

        $final="";                         ## final output string of a read
        @pseq=();                             ## precursor sequence  
        @rseq=();                             ## read sequence

        $totalreads = 0;

        %assign_str = ();
        %assign_str_exp = ();

        %struct = ();


        $bpo1=-10;                             ## left nt pos in first bp 
        $bpo2=-10;                             ## right nt pos in first bp 
        $bpo1r=-10;                            ## left nt pos in second bp 
        $bpo2r=-10;                            ## right nt pos in second bp 


        $ffe=0;                                ## first flank end position
        $ff2b=0;                               ## second flank begin position

        @sorted=();                               ## array that stores sorted order of fl1,m,l,s,fl2 
        $y=$yorig;                             ## y coordinate


        ($minx,$miny,$maxx,$maxy)=(0,0,0,0);             ## min and max x,y coordinates of rna sequence
        @rna=();                                  ## rna sequence

        %xc = ();
        %yc = ();


#### main program;       
        $pri_seq =  $hash{$sid}{"pri_seq"};
        @rna = split(//,$pri_seq);

        chomp $pri_seq;

        my @h2m=split(",",$hairpin2mature{$sid});

        ## determine max expressed and less expressed
        my $maxc=0;
        my $maxid='';
        foreach my $h(@h2m){
            if($hash{$sid}{'mapped'}{$h} > $maxc){ $maxc=$hash{$sid}{'mapped'}{$h};$maxid=$h;}
        }
        if($maxid =~ /3p/){
            $col_mature='#e60000';
            $col_star_obs="#0066ff";

        }else{
            $col_star_obs='#e60000';
            $col_mature="#0066ff";
        }


        my @desc2 =split(//,$hash{$sid}{'exp'});
        for (my $i=0;$i < scalar @desc2; $i++){
            if ($desc2[$i] eq "f"){   ## assign_str now starts at 0 not at one
                $assign_str_exp{$i} = "black";
                $assign_str{$i} = "black";
            } elsif ($desc2[$i] eq "l"){
                $assign_str_exp{$i} = $col_loop;
                $assign_str{$i} = $col_loop;

            } elsif ($desc2[$i] =~ /[S3]+/){
                $assign_str_exp{$i} = $col_star_obs;
                $assign_str{$i} = $col_star_obs;
            } elsif ($desc2[$i] =~/[M5]+/){
                $assign_str_exp{$i} = $col_mature;
                $assign_str{$i} = $col_mature;
            } elsif($desc2[$i] eq "P"){
                $assign_str_exp{$i} = 'grey';
                $assign_str{$i} = 'grey';

            }else{
                $assign_str_exp{$i} = "black";
                $assign_str{$i} = "black";
            }
        }



        $struct = $hash{$sid}{"pri_struct"};
        $lstruct_multi = ((length($struct)+2)*$multiplier);


        my ($lv1,$lv2,$llv);	
        ## new approach
        for my $tag(keys %{$hash{$sid}{"reads"}}){ ## for each tag
            for my $read(sort keys %{$hash{$sid}{"reads"}{$tag}}){
                next if($read =~ /x/);
#                print "$tag x \t $read r \t $sid s \n";
                if ($hash{$sid}{"reads"}{$tag}{$read}{"seq"} =~ /^(\.*)(\w+)\.*$/){
                    my $v1=$1;
                    my $v2=$2;
                    $lv1=length($v1);
                    $lv2=length($v2);
                    $llv=$lv1+$lv2;
                    $hash2{$read}=$lv1; ## begin of read in precursor
                    if($hash{$sid}{"reads"}{$tag}{$read}{"rid"} =~ /_x(\d+)/){
                        my $dc = $1;
                        if($options{'W'}){
                            $dc/=$weighted{$hash{$sid}{"reads"}{$tag}{$read}{"rid"}};
                            #die "here  ===  $hash{$sid}{'reads'}{$tag}{$read}{'rid'} $dc\n";
                        }else{
                            ## check how often the read id occurs in precursor
                            #$dc/=$hash{$sid}{"reads"}{$tag}{$hash{$sid}{"reads"}{$tag}{$read}{"rid"}}; ## this works until here but is not output properly ...
                        }

                        $totalreads+= $dc;
                        #$hash2c{$tag}{$v2}+=$1;     ## number of reads with same sequence
                        $hash2c{$tag}{$v2}+=$dc;
                        $hash2key{$read}=$v2;
                        $hash2order{$read} = $read;
                        $hash2mm{$read}=$hash{$sid}{"reads"}{$tag}{$read}{"mm"}; ## number of mismatches with precursor sequence
                        $hash2seq{$read} = $hash{$sid}{"reads"}{$tag}{$read}{"seq"};


                        $hash2sample{$read} = $tag;

                        for ($i=$lv1; $i < $llv; $i++)
                        {            
                            $struct{$i}+= $dc; ## saves how often a nt in precursor is covered by a read
                        }
                    }

                }
            }
        }                           ## end of reading aligned sequences

        $y = $yorig-$downy;
        chdir "./pdfs_$time";
        CreatePDFQuantifier(\%hash);

##############################################################################################
##############################################################################################
## insert secondary structure now
        if(not $options{'A'}){	
            DrawStructure($filename);
        }
        if($totalreads ne '0'){
            CreateHistogramQuantifier();
            #ClosePDF($filename);
            #exit;
        }
        ## by ref
        CreateAlignmentQuantifier(\%hash);
        $y -=20;

### here the Frequency histogram is drawn



        ClosePDF($filename);
        chdir "..";
        $tk--;
        print STDERR "$tk left; creating pdf for $filename finished\r";
    }
    print STDERR "\n";
}

sub CreateHistogramQuantifier{
#    my ($hash) = @_;
    $dline->linewidth(2);

    $y = $yorig-$downy+90;

    $y = 560 if($options{'A'});

    $dline = $page->gfx;

    ##draw axes
    $dline->strokecolor('black');


    $dline->move($xposshift+20,$y+160);
    $dline->line($xposshift+20,$y+50-1);
    $dline->stroke;
    $dline->move($xposshift+20,$y+50);
    $dline->strokecolor('grey');
    $dline->line($xposshift+20+$lstruct_multi,$y+50);
    $dline->stroke;
    $dline->strokecolor('black');
    $dline->move($xposshift+17+$lstruct_multi,$y+53);
    $dline->line($xposshift+20+$lstruct_multi,$y+50);
    $dline->line($xposshift+17+$lstruct_multi,$y+47);

    $dline->move($xposshift+17,$y+157);
    $dline->line($xposshift+20,$y+160);
    $dline->line($xposshift+23,$y+157);

    $dline->move($xposshift+17,$y+150); 
    $dline->line($xposshift+23,$y+150);


    $gfx->textlabel($xposshift+12,$y+165,$trb,6,"freq." ,-color=>'black'); 
    $gfx->textlabel($xposshift+$lstruct_multi,$y+40,$trb,8,"length" ,-color=>'black'); 

    $gfx->textlabel($xposshift+10,$y+148,$trb,6,"1" ,-color=>'black'); 


    $dline->move($xposshift+17,$y+125); ##.75
    $dline->line($xposshift+23,$y+125);
    $gfx->textlabel($xposshift+2,$y+122,$trb,6,"0.75" ,-color=>'black'); 

    $dline->move($xposshift+17,$y+100); ## .5
    $dline->line($xposshift+23,$y+100);
    $gfx->textlabel($xposshift+6,$y+98,$trb,6,"0.5" ,-color=>'black'); 

    $dline->move($xposshift+17,$y+75); ## .25
    $dline->line($xposshift+23,$y+75); ## .25
    $gfx->textlabel($xposshift+2,$y+73,$trb,6,"0.25",-color=>'black'); 

    $gfx->textlabel($xposshift+12,$y+48,$trb,6,"0" ,-color=>'black'); 
    $dline->stroke;

    ## draw flank1
    $dline->strokecolor('black');
    $dline->move($xposshift+20,$y+50);

    $lastx = $position_hash{0};
    if(not $struct{0}){$struct{0}=0;}## that should be made in a clean way but for now it works


    $lasty = (($struct{0}/$totalreads)*100)+$y+50;

    $dline->strokecolor($assign_str{0});
    $dline->move($lastx,$lasty);
    for($i = 0; $i < length($struct); $i++){
        $dline->strokecolor($assign_str{$i});
        $dline->move($lastx,$lasty);
        $struct{$i} = 0 if(not $struct{$i}); ## that should be made in a clean way but for now it works
        $dline->line($position_hash{$i+1},(($struct{$i}/$totalreads)*100)+$y+50); ## .25
        $dline->stroke;
        $lastx = $position_hash{$i+1};
        $lasty = (($struct{$i}/$totalreads)*100)+$y+50;
    }

}


sub CreatePDFQuantifier{
    my ($hash) = @_;
    $pdf=PDF::API2->new; 

    $spacer = length($sid);
    $pdf->mediabox('A4');
    $page=$pdf->page;
    $gfx=$page->gfx;
    $text=$gfx;
    $trb=$pdf->corefont('Times-Roman', -encode=>'latin1');


    ## move everything except the structure downwards if $mirbase is set
    my $madd = 60;

    $gfx->textlabel($xposshift+20,$y+300+$downy,$trb,8,"Precursor",-color=>'black');
    $gfx->textlabel($xposshift+110,$y+300+$downy,$trb,8,": $sid",-color=>'black');


    #$spaces = " " x ($spacer - length($$hash{$sid}{"freq_total"}));
    $gfx->textlabel($xposshift+20,$y+230+$madd+$downy,$trb,8,"Total read count",-color=>'black');       
    $gfx->textlabel($xposshift+110,$y+230+$madd+$downy,$trb,8,": $$hash{$sid}{'freq_total'}",-color=>'black');

    ## here should be written how many annotated stuff is actually there and how many not
    my $jk =10;
    # old
    #for(sort {$$hash{$sid}{'mapped'}{$b} <=> $$hash{$sid}{'mapped'}{$a}} keys %{$$hash{$sid}{'mapped'}}){
    #new
    my @h2m=split(",",$hairpin2mature{$sid});

    my $xpadd=100;

    ## determine max expressed and less expressed
    my $maxc=0;
    my $maxid='';
    foreach my $h(@h2m){
        if($$hash{$sid}{'mapped'}{$h} > $maxc){ $maxc=$$hash{$sid}{'mapped'}{$h};$maxid=$h;}
#		print "$h\t$$hash{$sid}{'mapped'}{$h}\n$maxc\t$maxid\n";
    }

    foreach my $h(@h2m){
        next if($h =~ /^\s*$/);
        $xpadd=110;	
        if($options{'t'}){
            next if($_ !~ $options{'m'});
        }
        my $idlen=length($_);
        $spaces = " " x ($spacer - $idlen);
        $gfx->textlabel($xposshift+20,$y+230-$jk+$madd+$downy,$trb,8,"$h counts",-color=>'black');
        if($idlen > 16){$xpadd+= (  ($idlen-14)*5);}
        $gfx->textlabel($xposshift+$xpadd,$y+230-$jk+$madd+$downy,$trb,8,": $$hash{$sid}{'mapped'}{$h}",-color=>'black');
        $jk+=10;
    }
    #$spaces = " " x ($spacer - length("remaining reads"));
    $gfx->textlabel($xposshift+20,$y+230-$jk+$madd+$downy,$trb,8,"remaining reads",-color=>'black');      
    $gfx->textlabel($xposshift+110,$y+230-$jk+$madd+$downy,$trb,8,": $$hash{$sid}{'remaining_rc'}",-color=>'black');
    $jk+=10;
    $trb=$pdf->corefont('Courier', -encode=>'latin1');
    $dline = $page->gfx;
}




sub ClosePDF{
    my $file = shift;
    $file = "output" if($file eq"");
    $pdf->saveas("$cwd/pdfs_$time/$file.pdf");
}



## draw a line in PDF 
sub Line{
    my ($x1,$y1,$x2,$y2,$col,$width) = @_;
    #return 0;
    #$dline->linewidth($width);
    #$dline->strokecolor($col);
    $dline->move($x1,$y1);
    $dline->line($x2,$y2);
    $dline->stroke;
}

## draw a letter in PDF
sub Base{
    my ($x1,$y1,$base,$col,$size,$trb) = @_;
    #return 0;
    $gfx->textlabel($x1,$y1,$trb,$size,$base,-color=>$col);
}	

sub CreateAlignmentQuantifier{
    my ($hash) = @_;
    $y+=20;

    if($trna and $options{'i'} !~ /dummy/){

        my @k= sort keys %{$mat_pre_arf{$sid}};
        $gfx->textlabel($xposshift+(20+(($mat_pre_arf{$sid}{$k[0]}{'pos'}-3)*$multiplier)),$y,$trb,6,"Anticodon",-color=>'orange');
        $y-=10;
    }

    for my $k1(sort keys %{$mat_pre_arf{$sid}}){
        if($k1 =~ /\*/){
            $gfx->textlabel($xposshift+(20+(($mat_pre_arf{$sid}{$k1}{'pos'})*$multiplier)),$y,$trb,6,"$k1",-color=>$col_star_obs);
        }else{
            if($k1 =~ /trna\S+3p/){
                $gfx->textlabel($xposshift+(20+(($mat_pre_arf{$sid}{$k1}{'pos'})*$multiplier)),$y,$trb,6,"$k1",-color=>$col_star_obs);
            }else{
                ### here
                next if($k1 =~ /mb/ or $k1 =~ /me/ or $k1 =~ /sb/ or $k1 =~ /se/);   
                $gfx->textlabel($xposshift+(20+(($mat_pre_arf{$sid}{$k1}{'pos'})*$multiplier)),$y,$trb,6,"$k1",-color=>$assign_str{$mat_pre_arf{$sid}{$k1}{'pos'}});
            }
        }
        $y-=10;
    }

    for(my $i=0; $i < scalar @rna; $i++){
        $gfx->textlabel($position_hash{$i},$y,$trb,6,$rna[$i],-color=>$assign_str{$i});   
    }

    $gfx->textlabel($xposshift+25+ $lstruct_multi,$y,$trb,6,'-3\'' ,-color=>'black');
    $gfx->textlabel($xposshift+10,$y,$trb,6,'5\'-' ,-color=>'black');

    if($$hash{$sid}{'obs'}){
        $gfx->textlabel($xposshift+50+ $lstruct_multi,$y,$trb,6,'obs' ,-color=>'black');
    }else{
        $gfx->textlabel($xposshift+50+ $lstruct_multi,$y,$trb,6,'exp' ,-color=>'black');
    }

    if($$hash{$sid}{'obs'}){
        $y -= 10;
        for(my $i=0; $i < scalar @rna; $i++){
            $gfx->textlabel($position_hash{$i},$y,$trb,6,$rna[$i],-color=>$assign_str_exp{$i}); 
        }
        $gfx->textlabel($xposshift+50+ $lstruct_multi,$y,$trb,6,'exp' ,-color=>'black');

    }

    return if($options{'G'});


    $y -= 10;

    my @structx = split(//,$struct);
    my $sadd = 0;
    $gfx->textlabel($position_hash{0},$y,$trb,6,$struct,-color=>'black');
    $gfx->textlabel($xposshift+30+ $lstruct_multi,$y,$trb,6,'reads' ,-color=>'black');
    $gfx->textlabel($xposshift+70+ $lstruct_multi,$y,$trb,6,'mm' ,-color=>'black');
    $gfx->textlabel($xposshift+110+ $lstruct_multi,$y,$trb,6,'sample' ,-color=>'black');


    $y -= 10;    
    if($options{'o'}){
        for my $tag(keys %{$$hash{$sid}{"reads"}}){
            for my $k(sort { $hash2order{$a} <=> $hash2order{$b} } keys %hash2order){
                next if($hash2sample{$k} ne $tag);

                $gfx->textlabel($position_hash{0},$y,$trb,6,$hash2seq{$k},-color=>'black');

                ## matches and read numbers
                $gfx->textlabel($xposshift+30+ $lstruct_multi,$y,$trb,6,int($hash2c{$tag}{$hash2key{$k}}) ,-color=>'black');
                $gfx->textlabel($xposshift+70+ $lstruct_multi,$y,$trb,6,$hash2mm{$k} ,-color=>'black');
                $gfx->textlabel($xposshift+110+ $lstruct_multi,$y,$trb,6,$hash2sample{$k} ,-color=>'black');

                $y -= 10;
                if($y < 100){
                    $page=$pdf->page();
                    $pdf->mediabox('A4');
                    $text2=$page->text();
                    $gfx = $text2;
                    $y=800;
                    for(my $i=0; $i < scalar @rna; $i++){
                        $gfx->textlabel($position_hash{$i},$y,$trb,6,$rna[$i],-color=>$assign_str{$i});
                    }
                    $y+=20;
                    for my $k1(keys %{$mat_pre_arf{$sid}}){
                        next if($k1 =~ /mb/ or $k1 =~ /me/ or $k1 =~ /sb/ or $k1 =~ /se/);   
                        if($k1 =~ /\*/){
                            $gfx->textlabel($xposshift+(20+(($mat_pre_arf{$sid}{$k1}{'pos'})*$multiplier)),$y,$trb,6,"$k1",-color=>$col_star_obs);
                        }else{

                            next if($k1 =~ /mb/ or $k1 =~ /me/ or $k1 =~ /sb/ or $k1 =~ /se/);   
                            $gfx->textlabel($xposshift+(20+(($mat_pre_arf{$sid}{$k1}{'pos'})*$multiplier)),$y,$trb,6,"$k1",-color=>$assign_str{$mat_pre_arf{$sid}{$k1}{'pos'}});
                            #$gfx->textlabel($xposshift+(20+(($mat_pre_arf{$sid}{$k1}{'pos'})*$multiplier)),$y,$trb,6,"$k1",-color=>$col_mature);
                        }
                        $y-=10;
                    }
                    $y-=20;
                }
            }
            $y-=10;
        }
    }else{
        for my $k(sort { $hash2order{$a} <=> $hash2order{$b} } keys %hash2order){
            my $tag = $hash2sample{$k};
            $gfx->textlabel($position_hash{0},$y,$trb,6,$hash2seq{$k},-color=>'black');

            ## matches and read numbers
            $gfx->textlabel($xposshift+30+ $lstruct_multi,$y,$trb,6,int($hash2c{$tag}{$hash2key{$k}}) ,-color=>'black');
            $gfx->textlabel($xposshift+70+ $lstruct_multi,$y,$trb,6,$hash2mm{$k} ,-color=>'black');
            $gfx->textlabel($xposshift+110+ $lstruct_multi,$y,$trb,6,$hash2sample{$k} ,-color=>'black');

            $y -= 10;
            if($y < 100){
                $page=$pdf->page();
                $pdf->mediabox('A4');
                $text2=$page->text();
                $gfx = $text2;
                $y=800;
                $y+=20;
                for my $k1(keys %{$mat_pre_arf{$sid}}){
                    next if($k1 =~ /mb/ or $k1 =~ /me/ or $k1 =~ /sb/ or $k1 =~ /se/);   
                    if($k1 =~ /\*/){
                        $gfx->textlabel($xposshift+(20+(($mat_pre_arf{$sid}{$k1}{'pos'})*$multiplier)),$y,$trb,6,"$k1",-color=>$col_star_obs);
                    }else{
                        next if($k1 =~ /mb/ or $k1 =~ /me/ or $k1 =~ /sb/ or $k1 =~ /se/);   
                        $gfx->textlabel($xposshift+(20+(($mat_pre_arf{$sid}{$k1}{'pos'})*$multiplier)),$y,$trb,6,"$k1",-color=>$assign_str{$mat_pre_arf{$sid}{$k1}{'pos'}});
                        #$gfx->textlabel($xposshift+(20+(($mat_pre_arf{$sid}{$k1}{'pos'})*$multiplier)),$y,$trb,6,"$k1",-color=>$col_mature);
                    }
                    $y-=10;
                }
                for(my $i=0; $i < scalar @rna; $i++){
                    $gfx->textlabel($position_hash{$i},$y,$trb,6,$rna[$i],-color=>$assign_str{$i});
                }
                $y-=10;
            }
        }
    }
}


sub Shifting{
    $minx = min(values %xc);
    $maxx = max(values %xc);
    $miny = min(values %yc);
    $maxy = max(values %yc);


    ## now place sec-structure in upper right corner
    my $shiftx = abs($minx)+10;
    my $shifty = abs($miny)+10;

    #shift everything to printable area;
    if($minx < 0){
        for(my $i=0; $i < scalar @rna_d; $i++){
            $xc{$i}+=$shiftx;
        }
    }
    if($miny < 0){
        for(my $i=0; $i < scalar @rna_d; $i++){
            $yc{$i}+=$shifty;
        }
    }

    if($maxx > 600){
        for(my $i=0; $i < scalar @rna_d; $i++){
            $xc{$i}-=$minx+10;
        }

    }
    if($maxy > 800){
        for(my $i=0; $i < scalar @rna_d; $i++){
            $yc{$i}-=$miny+10;
        }
    }	
}


sub DrawStructure{
    my $filename = shift;
    #return 
    ## run RNAplot to create secondary structure ps file
    my $cw =cwd;
    #print STDERR "$sid\tcwd\t$cw\n";

    #exit;
    my $in_pos=0;
    my $in_pairs=0;

    my $count=0;                ## counter


    my %bp;                     ## which nt pair

    my $line;                   ## input line of file
    my $fh;
    if(! -f "$ENV{'HOME'}/mirbase/21/cache/$filename.ps"){ 
        system("echo 5${pri_seq}3| RNAfold -d 0 > /dev/null;mv rna.ps $ENV{'HOME'}/mirbase/21/cache/$filename.ps");
    }
    if(-f "$ENV{'HOME'}/mirbase/21/cache/$filename.ps"){}else{die "File $ENV{'HOME'}/mirbase/21/cache/$filename.ps not found here\n";}
    open $fh,"$ENV{'HOME'}/mirbase/21/cache/$filename.ps" or die "File $ENV{'HOME'}/mirbase/21/cache/$filename.ps not found in cache\n";



    my ($minx, $miny) = 10000;
    my ($maxx,$maxy)  = 0;

    my $centering_x= 0; ## base in precursor sequences that is choosen as center point for rotation
    my $twisted=0; ## if mature sequence comes after loop twisted is 1
    my $sums=0;                 ## diff between 2bp nts

    while (<$fh>)
    {
        if (/\/sequence/)       ## if sequence matched in rna.ps
        {
            $line = <$fh>;       ## read in rna sequence
            chomp $line;
            #$line =~ s/U/T/g;
            if($line =~ /\\$/){
                chop $line;
            }
            @rna_d=split(//,$line); ## read in to @rna_d
            next;
        }

        if (/\/coor\s*\[/) ## if nt coordinates section in ps comes now
        {
            $in_pos = 1;
            next;
        }
        if ($in_pos and /\[(\S+)\s+(\S+)\]/) 
        {
            $xc{$count} = $1;   ## x cooridnate
            $yc{$count} = $2;   ## y coordinate
            $count++;
            next;
        }

        if (/\/pairs/)          ## read in base pairs
        {
            $count=0;
            $in_pos=0;
            $in_pairs = 1;
            next;
        }

        if(not defined $mb){die $filename;}
        $twisted = 1 if($mb > $lb ); ## mature begin is after loop begin
#		die "$mb > $lb \n";


        if ($in_pairs and /\[(\S+)\s+(\S+)\]/){
            if ($twisted){
                $bpo2r=$1;
                $bpo1r=$2;
            } else {
                $bpo2r=$2;
                $bpo1r=$1;
            }


            ## mm difficult

            ## determine two subsequent bases in mature having subsequent paired bases
            if ($bpo1r >= $mb-$offset and $bpo1r < $me-$offset and $centering_x==0){
                if ($twisted){
                    $sums = -1;
                } else {
                    $sums=1;
                }
                if (($bpo1r-$bpo1) == $sums and ($bpo2-$bpo2r ) == $sums) {
                    if($twisted){
                        $centering_x=$bpo1r-2;
                    }else{
                        $centering_x=$bpo1r;
                    }
                }
            }
#			print "$bpo1r >= $mb - $offset and $bpo1r < $me - $offset and $centering_x ==0\n";

            $bpo1= $bpo1r;
            $bpo2 = $bpo2r;


            $bp{$bpo1r-1}=$bpo2r-1; ## saving nt pairs in hash %bp
            next;
        }
        if ($in_pairs and /\]\s*def/) ## end of nt pairs in ps file 
        {
            $in_pairs = 0;
            last;
        }
    }
    close $fh;

    Shifting();

    $minx = min(values %xc);
    $maxx = max(values %xc);
    $miny = min(values %yc);
    $maxy = max(values %yc);

    ##determine if mirror or not
    my $mir=0;

    $mir = 1 if($twisted);

    my $yshift=0;
    ########### mirror sequence so that loop is on the right hand side

    my $cshift = 3; ## determines the nt in mature sequence which should be the rotating center
    if ($mir){
        for (my $i=0; $i < scalar @rna_d; $i++)
        {
            $xc{$i} = -$xc{$i};
        }
    }
    $minx = min(values %xc);
    $maxx = max(values %xc);
    $miny = min(values %yc);
    $maxy = max(values %yc);


    #translate back to positive area
    if ($mir)
    {
        for (my $i=0; $i < scalar @rna_d; $i++)
        {
            $xc{$i} += abs($minx)+10;
        }
    }

    $minx = min(values %xc);
    $maxx = max(values %xc);
    $miny = min(values %yc);
    $maxy = max(values %yc);


    my $ax = $xc{$centering_x};
    my $ay = $yc{$centering_x};


    ## point relativ to center
    my $bx;
    my $by;
    #print "twisted $twisted\n";
    if ($twisted)
    {
        $bx = $xc{$centering_x-1};
        $by = $yc{$centering_x-1};
    } else {
        $bx = $xc{$centering_x+1};
        $by = $yc{$centering_x+1};
    }
    if(not $by){ die " here line 1015 $filename $centering_x $twisted\n";}

    my $gk = $by-$ay;
    my $ak = $bx-$ax;


    my $r = sqrt(($ak**2)+($gk**2));       

    my $phi = asin($gk/$r);

    if ($bx < $ax and $by > $ay)
    {
        $phi = 3.141593-$phi;
    }
    if ($bx <= $ax and $by <= $ay)
    {
        $phi *= (-1);
        $phi += 3.141593;
    }

    my $alpha;
    my $do_rot = 1;
    if ($do_rot)
    {
        my $last = $xc{0};
        ### now rotate every point in a designated angle of phi
        for (my $i=0; $i < scalar @rna_d; $i++)
        {
            next if ($i == $centering_x);

            $bx = $xc{$i};
            $by = $yc{$i};

            $gk = $by-$ay;
            $ak = $bx-$ax;             
            $r = sqrt(($ak**2)+($gk**2));

            $alpha = asin($gk/$r);

            if ($bx < $ax and $by > $ay)
            {
                $alpha = 3.141593-$alpha;
            }
            if ($bx <= $ax and $by <= $ay)
            {
                $alpha *= (-1);
                $alpha += 3.141593;
            }

            $alpha -= $phi;

            $xc{$i} = $ax + $r*cos($alpha);
            $yc{$i} = $ay + $r*sin($alpha);

            my $dif =  ($xc{$i} - $last);
            $last=$xc{$i};
        }
    }


    my $reduce = 0;
    my $red_dist= abs($xc{$mb+$cshift}-$xc{$mb-1+$cshift});

    Shifting();

    ## check if to mirror horizontally again because RNAfold does not take care of sequence input direction
    if (not $twisted)
    {
        my @bpkeys = sort keys %bp;
        $maxy = max(values %yc);
        if ($yc{$bpkeys[0]} < $yc{$bp{$bpkeys[0]}})
        {
            for (my $i=0; $i < scalar @rna_d; $i++)
            {
                $yc{$i} = -$yc{$i}+$maxy+10;
            }
        }
    }



    ## draw structure
    $minx = min(values %xc);
    $maxx = max(values %xc);
    $miny = min(values %yc);
    $maxy = max(values %yc);

    $y = $yorig+300;
    my $x = 550;

    my $rx = 300;

    if ($maxx < $x) {
        $x = $x-$maxx; 
    } else {
        $x = -abs($x-$maxx); 
    }

    if ($maxy < $y) {
        $y = $y-$maxy; 
    } else {
        $y = -abs($y-$maxy); 
    }

    #   print "$minx+$x\t$maxx+$x\n";

    ## scaling now structure again

    my ($scx,$scy) = (250,250); ## scaling center coordinates on page
    my ($tx,$ty);               ## translated x and y coordinates
    my $scfactor = ($rx/($maxx-$minx));  ## scaling factor
    #   print $scfactor,"\n";

    if($scfactor < 1){

        for (my $i=0; $i < scalar @rna_d; $i++)
        {
            $tx = $xc{$i}-$scx;
            $ty = $yc{$i}-$scy;
            $xc{$i} = $tx*$scfactor + $scx;
            $yc{$i} = $ty*$scfactor + $scx;
        }
    }

    $minx = min(values %xc);
    $maxx = max(values %xc);
    $miny = min(values %yc);
    $maxy = max(values %yc);

    $y = $yorig+280;#300;
    $x = 550;

    if ($maxx < $x) {
        $x = $x-$maxx; 
    } else {
        $x = -abs($x-$maxx); 
    }

    if ($maxy < $y) {
        $y = $y-$maxy; 
    }else {
        $y = -abs($y-$maxy); 
    }


    if(($minx+$x) < ($x-$rx)){ ## left most x-coord for a base 
        for(keys %xc){
            $xc{$_} += ($x-$rx)-($minx+$x);
        }
    }


    $trb=$pdf->corefont('Courierbold', -encode=>'latin1');
    $dline->linewidth(0.5);
    $dline->strokecolor("grey");


    for (my $i=1; $i < scalar @rna_d; $i++) {
        Line($xc{$i-1}+$x+2.5,$yc{$i-1}+2+$y,$xc{$i}+2.5+$x,$yc{$i}+2+$y,"grey",0.5);
    }

    Base($xc{0}+$x,$yc{0}+$y,"$rna_d[0]'",'black',8,$trb);
    for (my $i=1; $i < (scalar @rna_d)-1; $i++) {

        if($star_exp_hit_pos{$i-1+$offset}){
            Base($xc{$i}+$x,$yc{$i}+$y,$rna_d[$i],$col_star_exp,8,$trb);
        }else{
            Base($xc{$i}+$x,$yc{$i}+$y,$rna_d[$i],$assign_str{$i-1+$offset},8,$trb);
        }
    }
    Base($xc{(scalar @rna_d)-1}+$x,$yc{(scalar @rna_d)-1}+$y,"$rna_d[(scalar @rna_d)-1]'",'black',8,$trb);
    #die $rna_d[96];

    ## drawing bp lines
    my $scfactorl= 0.4;         #scaling factor

    my $fx;                     ## from x coordinate
    my $tox;                    ## from y coordinate

    my $fy;                     ## from y cooridinate
    my $toy;                    ## to y coordinate

    my $dx;                     ## xlength
    my $dy;                     ## y length

    my $dx1;   ## difference between orig x length and scaled x length
    my $dy1;   ## difference between orig y length and scaled y length


    for (keys %bp)
    {
        $dx = abs($xc{$_} - $xc{$bp{$_}});
        $dy = abs($yc{$_} - $yc{$bp{$_}});

        $dx1 = ($dx-$scfactorl*$dx)/2;
        $dy1 = ($dy-$scfactorl*$dy)/2;


        if ($xc{$_} > $xc{$bp{$_}})
        {
            $fx = $xc{$_} - $dx1;
            $tox = $xc{$bp{$_}} + $dx1;
        } else {
            $fx = $xc{$_} + $dx1;
            $tox = $xc{$bp{$_}} - $dx1;

        }

        if ($yc{$_} > $yc{$bp{$_}})
        {
            $fy = $yc{$_} - $dy1;
            $toy = $yc{$bp{$_}} + $dy1;
        } else {
            $fy = $yc{$_} + $dy1;
            $toy = $yc{$bp{$_}} - $dy1;
        }
        Line($fx+2.5+$x,$fy+2+$y,$tox+2.5+$x,$toy+2+$y,"black",0.5);	
    }
}

sub CreateHTML{
## print html
    my $img=insert_image();
    print HTML <<EOF;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
        "http://www.w3.org/TR/html4/strict.dtd">
        <html>
        <head>
        <title>Quantifier</title>
        <!-- CSS code -->
            <link rel="stylesheet" type="text/css" href="mystyle.css">
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.2.1/jquery.min.js"></script>
    <script>
    \$(document).ready(function(){
        \$('#evs').hide();
        \$('#scor').hide();
        \$('#stats').hide();

    \$("#b1").click(function(){
        \$('#scor').hide();
        \$('#stats').hide();
        \$("#evs").toggle(function(){
        \$(this).load('sendplot_mod.html');
        });
    });

    \$("#b2").on("click", function(){
        \$('#evs').hide();
        \$('#stats').hide();
        \$("#scor").toggle();
    });

    \$("#b3").click(function(){
        \$('#scor').hide();
        \$('#evs').hide();
        \$('#stats').toggle(function(){
            \$(this).load('../cell_processing_statistics.html');
        });
    });

    \$("#b4").click(function(){
        \$('.zero').toggle();
    });

    });
    </script>

</head>
<body>
<table border="0" width="100%">
<colgroup>
<col width="5*">
<col width="5*">
</colgroup>
<tr height="200" valign="top" >
<td>
<!-- <img src="./quantifier.png" style="border-style: none" name="quant" title="Quantifier" align=left alt="Quantifier"/></a></td> -->
$img
</td>

<br>
<br>
</font></td>
</tr>
</table>
EOF
if(not $trna){
    print HTML <<EOF;
<button id="b3">preprocessing statistics</button>
<button id="b1">Expresssion vs noise</button>
<button id="b2">miRNA correlation</button>
<div id="evs">
</div>

<div id="scor">
        <img src = "Rplots.png"</img>
</div>
<div id="stats">
</div>
EOF
}
print HTML "
<button id=\"b4\">0entries</button>
";

}

sub CloseHTML{
    print HTML <<EOF;
</table>
</body>
</html>
EOF
}


sub PrintHtmlTableHeader{
    my ($hl,$csv) = @_;
    my %h;	

    ## divide string by linebreaks every x characters
    my $p1 ='<th><a href="" class="tooltip3">';
    my $p11='<th><a href="" class="tooltip4">';

    my $p2='<span>';
    my $q ='</span></a></th>';

    if($hl =~ /novel/i){
        $h{1}{1} = 'provisional id';
        $h{1}{2} = 'this is a provisional miRNA name assigned by miRDeep2. The first part of the id designates the chromosome or genome contig on which the miRNA gene is located. The second part is a running number that is added to avoid identical ids. The running number is incremented by one for each potential miRNA precursor that is excised from the genome. Clicking this field will display a pdf of the structure, read signature and score breakdown of the reported miRNA.';
        $h{2}{1} = 'miRDeep2 score';
        $h{2}{2} = 'the log-odds score assigned to the hairpin by miRDeep2';
        $h{3}{1} = 'estimated probability that the miRNA candidate is a true positive';
        $h{3}{2} = 'the estimated probability that a predicted novel miRNA with a score of this or higher is a true positive. To see exactly how this probability is estimated, mouse over the \'novel miRNAs, true positives\' in the table at the top of the webpage.';
        $h{4}{1} = 'rfam alert';
        $h{4}{2} = 'this field indicates if the predicted miRNA hairpin has sequence similarity to reference rRNAs or tRNAs. Warnings in this field should overrule the estimated probability that a reported miRNA is a true positive (previous field).';
        $h{5}{1} = 'total read count';
        $h{5}{2} = 'this is the sum of read counts for the predicted mature, loop and star miRNAs.';
        $h{6}{1} = 'mature read count';
        $h{6}{2} = 'this is the number of reads that map to the predicted miRNA hairpin and are contained in the sequence covered by the predicted mature miRNA, including 2 nts upstream and 5 nts downstream.';
        $h{7}{1} = 'loop read count';
        $h{7}{2} = 'this is the number of reads that map to the predicted miRNA hairpin and are contained in the sequence covered by the predicted miRNA loop, including 2 nts upstream and 5 nts downstream.';
        $h{8}{1} = 'star read count';
        $h{8}{2} = 'this is the number of reads that map to the predicted miRNA hairpin and are contained in the sequence covered by the predicted star miRNA, including 2 nts upstream and 5 nts downstream.';
        $h{9}{1} = 'significant randfold p-value';
        $h{9}{2} = 'this field indicates if the estimated randfold p-value of the excised potential miRNA hairpin is equal to or lower than 0.05 (see Bonnet et al., Bioinformatics, 2004).';
        $h{10}{1} = 'miRBase miRNA';
        $h{10}{2} = 'this field displays the ids of any reference mature miRNAs for the species that map perfectly (full length, no mismatches) to the reported miRNA hairpin. If this is the case, the reported miRNA hairpin is assigned as a known miRNA. If not, it is assigned as a novel miRNA. If more than one reference mature miRNA map to the miRNA hairpin, then only the id of the miRNA that occurs last in the input file of reference mature miRNAs for the species is displayed.';
        $h{11}{1} = 'example miRBase miRNA with the same seed';
        $h{11}{2} = 'this field displays the ids of any reference mature miRNAs from related species that have a seed sequence identical to that of the reported mature miRNA. The seed is here defined as nucleotides 2-8 from the 5\' end of the mature miRNA. If more than one reference mature miRNA have identical seed, then only the id of the miRNA that occurs last in the input file of reference mature miRNAs from related species is displayed.';
        $h{12}{1} = 'UCSC browser';
        $h{12}{2} = 'if a species name was input to miRDeep2, then clicking this field will initiate a UCSC blat search of the consensus precursor sequence against the reference genome.';
        $h{13}{1} = 'NCBI blastn';
        $h{13}{2} = 'clicking this field will initiate a NCBI blastn search of the consensus precursor sequence against the nr/nt database (non-redundant collection of all NCBI nucleotide sequences).';
        $h{14}{1} = 'consensus mature sequence';
        $h{14}{2} = 'this is the consensus mature miRNA sequence as inferred from the deep sequencing reads.';
        $h{15}{1} = 'consensus star sequence';
        $h{15}{2} = 'this is the consensus star miRNA sequence as inferred from the deep sequencing reads.';
        $h{16}{1} = 'consensus precursor sequence';
        $h{16}{2} = 'this is the consensus precursor miRNA sequence as inferred from the deep sequencing reads. Note that this is the inferred Drosha hairpin product, and therefore does not include substantial flanking genomic sequence as does most miRBase precursors.';
        if($options{'a'}){
            $h{17}{1} = 'genomic position';
        }


    }elsif($hl =~ /miRBase miRNAs in dataset/i){
        $h{1}{1} = 'tag id';
        $h{1}{2} = 'this is a tag id assigned by miRDeep2. The first part of the id designates the chromosome or genome contig on which the miRNA gene is located. The second part is a running number that is added to avoid identical ids. The running number is incremented by one for each potential miRNA precursor that is excised from the genome. Clicking this field will display a pdf of the structure, read signature and score breakdown of the miRNA.';
        $h{2}{1} = 'miRDeep2 score';
        $h{2}{2} = 'the log-odds score assigned to the hairpin by miRDeep2';
        $h{3}{1} = 'estimated probability that the miRNA is a true positive';
        $h{3}{2} = 'the estimated probability that a predicted miRNA with a score of this or higher is a true positive. To see exactly how this probability is estimated, mouse over the \'novel miRNAs, true positives\' in the table at the top of the webpage. For miRBase miRNAs, this reflects the support that the data at hand lends to the miRNA.';
        $h{4}{1} = 'rfam alert';
        $h{4}{2} = 'this field indicates if the miRNA hairpin has sequence similarity to reference rRNAs or tRNAs. Warnings in this field should overrule the estimated probability that a reported miRNA is a true positive (previous field).';
        $h{5}{1} = 'total read count';
        $h{5}{2} = 'this is the sum of read counts for the mature, loop and star miRNAs.';
        $h{6}{1} = 'mature read count';
        $h{6}{2} = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the consensus mature miRNA, including 2 nts upstream and 5 nts downstream.';
        $h{7}{1} = 'loop read count';
        $h{7}{2} = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the consensus miRNA loop, including 2 nts upstream and 5 nts downstream.';
        $h{8}{1} = 'star read count';
        $h{8}{2} = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the consensus star miRNA, including 2 nts upstream and 5 nts downstream.';
        $h{9}{1} = 'significant randfold p-value';
        $h{9}{2} = 'this field indicates if the estimated randfold p-value of the miRNA hairpin is equal to or lower than 0.05 (see Bonnet et al., Bioinformatics, 2004).';
        $h{10}{1} = 'miRBase miRNA';
        $h{10}{2} = 'this field displays the ids of any reference mature miRNAs for the species that map perfectly (full length, no mismatches) to the reported miRNA hairpin. If this is the case, the reported miRNA hairpin is assigned as a known miRNA. If not, it is assigned as a novel miRNA. If more than one reference mature miRNA map to the miRNA hairpin, then only the id of the miRNA that occurs last in the input file of reference mature miRNAs for the species is displayed. Clicking this field will link to miRBase.';
        $h{11}{1} = 'example miRBase miRNA with the same seed';
        $h{11}{2} = 'this field displays the ids of any reference mature miRNAs from related species that have a seed sequence identical to that of the reported mature miRNA. The seed is here defined as nucleotides 2-8 from the 5\' end of the mature miRNA. If more than one reference mature miRNA have identical seed, then only the id of the miRNA that occurs last in the input file of reference mature miRNAs from related species is displayed.';
        $h{12}{1} = 'UCSC browser';
        $h{12}{2} = 'if a species name was input to miRDeep2, then clicking this field will initiate a UCSC blat search of the consensus precursor sequence against the reference genome.';
        $h{13}{1} = 'NCBI blastn';
        $h{13}{2} = 'clicking this field will initiate a NCBI blastn search of the consensus precursor sequence against the nr/nt database (non-redundant collection of all NCBI nucleotide sequences).';
        $h{14}{1} = 'consensus mature sequence';
        $h{14}{2} = 'this is the consensus mature miRNA sequence as inferred from the deep sequencing reads.';
        $h{15}{1} = 'consensus star sequence';
        $h{15}{2} = 'this is the consensus star miRNA sequence as inferred from the deep sequencing reads.';
        $h{16}{1} = 'consensus precursor sequence';
        $h{16}{2} = 'this is the consensus precursor miRNA sequence as inferred from the deep sequencing reads. Note that this is the inferred Drosha hairpin product, and therefore does not include substantial flanking genomic sequence as does most miRBase precursors.';
        if($options{'a'}){
            $h{17}{1} = 'genomic position';
        }

    }else{
        if($trna){
            $h{1}{1} = 'tRNA id';
            $h{1}{2} = 'Clicking this field will display a pdf of the structure and read signature of the tRNA.';
        }elsif($snrna){
            $h{1}{1} = 'snoRNA id';
            $h{1}{2} = 'Clicking this field will display a pdf of the structure and read signature of the snoRNA.';
        }else{
            $h{1}{1} = 'miRBase precursor id';
            $h{1}{2} = 'Clicking this field will display a pdf of the structure and read signature of the miRNA.';
        }
        $h{2}{1} = '-';
        $h{2}{1} .="&#160" x 5;
        $h{2}{2} = '-';
        $h{3}{1} = '-';
        $h{3}{1} .="&#160" x 5;
        $h{3}{2} = '-';
        $h{4}{1} = '-';
        $h{4}{1} .="&#160" x 5;
        $h{4}{2} = '-';
        $h{5}{1} = 'total read count';

        if($options{'z'}){	
            $h{6}{1} = 'read counts per sample';
        }else{
            if($options{'P'}){
                $h{5}{2} = 'this is the sum of read counts for the 5p and 3p sequences.';
                $h{6}{1} = '5p read counts';
                $h{6}{2} = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the 5p sequence, including 2 nts upstream and 5 nts downstream. In parenthesis are normalized read counts shown.';

                if($options{'Z'}){
                    $h{6}{1} = '5p/3p read counts';
                    $h{6}{2} = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the 5p/3p sequence, including 2 nts upstream and 5 nts downstream. In parenthesis are normalized read counts shown.';

                }

                $h{8}{1} = '3p read counts';
                $h{8}{2} = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the 3p sequence, including 2 nts upstream and 5 nts downstream. In parenthesis are normalized read counts shown.';
                $h{9}{1} = 'remaining reads';
                $h{9}{2} = 'this is the number of reads that did not map to any of the 5p or 3p sequences';
            }else{
                $h{5}{2} = 'this is the sum of read counts for the mature and star miRNAs.';
                $h{6}{1} = 'mature read count(s)';
                $h{6}{2} = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the mature miRNA, including 2 nts upstream and 5 nts downstream. If more than one mature sequence is given this will be a comma separated list. In parenthesis are normalized read counts shown.';
                $h{8}{1} = 'star read count';
                $h{8}{2} = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the star miRNA, including 2 nts upstream and 5 nts downstream. This field is empty unless a reference star miRNA was given as input to qunatify. If more than one mature sequence is given this will be a comman separated list';
                $h{9}{1} = 'remaining reads';
                $h{9}{2} = 'this is the number of reads that did not map to any of the mature and star sequences';
            }
        }

        $h{7}{1} = '-    ';
        $h{7}{2} = '-    ';



        $h{10}{1} = '-'; #'miRBase mature id';
        $h{10}{2} = '-'; #'Clicking this field will link to miRBase.';
        $h{11}{1} = '-';
        $h{11}{2} = '-';
        $h{12}{1} = 'UCSC browser';
        $h{12}{2} = 'if a species name was input to quantify, then clicking this field will initiate a UCSC blat search of the miRNA precursor sequence against the reference genome.';
        $h{13}{1} = 'NCBI blastn';
        $h{13}{2} = 'clicking this field will initiate a NCBI blastn search of the miRNA precursor sequence against the nr/nt database (non-redundant collection of all NCBI nucleotide sequences).';
        if(not $options{'S'}){
            if(not $options{'z'}){
                if($options{'P'}){
                    $h{14}{1} = 'miRBase 5p sequence(s)';
                    $h{14}{2} = 'this is/are the 5p miRNA sequence(s) input to qunatify.';
                    $h{15}{1} = 'miRBase 3p sequence(s)';
                    $h{15}{2} = 'this is/are the 3p miRNA sequence(s) input to qunatify.';
                }else{
                    $h{14}{1} = 'miRBase mature sequence(s)';
                    $h{14}{2} = 'this is/are the mature miRNA sequence(s) input to qunatify.';
                    $h{15}{1} = 'miRBase star sequence(s)';
                    $h{15}{2} = 'this is/are the star miRNA sequence(s) input to qunatify. This field is empty unless a reference star miRNA was given as input to qunatify.';
                }
            }

            $h{16}{1} = 'miRBase precursor sequence';
            $h{16}{2} = 'this is the precursor miRNA sequence input to quantify.';
            if($trna){
                $h{16}{1} = 'tRNA precursor sequence';
                $h{16}{2} = 'this is the precursor tRNA sequence input to quantify.';
            }
            if($options{'a'}){
                $h{17}{1} = 'genomic position';
            }
        }

    }


    if($csv){
        my $f=1;
        print CSV "$hl\n";
        for(sort {$a <=> $b} keys %h){
            if($f){
                $f =0;
                print CSV "$h{$_}{1}";
            }else{
                print CSV "\t$h{$_}{1}";
            }
        }
        print CSV "\n";
        return;
    }


    print HTML <<EOF;

<br>
        <br>
        <br><h2>$hl</h2><br>
        <font face="Times New Roman" size="2">
        <table class="all">
EOF
    for(sort {$a <=> $b} keys %h){

        next if($_ == 2);
        next if($_ == 3);
        next if($_ == 4);
        next if($_ == 7);
        next if($_ == 10);
        next if($_ == 11);
        next if($_ == 8 and $options{'Z'});


        if($_ ne 16){
            if($_ eq 9 and $hl !~ /not/){

                print HTML "<th><a href=\"http://www.ncbi.nlm.nih.gov/entrez/utils/fref.fcgi?PrId=3051&itool=AbstractPlus-def&uid=15217813&nlmid=9808944&db=pubmed&url=http://bioinformatics.oxfordjournals.org/cgi/pmidlookup?view=long&pmid=15217813\" target=\"blank\"class=\"tooltip3\">$h{$_}{1}$p2$h{$_}{2}$q\n\n";

            }else{ 

                if($options{'z'} and ($_ == 5 or $_ ==6)){# we do not have entries for a second item when using options{'z'}	
                    print HTML "$p1$h{$_}{1}\n\n";
                }else{
                    print HTML "$p1$h{$_}{1}$p2$h{$_}{2}$q\n\n";
                }
            }

        }else{
            print HTML "$p11$h{$_}{1}$p2$h{$_}{2}$q\n\n";
        }
    }
}

sub PrintQuantifier{
    my %not_seen; 
    my %signature;
    my $reads;

    ## create HTML for quantifier module
    open HTML,">$exdir/expression_${time}.html" or die "cannot create $exdir/expression_${time}.html\n";
    CreateHTML(); ##
    PrintHtmlTableHeader("",0);

    open UT,">$exdir/precursor_expression_${time}.csv" or die "Could not create file with precursor expression\n";
    open UTB,">$exdir/precursor_expression_${time}_bin.csv" or die "Could not create file with precursor expression\n";


    ## now read in the mature_ref_this_species mapped against precursors from the quantifier module 
    ## store ids in hashes

    my %exprs;
    open IN,"<$exdir/miRNA_expressed.csv" or die "Error: File $exdir/miRNA_expressed.csv not found\n";

    ## this is the replacement for exprs hash
    my %exprs2;

    while(<IN>){
        chomp;
        next if(/precursor/);
        my @line = split(/\t/);

        ## here comes up the trouble when two mature map to same precursor
        $mature2hairpin{$line[0]}=$line[2]; 
        $hairpin2mature{$line[2]}.="$line[0],";
        $exprs{$line[0]} = $line[1];
        $exprs2{$line[2]}{$line[0]} = $line[1];
    }
    close IN;


    my %exprs_sample;
    ## read in sample stuff to expression_value hash
    foreach(@files_mirnaex){
        open IN,"<$_" or die "File $_ not asdfsd found\n";
#        my $sample = $1 if(/_(\S\S\S).csv$/);
        my @samples;
        while(<IN>){
            chomp;
            my @line = split(/\t/);
            if(/precursor/){
                @samples=@line;
            }
            my $sample;
            ## we set precision of print results to one decimal digit
            next if(/^#/);
            for(my $index=4;$index <= $#line;$index++){
                $exprs_sample{$samples[$index]}{$line[2]}{$line[0]} = sprintf("%.f",$line[$index]);
            }
        }
        close IN;
    }


    my $skn=keys %exprs_sample;
    if($skn > 1){
        print UT "#precursor\t#total_reads";
        print UTB "#precursor\t#total_reads";
        for my $sample(sort keys %exprs_sample){
            print UT "\t$sample";
            print UTB "\t$sample";
        }
        for my $sample(sort keys %exprs_sample){
            print UT "\t$sample";
            print UTB "\t$sample";
        }
        print UT "\trest\n";
        print UTB "\trest\n";
    }else{
        print UT "#precursor\t#total_reads\t5p_reads\t3p_reads\trest\n";
        print UTB "#precursor\t#total_reads\t5p_reads\t3p_reads\trest\n";
    }




    ## read in mature,precursor and star seuqnces;


    open IN,"<$exdir/mature.converted" or die "file mature.converted not found\n";

    my ($seq,$id);
    while(<IN>){
        if(/>(\S+)/){
            $id = $1;
            $seq = <IN>;
            chomp $seq;
            $seq  =lc($seq);
            $seq =~ tr/t/u/;

            if($options{'P'}){
                if($id =~ /-5p/){
                    $hm{$id}=$seq;
                }elsif(/-3p/){
                    $hs{$id}=$seq;
                }else{
#					print STDERR "option P is used but no 3p or 5p identifier found in $id";
                    $hs{$id}=$seq;
                    $hm{$id}=$seq;
                }

            }else{
                $hm{$id} = $seq;
            }

#            die "$id\t$seq\n";
        }
    }


    my $width = 0;
    for (keys %hm){                                            #$exprs_sample{$sample}{$line[0]} = $line[1];
        $width  = length($_) if(length($_) > $width);
    }

    $width *=5;

    close IN;
    if(-f "$exdir/star.converted"){
        open IN,"<$exdir/star.converted" or die "file star.converted not found\n";


        while(<IN>){
            if(/>(\S+)/){
                $id = $1;
                $seq = <IN>;
                chomp $seq;
                $seq  =lc($seq);
                $seq =~ tr/t/u/;
                $hs{$id} = $seq;
            }
        }
        close IN;
    }
    open IN,"<$exdir/precursor.converted" or die "file precursor.converted not found\n";
    while(<IN>){
        if(/>(\S+)/){
            $id = $1;
            $seq = <IN>;
            chomp $seq;
            $hp{$id} = $seq;
        }
    }
    close IN;

    open IN,"<$exdir/mature2hairpin" or die "Error: File $exdir/mature2hairpin not found\n";
    my %hairpin2mature2;
    while(<IN>){
        chomp;
        my @line = split(/\t/);
        $hairpin2mature2{$line[0]}=$line[1];
    }
    close IN;


    ## read in the expression_analysis signature.arf of mature_sequences mapped to the precursors to get start and end positions
    open IN,"<$options{'i'}" or die "Error:cannot open $options{'i'} file or option -i not given\n\n";
    while(<IN>){
        chomp;
        my @line = split("\t");
        my $id1h = $line[0]; ## this is the mature ID
        my $id2h = $line[5]; ## this is the precursor ID

        ## remove multiple endings if ambigous just for matching with precursor
        $id1h =~ s/\-5p//g;
        $id1h =~ s/\-3p//g;

        ## here is assumed that multiple precursor ids have 3 - in their id, seems to be ok so far
        if($id2h =~/^(\w+\-\w+\-\w+)\-\d+$/){
            $id2h = $1;
        }
        #print STDERR "$id1h\t$line[5]\t$id2h\t$line[0]\n";
        next if($options{'l'} and $id1h !~ /$id2h/i and $id2h !~ /$id1h/i); ## stringent mapping let7a only allowed to map pre-let7a if k is given

        if($id1h =~ /trna(\d+)/){
            $trna=1;
            my $t1=$1;
            if($id2h =~ /trna(\d+)/){
                my $t2=$1;
                next if($t2 != $t1 and $options{'l'});
            }
        }


        $mat_pre_arf{$line[5]}{$line[0]}{'pos'} = $line[7];
        $mat_pre_arf{$line[5]}{$line[0]}{'end'} = $line[8];
        $mat_pre_arf{$line[5]}{$line[0]}{'seq'} = $line[9];

        ## shift by three if we have the 3p of a tRNA, the first 3 nt are the anticodon
        if($line[0] =~ /trna/ and $line[0] =~ /3p$/){
            $mat_pre_arf{$line[5]}{$line[0]}{'pos'}+=3;
        }
        if($line[0] =~ /\*/){
            $mat_pre_arf{$line[5]}{$line[0]}{'star'} = 1;
            if(not $mat_pre_arf{$line[5]}{'sb'}{'pos'}){
                $mat_pre_arf{$line[5]}{'sb'}{'pos'} = $line[7];
                $mat_pre_arf{$line[5]}{'se'}{'pos'} = $line[8];
            }
        }else{
            $mat_pre_arf{$line[5]}{$line[0]}{'star'} = 0;
            if(not $mat_pre_arf{$line[5]}{'mb'}{'pos'}){
                $mat_pre_arf{$line[5]}{'mb'}{'pos'} = $line[7];
                $mat_pre_arf{$line[5]}{'me'}{'pos'} = $line[8];
            }
        }


        #if($line[0] =~ /\*/){
        #$mat_pre_arf{$line[5]}{$line[0]}{'star'} = 1;
        #	}else{
        #$mat_pre_arf{$line[5]}{$line[0]}{'star'} = 0;
        #}
    }
    close IN;

    ## think about a stringent option that designates if to show only mappings where the ID of precursor and mature is the same !!!
    ## read in star mapped to mature sequence if a star file was given
    if($options{'j'} and -f $options{'j'}){
        open IN,"<$options{'j'}" or die "Error:cannot open $options{'j'} file or option -j not given\n\n";
        while(<IN>){
            chomp;
            my @line = split("\t");
            my $id1h = $line[0]; ## this is the mature ID
            my $id2h = $line[5]; ## this is the precursor ID

            ## remove multiple endings if ambigous just for matching with precursor
            $id1h =~ s/\-5p//g;
            $id1h =~ s/\-3p//g;

            ## here is assumed that multiple precursor ids have 3 - in their id, seems to be ok so far
            #if($id2h =~/^(\w+\-\w+\-\w+)\-\d+$/){
            #    $id2h = $1;
            #}
            next if($options{'l'} and $id1h !~ /$id2h/i and $id2h !~ /$id1h/i); ## stringent mapping let7a only allowed to map pre-let7a if k is given
            if($id1h =~ /trna(\d+)/){
                my $t1=$1;
                if($id2h =~ /trna(\d+)/){
                    my $t2=$1;
                    next if($t2 != $t1 and not $options{'l'});
                }
            }


            # print STDERR "$id2h\t$id1h\t$line[5]\t$line[0]\n";

            $mat_pre_arf{$line[5]}{$line[0]}{'pos'} = $line[7];
            $mat_pre_arf{$line[5]}{$line[0]}{'end'} = $line[8];
            $mat_pre_arf{$line[5]}{$line[0]}{'seq'} = $line[9];
            if($line[0] =~ /\*/){
                $mat_pre_arf{$line[5]}{$line[0]}{'star'} = 1;
                if(not $mat_pre_arf{$line[5]}{'sb'}{'pos'}){
                    $mat_pre_arf{$line[5]}{'sb'}{'pos'} = $line[7];
                    $mat_pre_arf{$line[5]}{'se'}{'pos'} = $line[8];
                }
            }else{
                $mat_pre_arf{$line[5]}{$line[0]}{'star'} = 0;
                if(not $mat_pre_arf{$line[5]}{'mb'}{'pos'}){
                    $mat_pre_arf{$line[5]}{'mb'}{'pos'} = $line[7];
                    $mat_pre_arf{$line[5]}{'me'}{'pos'} = $line[8];
                }
            }
        }
        close IN;
    }


    ## open miRBase.mrd from quantifier module ## precursor ids as hash
    my $oid;
    open IN,"<$options{'q'}" or die "Error: cannot open $options{'q'}\n";
    while(<IN>){
        if(/^\>(\S+)/){
            $id = $1;
            $oid = $1;
            $id =~ tr/|/_/;
            $hash_q{$id}{"oid"}=$oid;
            $hash_q{$id}{"id"}=$id;
            $counter =0;
        }
        elsif(/^remaining read count\s*(\d+)/){
            $hash_q{$id}{'remaining_rc'}=$1;

        }

        elsif(/^total read count\s*(\d*)/){
            $hash_q{$id}{"freq_total"}=$1;
        }

        ## read in here everything that mapped to the precursor
        elsif(/^(\S+) read count\s*(\d+)/){ ## everything else is just read in as id with real name
            $hash_q{$id}{'mapped'}{$1} = $2;
#			print STDERR "yeah $1\n";
        }
        elsif(/^pri_seq\s+(\S+)/){
            my $mseq;
            my $sseq;
            my $sseq_obs;
            $hash_q{$id}{"pri_seq"}=$1;
            my @d;
            if($hash_q{$id}{'obs'}){
                if($hash_q{$id}{'obs'}=~ /(M+)/){
                    my $ind=$-[0];
                    $mseq=substr($hash_q{$id}{"pri_seq"},$ind,length($1));
                }elsif($hash_q{$id}{'obs'}=~ /(5+)/){
                    my $ind=$-[0];
                    $mseq=substr($hash_q{$id}{"pri_seq"},$ind,length($1));
                }elsif($hash_q{$id}{'obs'}=~ /(3+)/){
                    my $ind=$-[0];
                    $mseq=substr($hash_q{$id}{"pri_seq"},$ind,length($1));
                }elsif($hash_q{$id}{'obs'}=~ /(S+)/){
                    my $ind=$-[0];
                    $sseq_obs=substr($hash_q{$id}{"pri_seq"},$ind,length($1));
                }
            }else{
                if(not $hash_q{$id}{'exp'}){

                }elsif($hash_q{$id}{'exp'}=~ /(M+)/){
                    my $ind=$-[0];
                    $mseq=substr($hash_q{$id}{"pri_seq"},$ind,length($1));
                }elsif($hash_q{$id}{'exp'}=~ /(5+)/){
                    my $ind=$-[0];
                    $mseq=substr($hash_q{$id}{"pri_seq"},$ind,length($1));
                }elsif($hash_q{$id}{'exp'}=~ /(3+)/){
                    my $ind=$-[0];
                    $mseq=substr($hash_q{$id}{"pri_seq"},$ind,length($1));

                }elsif($hash_q{$id}{'exp'}=~ /(S+)/){
                    my $ind=$-[0];
                    $sseq=substr($hash_q{$id}{"pri_seq"},$ind,length($1));
                }
            }


            $hash_q{$id}{"mat_seq"} = $mseq;
            $hash_q{$id}{"star_seq"} = $sseq;
            $hash_q{$id}{"star_seq_obs"} = $sseq_obs;
        }
        elsif(/^exp\s+(\S+)/){
            $hash_q{$id}{'exp'}=$1;
        }
        elsif(/^obs\s+(\S+)/){
            $hash_q{$id}{'obs'}=$1;
        }
        elsif(/^pri_struct\s+(\S+)/){
            $hash_q{$id}{"pri_struct"}=$1;
            $reads=1;
            $counter = 0;

            next;
        }
        elsif(/^([a-zA-Z-\d]+)(\S+)\s+(\S+)\s+(\S+)$/ and $reads){  ## do not read in known miRNAs
            $counter++;
            $hash_q{$id}{"reads"}{$1}{$counter}{"rid"} = "$1$2";
            $hash_q{$id}{"reads"}{$1}{$counter}{"seq"} = $3;
            $hash_q{$id}{"reads"}{$1}{$counter}{"mm"}  = $4;
            $hash_q{$id}{"reads"}{$1}{"$1$2"}++;
        }else{}
    }
    close IN;
    print STDERR "parsing miRBase.mrd file finished\n";

    my $sf;
    my $mf;


    #for my $id(sort {$hash_q{$b}{"freq_total"} <=> $hash_q{$a}{"freq_total"}} keys %exprs2){
    for my $id(sort {$hash_q{$b}{"freq_total"} <=> $hash_q{$a}{"freq_total"}} keys %hash_q){
        my %mature = (); ## make hash for mature and star
        my %star = ();
        my $ind;


        my $s_star= $hash_q{$id}{'star_seq'};
        ## now go over everything that is expressed

        for my $k2(keys %{$exprs2{$id}}){ ##k1 == precursor, k2 == mirna mapped to it;

            if($options{'P'}){
                if($k2 =~ /-3p/){
                    $star{$k2} = $exprs2{$id}{$k2};
                }elsif($k2 =~ /-5p/){
                    $mature{$k2} = $exprs2{$id}{$k2};
                }else{
                    $ind= index($hash_q{$id}{"pri_seq"},$hash_q{$id}{"mat_seq"});
#					print STDERR "===  $ind  $id $k2  ", length($hash_q{$id}{"pri_seq"})/2,"\n";

                    if($ind >  length($hash_q{$id}{"pri_seq"})/2 ){
                        $star{$k2} = $exprs2{$id}{$k2};
                        #$mature{$k2} = 0;
                    }elsif($ind >= 0){
                        $mature{$k2} = $exprs2{$id}{$k2};
                        #$star{$k2} = 0;
                    }else{
                        print STDERR "Could not determine where $k2 sits in precursor\nPutting it to the 5p species hash\n";
                        print STDERR "$hash_q{$id}{'pri_seq'}\n$hash_q{$id}{'mat_seq'}\n";
                        $mature{$k2} = $exprs2{$id}{$k2};
                    }
                }
            }else{
                if($k2 =~ /\*/){
                    $star{$k2} = $exprs2{$id}{$k2};
                }else{
                    $mature{$k2} = $exprs2{$id}{$k2};
                }
            }

        }



        $hash_q{$id}{"pdf"} = 1;
        $sig++;
        ## set blat link for pre-miRNA sequence
        if(not $org or $org eq ""){
            $blat = '<td> - </td>';
        }else{
            $blat = "<td><a href=\"http://genome.ucsc.edu/cgi-bin/hgBlat?org=$org&type=BLAT's guess&userSeq=$hash_q{$id}{'pri_seq'}\" target=\"_blank\">blat</a></td>";
        }

        $s_star=$hash_q{$id}{'star_seq_obs'} if($hash_q{$id}{'star_seq_obs'});
        my $s_mat = $hash_q{$id}{'mat_seq'};

        ##here the belonging precursor is shown
        $known="<td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=$id\" target=\"_blank\">$oid</a></td>";

        $sf = $hash_q{$id}{"freq_star"};
        $mf = $hash_q{$id}{"freq_mature"};

        if(not $hash_q{$id}{"freq_mature"} and not $exprs{$oid}){
            $mf=0;
        }else{

            if(not $hash_q{$id}{"freq_mature"} or $hash_q{$id}{"freq_mature"} ne $exprs{$oid} ){
                $mf = $exprs{$oid};
                if($sf){
                    $sf = $hash_q{$id}{"freq_mature"};
                }
            }
        }

        ## do not print precursors with 0 reads mapped to it
        next if(not $hash_q{$id}{"freq_total"});
        ## make links or no links html file depending if option d is given
        if(($options{'p'} and $options{'p'} eq $id) or -f "pdfs_$time/$id.pdf"){
            print HTML "<tr><td nowrap=\"nowrap\"><a href=\"pdfs_$time/$id.pdf\">$id</a></td>";
        }elsif($options{'d'} or not $hash_q{$id}{"freq_total"}){ ## no link to pdf
            print HTML "<tr><td nowrap=\"nowrap\"><div style=\"visibility: hidden\">$id.pdf </div>$id</a></td>";
        }else{
            print HTML "<tr><td nowrap=\"nowrap\"><a href=\"pdfs_$time/$id.pdf\">$id</a></td>";
        }
        print HTML <<EOF;

                        <td><br>$hash_q{$id}{"freq_total"}</td>            
                        <td>

EOF


        ## print to UT file
        print UT "$id\t$hash_q{$id}{'freq_total'}";
        print UTB "$id\t$hash_q{$id}{'freq_total'}";

        $mf = "<table>";
        $s_mat = "<table>";

#$mf .= "<tr><td> miRNA </td><td>total</td>";
        $mf .= "<tr><td WIDTH=$width>sample</td>";		

        my $hmf='';
        my $vmf='';
        for my $sample(sort keys %exprs_sample){                                            #$exprs_sample{$sample}{$line[0]} = $line[1];
            #$mf .= "<td>$sample</td>";
        }

        my ($mkey,$skey);
        my ($mk2,$sk2);
        my $countms=0;
        if(scalar keys %mature != 0){
            for my $mx(keys %mature){ 
                $countms++;
#    $mf .= "<tr><td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=$id\" target=\"_blank\">$mx</a></td><td> $mature{$mx} </td>";
                $vmf .= "\n<tr><td  nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=$id\" target=\"_blank\">$mx</a></td>";
                $mkey='';
                $mk2='';
                for my $sample(sort keys %exprs_sample){                                            #$exprs_sample{$sample}{$line[0]} = $line[1];
                    my @k=keys %{$exprs_sample{$sample}{$id}};
                    my $class='b';
                    foreach my $kk(@k){
                        next if($kk eq $mx);
                        if($exprs_sample{$sample}{$id}{$kk} > $exprs_sample{$sample}{$id}{$mx}){
                            $class="r";
                        }
                    }	
                    if($exprs_sample{$sample}{$id}{$mx} == 0){
                        $class="w";
                    }
                    $mkey.="\t$exprs_sample{$sample}{$id}{$mx}";
                    $mk2.="\t1" if($class eq "b");
                    $mk2.="\t0" if($class eq  "w");
                    $mk2.="\t-1" if($class eq "r");
                    #my $cl='class="zero"';
                    my $cl='';
                    if($exprs_sample{$sample}{$id}{$mx} != 0){
                        $cl='';
                    }
                    $hmf.="<td $cl>$sample</td>" if($countms ==1);
                    $vmf .= "<td $cl><nobr><span class=\"$class\">$exprs_sample{$sample}{$id}{$mx}</span></nobr></td>";
#		die $exprs_sample{$sample}{$id}{$mx};
                }

                $vmf .= "</tr>\n";
                $s_mat .= "<tr><td>$hm{$mx}</td></tr>";
            }
        }else{
            $vmf .= "\n<tr><td  nowrap=\"nowrap\"><a>-</a></td>";
            for my $sample(sort keys %exprs_sample){                                            #$exprs_sample{$sample}{$line[0]} = $line[1];
                $vmf .= "<td><nobr>0</nobr></td>";
                $mkey.="\t0";
                $mk2.="\t0";
#		die $exprs_sample{$sample}{$id}{$mx};
            }

            $vmf .= "</tr>\n";
            $s_mat .= "<tr><td>-</td></tr>";
        }

        ## add rest reads here per sample
        $vmf .= "\n<tr><td  nowrap=\"nowrap\">${id}-rest</td>";
        for my $sample(sort keys %exprs_sample){
            my $idr="${id}-rest";
            $vmf .= "<td><nobr>$exprs_sample{$sample}{$id}{$idr}</nobr></td>";
        }
        $vmf .= "</tr>\n";

        ## fuse them
        $mf.="$hmf$vmf";


        if(not $options{'Z'}){	
            $mf .= "</table>\n";
        }
        $s_mat .= "</table>";

        print HTML "$mf</td>\n";

        if(not $options{'z'}){	
            ## if((scalar keys %star) > 0)
            $s_star = "<table>";
            if(not $options{'Z'}){
                $sf = "<table>";
            }
#$sf .= "<tr><td> miRNA </td><td>total    </td>";
            if(not $options{'Z'}){
                $sf .= "<tr><td WIDTH=$width>sample</td>";		
                for my $sample(sort keys %exprs_sample){                                            #$exprs_sample{$sample}{$line[0]} = $line[1];
                    #$sf .= "<td>$sample</td>";
                }
            }
            my $vsf='';
            my $hsf='';

            my $countss=0;
            if((scalar keys %star) != 0){
                for my $sx(keys %star){
                    $countss++;
#    $sf .= "<tr><td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=$id\" target=\"_blank\">$sx</a></td><td>    $star{$sx} </td>";

                    $vsf .= "\n<tr><td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=$id\" target=\"_blank\">$sx</a></td>";
                    $skey='';	
                    $sk2='';
                    for my $sample(sort keys %exprs_sample){                                            #$exprs_sample{$sample}{$line[0]} = $line[1];
                        my @k=keys %{$exprs_sample{$sample}{$id}};
                        my $class='b';
                        foreach my $kk(@k){
                            next if($kk eq $sx);
                            if($exprs_sample{$sample}{$id}{$kk} > $exprs_sample{$sample}{$id}{$sx}){
                                $class="r";
                            }
                        }	
                        if($exprs_sample{$sample}{$id}{$sx} == 0){
                            $class="w";
                        }
                        $skey.="\t$exprs_sample{$sample}{$id}{$sx}";
                        #my $cl='class="zero"';
                        my $cl='';
                        if($exprs_sample{$sample}{$id}{$sx} != 0){
                            $cl='';
                        }
                        #$hsf.="<td $cl>$sample</td>" if($countss ==1);
                        $vsf .= "<td><nobr><span class=\"$class\">$exprs_sample{$sample}{$id}{$sx}</span></nobr></td>";

                        $sk2.="\t1" if($class eq "b");
                        $sk2.="\t0" if($class eq "w");
                        $sk2.="\t-1" if($class eq "r");

                    }
                    $vsf .="</tr>\n";

                    $s_star .= "<tr><td>$hs{$sx}</td>\n</tr>";
                }
            }else{
                $vsf .= "\n<tr><td  nowrap=\"nowrap\"><a>na</a></td>";
                for my $sample(sort keys %exprs_sample){                                            #$exprs_sample{$sample}{$line[0]} = $line[1];
                    #$vsf .= "<td class=\"zero\"><nobr>0</nobr></td>";
                    $vsf .= "<td><nobr>0</nobr></td>";
                    $skey.="\t0";
                    $sk2.="\t0";
#		die $exprs_sample{$sample}{$id}{$mx};
                }

                $vsf .= "</tr>\n";
                $s_star .= "<tr><td>-</td></tr>";	
            }

            $vsf .= "</table>\n";
            $sf.="$hsf$vsf";
            $s_star .= "</table>";


            #print HTML "<td>$sf</td>\n";
            print HTML "<td>$sf</td>\n";

            if(length($mkey) == 0){$mkey='\t0';}
            if(length($skey) == 0){$skey='\t0';}
            if(length($mk2) == 0){$mk2='\t0';}
            if(length($sk2) == 0){$sk2='\t0';}


            print UT "$mkey$skey";
            print UTB "$mk2$sk2";


            if(not $hash_q{$id}{'remaining_rc'}){
                $hash_q{$id}{'remaining_rc'}=0;
            }
            print UT "\t$hash_q{$id}{'remaining_rc'}\n";
            print UTB "\t$hash_q{$id}{'remaining_rc'}\n";




            print HTML "<td nowrap=\"nowrap\">$hash_q{$id}{'remaining_rc'}</td>\n";
        }
        print HTML <<EOF;
                        $blat
                        <td><a href=${blast}$hash_q{$id}{'pri_seq'}&JOB_TITLE=$hash_q{$id}{"oid"}${blast_query} target="_blank">blast</a></td>
EOF
        if(not $options{'S'}){
            if(not $options{'z'}){
                print HTML "<td>$s_mat</td>\n<td>$s_star</td>\n";
            }
            print HTML "<td>$hash_q{$id}{'pri_seq'}</td>";
        }
        print HTML "</tr>\n";     

    }
    close UT;
    close UTB;
}## end of function


sub Usage{
    print STDERR "\n$script_id version $version
    This is the html and pdf creator script for quantify
    \n[usage]\n\tperl $script_id -q quantify_outfile.mrd -i mature_mapped.arf -M mature_expressed.csv[options]\n\n";
    print STDERR "[options]";
    print STDERR "-i \tarf file with mature sequences mapped to precursor sequences\n";
    print STDERR "-M \tmain quantify expression profile output file\n";
    print STDERR "-t spec\t specify the organism from which your sequencing data was obtained\n";

    print STDERR "-u \t print all available UCSC input organisms\n";
    print STDERR "-o \tsort signature by sample in pdf file, default is by beginning position\n";
    print STDERR "-d \t do not generate pdfs\n\n\n";
    print STDERR "-a \tprint genomic coordinates of mature sequence (still testing)\n";
    print STDERR "-l \tbe stringent when assigning miRNA-precursor connections like mmu-mir only is assigned to mmu-precursor\n";
    print STDERR "-A \tsecondary structure is not drawn in pdf\n";
    print STDERR "-m \tspecies 3 letter code from miRBase\n";
    print STDERR "-X \tdont output alignments in pdf files\n";
    print STDERR "-N \tredo pdf files and overwrite existing ones\n";
    print STDERR "-S \tdont output sequences for e.g. precursors and mature miRNAs\n";
    print STDERR "-Z \tput 5p and 3p arm expression in same table in html file\n";
    print STDERR "-W \toutput also weighed reads if file with those is given as the option argument\n";

    exit;
}

sub insert_image{
    my $img='<img src="data:image/png;base64,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" alt="Quantifier" style="border-style: none" name="quant" title="Quantifier" align=left />';	
    return $img;
    }


    __DATA__
    Human
    Chimp
    Orangutan
    Rhesus
    Marmoset
    Mouse
    Rat
    Guinea Pig
    Cat
    Dog
    Horse
    Cow
    Opossum
    Platypus
    Chicken
    Zebra finch
    Lizard
    X. tropicalis
    Zebrafish
    Tetraodon
    Fugu
    Stickleback
    Medaka
    Lamprey
    Lancelet
    C. intestinalis
    S. purpuratus
    C. elegans
    C. brenneri
    C. briggsae
    C. remanei
    C. japonica
    P. pacificus
    D. melanogaster
    D. simulans
    D. sechellia
    D. yakuba
    D. erecta
    D. ananassae
    D. pseudoobscura
    D. persimilis
    D. virilis
    D. mojavensis
    D. grimshawi
    A. gambiae
    A. mellifera
    S. cerevisiae
