#!/usr/bin/perl

use strict;
use warnings;

use File::Basename;
use Getopt::Long qw(:config posix_default no_ignore_case gnu_compat);

$|=1;

# optimize/deoptimize a coding sequence by maximizing/minimizing the mRNA folding energy
# while controling the codon adaptation index and the rare codon count
#
# usage: 
# codon_optimization_HC.pl translation_table codon_frequency in.fa utr.fa outdir 

# settings
our $EnsEnr = "EnsembleEnergy --sequence --silent";

# options
our $CaiThsh = 0.75;
our $RccThsh = 0;
our $AdptWThsh = 0.1;
my $GenSt = 0;
my $GenSp = 11;
my $NTop = 3;
my $TmpDir = "/tmp";
my $Ithr = 0;
my $Nthr = 1;

GetOptions(
  'cai_thsh=f' => \$CaiThsh,
  'rcc_thsh=i' => \$RccThsh,
  'adptw_thsh=f' => \$AdptWThsh,
  'gen_st=i' => \$GenSt,
  'gen_sp=i' => \$GenSp,
  'n_top=i' => \$NTop,
  'tmp_dir=s' => \$TmpDir,
  'ithr=i' => \$Ithr,
  'nthr=i' => \$Nthr,
);

# mandatory
our $TrnTbl = shift;
our $CdnUsg = shift;
my $InFile = shift;
my $UtrFile = shift;
my $OutDir = shift;


print STDERR "loading the translation table\n";
# $AA2Codon->{"M"}->{"ATG"}=1
# $AA2Codon->{"M"}->{"CTG"}=1
# $AA2Codon->{"M"}->{"TTG"}=1
# $Codon2AA->{"ATG"}="M"
my $AA2Codon = {};
my $Codon2AA = {};
&transl_table($AA2Codon, $Codon2AA, $TrnTbl);


print STDERR "loading the codon usage\n";
# $Codon2Usage->{"TTT"}=4.0;
# $Codon2AdptW->{"TTT"}=0.15625;
my $Codon2Usage = {};
my $Codon2AdptW = {};
&codon_usage($Codon2Usage, $CdnUsg);
&codon_adptw($Codon2AdptW, $Codon2Usage, $AA2Codon);


print STDERR "loading the UTR sequence\n";
our ($UtrName, $UtrSeq) = &load_seq($UtrFile);


print STDERR "loading the original sequence\n";
my ($InName, $InSeq) = &load_seq($InFile);


sub load_seq {
    my $infile = shift;

    open(IN, $infile) or die "cannot read $infile\n";
    my $mfa = join("", <IN>);
    close(IN);
    my @name;
    my @seq;
    while ($mfa =~ /(>.+)\n([A-Za-z\n]+)/g){
	my $n = $1;
	my $s = $2;
	$n =~ /^>\s*(.*?)$/;
	$n = $1;
	$s =~ s/\n//g;
	$s = uc($s);
	push(@name, $n);
	push(@seq, $s);
	print STDERR "loaded sequence $n\n";
    }
    @name==@seq or die "wrong file format\n";
    @name==1 or die "multiple fasta file\n";

    return ($name[0], $seq[0]);
}


print STDERR "predicting the original sequence folding free energy\n";
my $ffe = &ensemble_energy($InName, $InSeq, $GenSt, 3*$GenSp, "$TmpDir/$InName.$Ithr.$Nthr.fa");
print STDERR "$InName score $ffe\n";


print STDERR "generating synonymous variants\n";
my $synlist = &synonymous_generate_nt($InSeq, $GenSt+1, $GenSp, $AA2Codon, $Codon2AA);
my $nsyn = @$synlist;
print STDERR "generated $nsyn senonymous sequences\n";


print STDERR "predicting synonymous variants folding free energies\n";
my $ffehash = {};
for (my $i=0; $i<@$synlist; $i++) {
    ($i % $Nthr)==$Ithr or next;

    my $synname = $InName . "_syn" . $i;
    my $synseq = $synlist->[$i];

    my $ok_cai = &codon_adaptation_index_check($synname, $synseq, $GenSt, $GenSp, $Codon2AdptW);
    $ok_cai or next;
    my $ok_rcc = &rare_codon_count_check($synname, $synseq, $GenSt, $GenSp, $Codon2AdptW);
    $ok_rcc or next;

    my $synffe = &ensemble_energy($synname, $synseq, $GenSt, 3*$GenSp, "$TmpDir/$InName.$Ithr.$Nthr.fa");
    $ffehash->{$i} = $synffe;
}
my @sortidx = sort { $ffehash->{$b} <=> $ffehash->{$a} } keys(%$ffehash);


print STDERR "selecting the top $NTop canddidates\n";
for (my $i=0; $i<$NTop; $i++) {
    my $idx = $sortidx[$i];
    my $synname = $InName . "_syn" . $idx;
    my $synseq = $synlist->[$idx];
    my $synffe = $ffehash->{$idx};
    my $synfile = $OutDir . "/" . $synname . ".fa";

    print STDERR "$i candidate $synname score $synffe\n";
    open(OUT, "> $synfile") or die "cannot write $synfile\n";
    print OUT ">", $synname, "\n", $synseq, "\n";
    close(OUT);
}
for (my $i=@sortidx-$NTop; $i<@sortidx; $i++) {
    my $idx = $sortidx[$i];
    my $synname = $InName . "_syn" . $idx;
    my $synseq = $synlist->[$idx];
    my $synffe = $ffehash->{$idx};
    my $synfile = $OutDir . "/" . $synname . ".fa";

    print STDERR "$i candidate $synname score $synffe\n";
    open(OUT, "> $synfile") or die "cannot write $synfile\n";
    print OUT ">", $synname, "\n", $synseq, "\n";
    close(OUT);
}

print STDERR "finished\n";


sub transl_table {
    my ($aa2codon, $codon2aa, $file) = @_;
    my ($aal, $stl, $b1l, $b2l, $b3l) = ("", "", "", "", "");

    open(IN, $file) or die "cannot read $file\n";
    $aal = <IN>;
    $stl = <IN>;
    $b1l = <IN>;
    $b2l = <IN>;
    $b3l = <IN>;
    close(IN);

    $aal =~ /\=\s*(\S+)$/ or die "wrong line format\n";
    $aal = $1;
    $stl =~ /\=\s*(\S+)$/ or die "wrong line format\n";
    $stl = $1;
    $b1l =~ /\=\s*(\S+)$/ or die "wrong line format\n";
    $b1l = $1;
    $b2l =~ /\=\s*(\S+)$/ or die "wrong line format\n";
    $b2l = $1;
    $b3l =~ /\=\s*(\S+)$/ or die "wrong line format\n";
    $b3l = $1;

    my $len = length($aal);
    $len==length($stl) or die "wrong line length\n";
    $len==length($b1l) or die "wrong line length\n";
    $len==length($b2l) or die "wrong line length\n";
    $len==length($b3l) or die "wrong line length\n";

    for (my $i=0; $i<$len; $i++) {
	my $aa = substr($aal, $i, 1);
	my $st = substr($stl, $i, 1);
	my $b1 = substr($b1l, $i, 1);
	my $b2 = substr($b2l, $i, 1);
	my $b3 = substr($b3l, $i, 1);
	my $codon = $b1 . $b2 . $b3;
	$aa2codon->{$aa}->{$codon} = 1;
	$codon2aa->{$codon} = $aa;
	print STDERR "$codon translated to $aa\n";
    }
}


sub codon_usage {
    my ($codon2usage, $file) = @_;

    open(IN, $file) or die "cannot read $file\n";
    while (my $line=<IN>) {
	chomp $line;
	while ($line=~/([ACGU]{3})\s+(\S+?)\(/g) {
	    my $codon = $1;
	    my $usage = $2;
	    $codon =~ s/U/T/g;
	    $codon2usage->{$codon} = $usage;
	    print STDERR "$codon usage $usage\n";
	}
    }
    close(IN);
}


sub codon_adptw {
    my ($codon2adptw, $codon2usage, $aa2codon) = @_;

    foreach my $aa (keys %$aa2codon) {
	my $usage_max = 0.0;
	foreach my $codon (keys %{$aa2codon->{$aa}}) {
	    $usage_max<$codon2usage->{$codon} and $usage_max=$codon2usage->{$codon};
	}
	foreach my $codon (keys %{$aa2codon->{$aa}}) {
	    my $adptw = $codon2usage->{$codon} / $usage_max;
	    $codon2adptw->{$codon} = $adptw;
	    print STDERR "$codon adaptation weight $adptw\n"
	}
    }
}


sub seq_translate {
    my ($seqnt, $codon2aa) = @_; 
    my $seqaa = "";

    for (my $i=0; $i<length($seqnt); $i=$i+3) {
	my $codon = substr($seqnt, $i, 3);
	$seqaa = $seqaa . $codon2aa->{$codon};
    }
    print STDERR "$seqnt translated to $seqaa\n";

    return $seqaa;
}


sub synonymous_generate_aa {
    my ($seqaa, $aa2codon) = @_;
    my $synlist = [];

    my @stack = ();
    push(@stack, "");

    while (@stack) {
	my $syn = pop(@stack);
	if (length($syn)/3==length($seqaa)) {
	    push(@$synlist, $syn);
	}
	else {
	    my $aa = substr($seqaa, length($syn)/3, 1);
	    foreach my $codon (sort keys %{$aa2codon->{$aa}}) {
		push(@stack, $syn . $codon);
	    }
	}
    }

    return $synlist;
}


sub synonymous_generate_nt {
    my ($seqnt, $st, $sp, $aa2codon, $codon2aa) = @_;
    my $synlist = [];
    
    my $seqnt_h = substr($seqnt, 0, 3*$st);
    my $seqnt_m = substr($seqnt, 3*$st, 3*($sp-$st));
    my $seqnt_t = substr($seqnt, 3*$sp);
    my $seqnt2 = $seqnt_h . $seqnt_m . $seqnt_t;
    $seqnt eq $seqnt2 or die "wrong sequence segmentation\n";

    my $seqaa = &seq_translate($seqnt_m, $codon2aa); 
    $synlist = &synonymous_generate_aa($seqaa, $aa2codon);

    for (my $i=0; $i<@$synlist; $i++) {
	$synlist->[$i] = $seqnt_h . $synlist->[$i] . $seqnt_t;
    }

    return $synlist;
}


sub ensemble_energy {
    my ($name, $seq, $st, $sp, $tmpfile) = @_;

    my $sequh = $UtrSeq . substr($seq, $st, ($sp-$st));

    open(OUT, "> $tmpfile") or die "cannot write $tmpfile\n";
    print OUT ">", $name, "\n", $sequh, "\n";
    close(OUT);

    my $command = "$EnsEnr $tmpfile";
    open(IN, "$command | head -n1 | ") or die "cannot execute $command\n";
    my $resline = <IN>;
    close(IN);

    $resline =~ /\s*(\S+)\s*kcal\/mol$/ or die "wrong result line\n";
    my $ee = $1;

    print STDERR "$name ensemble energy $ee\n";

    return $ee;
}


sub codon_adaptation_index_check {
    my ($name, $seq, $st, $sp, $codon2adptw) = @_;

    my $seqm = substr($seq, 3*$st, 3*($sp-$st));

    my $cai = 1.0;
    my $cnt = 0;
    for (my $i=0; $i<length($seqm); $i=$i+3) {
	my $codon = substr($seqm, $i, 3);
	$cai *= $codon2adptw->{$codon};
	$cnt++;
    }
    $cai = $cai ** (1.0 / $cnt);

    my $ok = ($cai<$CaiThsh) ? 0 : 1;
    if ($ok) {
	print STDERR "$name codon adaptation index accepted $cai\n";
    }
    else {
	print STDERR "$name codon adaptation index rejected $cai\n";
    }

    return $ok;
}


sub rare_codon_count_check {
    my ($name, $seq, $st, $sp, $codon2adptw) = @_;

    my $seqm = substr($seq, 3*$st, 3*($sp-$st));

    my @codon_rare = sort {$codon2adptw->{$a} <=> $codon2adptw->{$b}} (keys %$codon2adptw);
    my %codon_cnt;
    foreach my $codon (keys %$codon2adptw) {
	$codon_cnt{$codon} = 0;
    }
    for (my $i=0; $i<length($seqm); $i=$i+3) {
	my $codon = substr($seqm, $i, 3);
	$codon_cnt{$codon}++;
    }

    my $ok = 1;
    my @tag = ();
    for (my $i=0; $i<@codon_rare; $i++) {
	$codon2adptw->{$codon_rare[$i]}<$AdptWThsh or last;
	$codon_cnt{$codon_rare[$i]}>$RccThsh and $ok=0;
	push(@tag, $codon_rare[$i]);
	push(@tag, $codon_cnt{$codon_rare[$i]});
    }

    if ($ok) {
	print STDERR "$name rare codon count accepted ", join(" ", @tag), "\n";
    }
    else {
	print STDERR "$name rare codon count rejected ", join(" ", @tag), "\n";
    }

    return $ok;
}
