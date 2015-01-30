#usage
# perl miRNA_32.pl inputFile outputFile

#input file format (fasta)
#>sequence name
#ribonucleic acid sequence
#second structure

use strict;

my $inputFile             = $ARGV[0];
my $outputFile           = $ARGV[1];

open( IN, "$inputFile" ) or die "can't open the input file : $!";
open( OUT, ">$outputFile" ) or die "can't open the input file : $!";
my ( $name, $seq, $stru);
while (<IN>) {
	chomp $_; # delete the carriage return
	if ( $_ =~ />(.+)/ ) {
		$name = $1;
		$name=~s/[\r\n]$//;		
	}
	$_=<IN>; #read next line
	chomp $_; 
	if ( $_ =~ /([ACUG]+)/ ) {
		$seq = $1;
	}
	$_=<IN>; #read next line
	chomp $_;
	if ( $_ =~ /([\(\)\.]+)/ ) {
		$stru = $1;
	}
#	print OUT ">".$name."\n"; print $seq."\n"; print $stru."\n";
	my @coding_table_32 = sequenceStructure( $seq, $stru );
	my $len_coding_table_32 = @coding_table_32;
	for ( my $i = 0 ; $i < $len_coding_table_32 ; $i++ ) {
		my $rounded = sprintf("%.5f", $coding_table_32[$i]);
		print OUT "$rounded ";
	}
	print OUT "\n";
}

close IN;
close OUT;


sub sequenceStructure {
	my ( $seq_letter, $sec_stru ) = @_;

	#find left start position of base-pair
	my $lpos_bp = index( $sec_stru, "(" );

	#find right start position of base-pair
	my $rpos_bp = rindex( $sec_stru, ")" );

	#delete free base
	$seq_letter = delete_free_base( $seq_letter, $lpos_bp, $rpos_bp );
	$sec_stru   = delete_free_base( $sec_stru,   $lpos_bp, $rpos_bp );

	while ( $sec_stru =~ /(\(\.+\))/g ) {

		#get loop position
		my $pos_loop = index( $sec_stru, $1 );

		#delete loop
		$seq_letter = delete_loop( $seq_letter, $pos_loop, length($1) );
		$sec_stru   = delete_loop( $sec_stru,   $pos_loop, length($1) );
	}

	my @coding_table = translate_to_coding( $seq_letter, $sec_stru );
	unite_probability( \@coding_table );

=pod
	#print feature
	my $len_coding_table = @coding_table;
	for ( my $i = 0 ; $i < $len_coding_table ; $i++ ) {
		print "$coding_table[$i] ";
	}
	print "\n";
=cut

	return @coding_table;
}

sub delete_free_base {
	my ( $seq, $left_pos, $right_pos ) = @_;

	#delete left part
	$seq = substr $seq, $left_pos, length($seq);

	#print "subFun".$seq."\n";
	#delete right part
	$seq = substr $seq, 0, ( $right_pos - $left_pos + 1 );

	#print "subFun".$seq."\n";
	return $seq;
}

sub delete_loop {
	my ( $seq, $pos, $length_loop ) = @_;

	#get first part
	my $seq1 = substr $seq, 0, ( $pos + 1 );

	#get second part
	my $seq2 = substr $seq, ( $pos + $length_loop - 1 ), length($seq);
	$seq = $seq1 . $seq2;

	return $seq;
}

sub get_letter_value {
	my ($char) = @_;
	my ( $l_char, $ret );

	$l_char = lc($char);
	$ret    = -1;
	if ( $l_char eq 'a' ) {
		$ret = 0;
	}
	elsif ( $l_char eq 'g' ) {
		$ret = 1;
	}
	elsif ( $l_char eq 'c' ) {
		$ret = 2;
	}
	elsif ( $l_char eq 'u' ) {
		$ret = 3;
	}
	elsif ( $l_char eq 't' ) {
		$ret = 3;
	}
	else {
		print "ERROR: has not A G C U \n";
	}

	return $ret;
}

sub get_dot_brackle_value {
	my ($char) = @_;
	my ($ret);

	$ret = -1;
	if ( $char eq '.' ) {
		$ret = 0;
	}
	elsif ( $char eq '(' ) {
		$ret = 1;
	}
	elsif ( $char eq ')' ) {
		$ret = 1;
	}
	else {
		print "ERROR: has neither . or ( or )\n";
	}

	return $ret;
}

sub get_struc_value {
	my ($line) = @_;
	my ( $value, $len_line, $i, $char, $value_char );

	$len_line = length $line;
	$value    = 0;
	for ( $i = $len_line - 1 ; $i >= 0 ; $i-- ) {
		$char       = substr( $line, $i, 1 );
		$value_char = get_dot_brackle_value($char);
		$value      = $value + 2**( $len_line - 1 - $i ) * $value_char;
	}

	return $value;
}

sub translate_to_coding {
	my ( $seq_letter, $sec_stru ) = @_;

	my @table;
	for ( my $i = 0 ; $i < 32 ; $i++ ) {
		push( @table, 0 );
	}

	if ( length($seq_letter) != length($sec_stru) ) {
		die "seq_letter is not equal to sec_stru";
	}

	for ( my $i = 0 ; $i < length($seq_letter) ; $i++ ) {
		my $letter = substr( $seq_letter, $i, 1 );
		my $middle = substr( $sec_stru,   $i, 1 );

		my ( $near_left, $near_right, $struc_comb );
		if ( $i == 0 ) {
			$near_left = ".";
			$near_right = substr( $sec_stru, ( $i + 1 ), 1 );
		}
		elsif ( $i == ( length($seq_letter) - 1 ) ) {
			$near_left = substr( $sec_stru, ( $i - 1 ), 1 );
			$near_right = ".";
		}
		else {
			$near_left  = substr( $sec_stru, ( $i - 1 ), 1 );
			$near_right = substr( $sec_stru, ( $i + 1 ), 1 );
		}

		#rectify the empty loop structure
		if ( ( $middle eq "(" ) and ( $near_right eq ")" ) ) {
			$near_right = ".";
		}
		if ( ( $middle eq ")" ) and ( $near_left eq "(" ) ) {
			$near_left = ".";
		}

		#calculate the combin value
		my $letter_value = get_letter_value($letter);
		$struc_comb = $near_left . $middle . $near_right;
		my $struc_value = get_struc_value($struc_comb);
		my $comb_value  = $letter_value * 8 + $struc_value;

		$table[$comb_value]++;
	}
	return @table;
}

sub unite_probability {
	my ($table) = @_;
	my ( $len_table, $i, $sum );

	#	local($show_table);

	$len_table = @$table;
	$sum       = 0;
	for ( $i = 0 ; $i < $len_table ; $i++ ) {
		$sum = $sum + $$table[$i];
	}

	#	print OUT_1 "sum is $sum\n";
	for ( $i = 0 ; $i < $len_table ; $i++ ) {
		$$table[$i] = $$table[$i] / $sum;
	}
}
1;
