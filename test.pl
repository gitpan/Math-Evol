#! /usr/bin/perl
#########################################################################
#        This Perl script is Copyright (c) 2002, Peter J Billam         #
#               c/o P J B Computing, www.pjb.com.au                     #
#                                                                       #
#     This script is free software; you can redistribute it and/or      #
#            modify it under the same terms as Perl itself.             #
#########################################################################

require './Evol.pm'; import Math::Evol;

sub minimise {
	my $sum = 0.0;
	foreach (@_) { $sum += ($_ * $_); }
	return $sum;
}
sub contain {
   if ($_[0] > 1.0) { $_[0] = 1.0;  # it's a greyscale value
   } elsif ($_[0] < 0.0) { $_[0] = 0.0;
   }
   if ($_[1] > 1.0) { $_[1] = 1.0;  # it's a greyscale value
   } elsif ($_[1] < 0.0) { $_[1] = 0.0;
   }
   if ($_[2] > 1.0) { $_[2] = 1.0;  # it's a greyscale value
   } elsif ($_[2] < 0.0) { $_[2] = 0.0;
   }
	return @_;
}
sub choosebetter { my ($a, $b, $c) = @_;
	return ($preference, $continue);
}
my $text = <<'EOT';
/w -3.1416 def % evol step 0.0001 min -3.17 max -3.12
/x 1.4 def % evol step 1
/y 2.17 def % evol step 0.01 min 2.1
/z 4.66 def % evol step 0.1 max 5
EOT

my @x  = (3.456, 1.234, -2.345, 4.567);
my @sm = (.8, .4, .6, 1.2);
my $fail = 0;

my @returns = &evol(\@x, \@sm, \&minimise, \&contain, 10);
my $fail1 = 0;
foreach (@{$returns[0]}) {
	if (abs $_ > 0.0001) { warn "\$_ = $_\n"; $fail++; $fail1++; }
}
if ($fail1) { warn "subroutine &evol failed to find the minimum\n"; }
foreach (@{$returns[1]}) {
	if (abs $_ > 0.0001) { warn "step size still $_\n"; $fail++; $fail1++; }
}
if ($fail1) {
	warn "evol returns:\n x = ", join(", ", @{$returns[0]}), "\n";
	warn "sm = ", join(", ", @{$returns[1]}), "\n";
	warn "objective = $returns[2]\n";
	warn "success   = $returns[3]\n";
} else {
	print "subroutine evol OK\n";
}
if (! $returns[3]) {
	warn "evol ran out of time; maybe you have a slow cpu ?\n";
}

if ($fail) { warn "failed $fail tests\n"; exit 1;
} else { warn "passed all tests\n"; exit 0;
}

&text_evol( $text, \&choosebetter, 2);
exit;

__END__

=pod

=head1 NAME

test.pl - Perl script to test Math::Evol.pm

=head1 SYNOPSIS

=head1 DESCRIPTION

This script tests Math::Evol.pm

=head1 SUBROUTINES

=over 3

=item I<subrname>( $arg1, $arg2 );

=back

=head1 AUTHOR

Peter J Billam <peter@pjb.com.au>

=head1 CREDITS

Based on

=head1 SEE ALSO

http://www.pjb.com.au/, perl(1).

=cut

