use ExtUtils::MakeMaker;

WriteMakefile(
	NAME      => 'Math::Evol',
	VERSION_FROM => 'Evol.pm',
	EXE_FILES => [ glob "bin/*" ],
	AUTHOR    => 'PJB, Peter Billam, computing@pjb.com.au',
	ABSTRACT  => 'Evolution Search Algorithm',
	dist      => {COMPRESS=>'gzip -9f', SUFFIX => 'gz'},
);
