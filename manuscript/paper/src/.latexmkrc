$out_dir = '../build';
$aux_dir = '../build';

$ENV{'TEXINPUTS'} = '../styles//:' . ($ENV{'TEXINPUTS'} // '');
$ENV{'BIBINPUTS'} = '../assets/reference//:' . ($ENV{'BIBINPUTS'} // '');
$ENV{'BSTINPUTS'} = '../styles//:' . ($ENV{'BSTINPUTS'} // '');
