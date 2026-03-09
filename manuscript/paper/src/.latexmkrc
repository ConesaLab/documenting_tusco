$out_dir = '../build';
$aux_dir = '../build';

$ENV{'TEXINPUTS'} = '../styles//:' . ($ENV{'TEXINPUTS'} // '');
$ENV{'BSTINPUTS'} = '../styles//:' . ($ENV{'BSTINPUTS'} // '');
