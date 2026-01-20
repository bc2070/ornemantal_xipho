<?php
//thin the aims file
$nThinTo = 10e6/1000; //about 1 per kb
$sAIMs = $argv[1]; //"sims/drift0.1/simulated_AIMs_for_AncestryHMM";
$sAIMCounts = $argv[2]; //"sims/drift0.1/simulated_parental_counts_for_AncestryHMM";

$sAIMThinned = $argv[3]; //"sims/drift0.1/simulated_AIMs_for_AncestryHMM_thinned";
$sAIMCountsThinned = $argv[4]; //"sims/drift0.1/simulated_parental_counts_for_AncestryHMM_thinned";

$nLn = exec("cat $sAIMs | wc -l");
$nLn2 = exec("cat $sAIMCounts | wc -l");

if ($nLn != $nLn2) {
	die("Error, line counts differ between $sAIMs and $sAIMCounts\n");
}

$nSampleProb = $nThinTo/$nLn;

if ($nSampleProb >=1) {
	echo("No need to thin, directly copying files\n");
	exec("cp $sAIMs $sAIMThinned; cp $sAIMCounts $sAIMCountsThinned;");
	die();
}

echo("Sample probability: $nSampleProb\n");

$h = fopen($sAIMs, 'r');
$hC = fopen($sAIMCounts, 'r');

$hO = fopen($sAIMThinned, 'w');
$hOC = fopen($sAIMCountsThinned, 'w');

$nRandMax = mt_getrandmax();
$sPrevChr = "";
$nAccumRecRate = 0;

while(true) {
	$sLn1 = fgets($h);
	$sLn2 = fgets($hC);

	if ($sLn1 === false || $sLn2===false) {
		if (!($sLn1 === false && $sLn2===false)) {
			die("Error, files unsynchronized\n");
		}

		break;
	}

	$arrF1 = explode("\t", trim($sLn1));
	$arrF2 = explode("\t", trim($sLn2));
	if ($arrF1[0]!=$arrF2[0] || $arrF1[1]!=$arrF2[1] ) {
		die("Files unsynchronized: $sLn1\n$sLn2\n");
	}

	$bSample = ((mt_rand() / $nRandMax)<$nSampleProb);
	if ($bSample) {
		if ($sPrevChr != $arrF1[0]) {
			fwrite($hO, $sLn1);
			fwrite($hOC, $sLn2);
			$nAccumRecRate = 0;
		} else {

			fwrite($hO, $sLn1);
			$arrF2[6] = strtolower($arrF2[6] + $nAccumRecRate);
			fwrite($hOC, implode("\t", $arrF2)."\n");
			$nAccumRecRate = 0;
		}
	} else {
		$nAccumRecRate += $arrF2[6];
	}

	$sPrevChr = $arrF1[0];
}
?>
