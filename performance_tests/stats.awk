#!/usr/bin/awk

BEGIN {}
{
	for (i = 1;i <= NF; i++) vals[i,NR] = $i;
}
END {
	N    = 10
	jmin = NR - N
	if (jmin < 2) exit; # need at least N entries, first line is names

	printf "\nComparing to average and standard deviation of last %d runs\n", N
	jmax     = NR - 1
	exitCode = 0 # all ok

	for (i=1; i <= NF ; i++) {
		y   = 0;
		val = vals[i, NR]

		sum = 0
		for (j = jmin; j <= jmax; ++j) sum += vals[i, j];
		av = sum/N

		for(j = jmin; j <= jmax; j++) y+=(vals[i,j] - av)^2;
		sd = sqrt(y/(N-1))

		printf "%20s (%8.2f +- %8.2f) : %8.2f", vals[i, 1], av, sd, val
		if (val - av > sd) {
			printf "; *** FAIL! *** \n", names[i];
			exitCode = 1;
		} 
		else {
			printf "; GOOD!\n", names[i];
		}
	}
	exit exitCode
}
