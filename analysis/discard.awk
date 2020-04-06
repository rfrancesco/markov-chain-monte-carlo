#!/usr/bin/awk -f

#Discards first N measures, preserving comments (header/footer)

BEGIN {
	if (N == "") {
		printf("Sample usage: ./thermalize.awk -v N=100 input\n")
		exit 1
	}
	header_len = 0;
}

{
	if (/^#/) {
		print 
		header_len++
	}
	else if ((NR - header_len) > N)
		print
}
