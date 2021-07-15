# build an executable named my2dmhd
all:dcmain.C dcfunctions.C
	g++ -c dcfunctions.C -o dcfunctions.o -O3 -g -Wall -Wextra -pedantic
	g++ -c dcmain.C -o dcmain.o  -O3 -g -Wall -Wextra -pedantic
	g++ dcfunctions.o  dcmain.o -o divmhd -O3 -g -Wall -Wextra -pedantic
	






