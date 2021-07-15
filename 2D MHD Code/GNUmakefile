# build an executable named my2dmhd
all:2dmhdmain.C 2dmhdfunctions.C
	g++ -c 2dmhdfunctions.C -o mhdfunctions.o -O3 -g -Wall -Wextra -pedantic
	g++ -c 2dmhdmain.C -o mhdmain.o  -O3 -g -Wall -Wextra -pedantic
	g++ mhdfunctions.o  mhdmain.o -o my2dmhd -O3 -g -Wall -Wextra -pedantic
	






