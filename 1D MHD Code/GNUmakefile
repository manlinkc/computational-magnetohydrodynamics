# build an executable named my 1dmhd
all:1dmhdmain.C 1dmhdfunctions.C 1dmhdtests.C
	g++ -c 1dmhdtests.C -o 1dmhdtests.o -g -Wall -Wextra -pedantic
	g++ -c 1dmhdfunctions.C -o 1dmhdfunctions.o -g -Wall -Wextra -pedantic
	g++ -c 1dmhdmain.C -o 1dmhdmain.o -g -Wall -Wextra -pedantic
	g++ 1dmhdtests.o  1dmhdfunctions.o   1dmhdmain.o -o my1dmhd -O3 -g -Wall -Wextra -pedantic