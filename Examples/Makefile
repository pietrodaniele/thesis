CC = g++

LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

main.exe : gaussian_fit.cpp
	$(CC) gaussian_fit.cpp -o main.exe ${INCS} ${LIBS}
