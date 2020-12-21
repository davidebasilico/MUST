LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

fai: MUST_v0.cxx
	g++ MUST_v0.cxx -o MUST_v0 -lMinuit ${INCS} ${LIBS}
clean:
	rm *.o output
