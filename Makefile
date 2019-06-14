
1:
	g++ src/genome.cpp -DDATA1 -g -o build/ass.exe

2:
	g++ src/genome.cpp -DDATA2 -g -o build/ass.exe

3:
	g++ src/genome.cpp -DDATA3 -g -o build/ass.exe

4:
	g++ src/genome.cpp -DDATA4 -g -o build/ass.exe


run:
	build/ass.exe

debug:
	gdb build/ass.exe

fun:
	g++ src/havefun.cpp -o build/fun.exe
	build/fun.exe
