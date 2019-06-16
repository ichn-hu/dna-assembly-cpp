
1:
	g++ src/main.cpp -DDATA1 -g -o build/ass.exe
	build/ass.exe

2:
	g++ src/main.cpp -DDATA2 -g -o build/ass.exe
	build/ass.exe

3:
	g++ src/main.cpp -DDATA3 -g -o build/ass.exe
	build/ass.exe

4:
	g++ src/main.cpp -DDATA4 -g -o build/ass.exe
	build/ass.exe

d1:
	g++ src/main.cpp -DDATA1 -g -o build/ass.exe
	gdb build/ass.exe

d2:
	g++ src/main.cpp -DDATA2 -g -o build/ass.exe
	gdb build/ass.exe

d3:
	g++ src/main.cpp -DDATA3 -g -o build/ass.exe
	gdb build/ass.exe

d4:
	g++ src/main.cpp -DDATA4 -g -o build/ass.exe
	gdb build/ass.exe


run:
	build/ass.exe

debug:
	gdb build/ass.exe

fun:
	g++ src/havefun.cpp -o build/fun.exe
	build/fun.exe
