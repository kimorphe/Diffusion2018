all: pack pore

pack: pack.o grid.o vec2.o set2d.o tcntrl.o
	g++ -o pack pack.o grid.o vec2.o set2d.o tcntrl.o 

pore: pore.o vec2.o set2d.o tcntrl.o pack.o grid.o
	g++ -o pore pore.o vec2.o set2d.o tcntrl.o grid.o

pore.o: pore.cpp set2d.h vec2.h tcntrl.h grid.h
	g++ -c pore.cpp
pack.o: pack.cpp set2d.h tcntrl.h grid.h
	g++ -c pack.cpp
vec2.o: vec2.cpp vec2.h
	g++ -c vec2.cpp
set2d.o: set2d.cpp set2d.h vec2.h
	g++ -c set2d.cpp
tcntrl.o: tcntrl.cpp tcntrl.h
	g++ -c tcntrl.cpp
grid.o: grid.cpp grid.h tcntrl.h set2d.h
	g++ -c grid.cpp 
