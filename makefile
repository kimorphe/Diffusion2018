all: pack wet rwk prwk

pack: pack.o grid.o vec2.o set2d.o tcntrl.o
	g++ -o pack pack.o grid.o vec2.o set2d.o tcntrl.o 

wet: wet.o pore.o vec2.o set2d.o tcntrl.o pack.o grid.o 
	g++ -o wet wet.o pore.o vec2.o set2d.o tcntrl.o grid.o
rwk: rwk.o pore.o vec2.o set2d.o tcntrl.o pack.o grid.o 
	g++ -o rwk rwk.o pore.o vec2.o set2d.o tcntrl.o grid.o

prwk: prwk.o pore.o vec2.o set2d.o tcntrl.o pack.o grid.o 
	g++ -o prwk prwk.o pore.o vec2.o set2d.o tcntrl.o grid.o

pore.o: pore.cpp set2d.h vec2.h tcntrl.h grid.h pore.h
	g++ -O2 -c pore.cpp
pack.o: pack.cpp set2d.h tcntrl.h grid.h
	g++ -O2 -c pack.cpp
vec2.o: vec2.cpp vec2.h
	g++ -O2 -c vec2.cpp
set2d.o: set2d.cpp set2d.h vec2.h
	g++ -O2 -c set2d.cpp
tcntrl.o: tcntrl.cpp tcntrl.h
	g++ -O2 -c tcntrl.cpp
grid.o: grid.cpp grid.h tcntrl.h set2d.h
	g++ -O2 -c grid.cpp 
wet.o: wet.cpp tcntrl.h grid.h set2d.h pore.h
	g++ -O2 -c wet.cpp
rwk.o: rwk.cpp tcntrl.h grid.h pore.h
	g++ -O2 -c rwk.cpp
prwk.o: prwk.cpp tcntrl.h grid.h pore.h
	g++ -O2 -c prwk.cpp
