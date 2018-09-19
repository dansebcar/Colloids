build:
	g++ -c main.cpp colloid.cpp cluster.cpp
	g++ main.o cluster.o colloid.o -o a -lsfml-graphics -lsfml-window -lsfml-system    
