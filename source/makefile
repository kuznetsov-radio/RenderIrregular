RenderIrregular	:	Render.o ExtInterface.o ExtMath.o
			g++ $^ -shared -o  RenderIrregular.so -fopenmp
Render.o	:	Render.cpp ExtMath.h Render.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX Render.cpp -fopenmp
ExtInterface.o	:	ExtInterface.cpp Render.h
			g++ -c -std=c++11 -O3 -fPIC -D LINUX ExtInterface.cpp -fopenmp
ExtMath.o	:	ExtMath.cpp
			g++ -c -std=c++11 -O3 -fPIC -D LINUX ExtMath.cpp -fopenmp