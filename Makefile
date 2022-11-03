CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6dev`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6dev`

chopper_test: chopper_test.o
	$(CC) -o chopper_test chopper_test.o $(LDFLAGS)

chopper_test.o: chopper_test.cpp
	$(CC) -c -o chopper_test.o chopper_test.cpp $(CXXFLAGS)

clean:
	$(RM) *~ *.o chopper_test

