all: raytracer

# Compile the binary by calling the compiler with cflags, lflags, and any libs (if defined) and the list of objects.
raytracer: app.o
	$(CXX) $(CXXFLAGS) -o raytracer app.o $(LFLAGS) $(LIBS)

# Get a .o from a .cpp by calling compiler with cflags and includes (if defined)
app.o : app.cpp egg.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<
