SOURCES := main.cpp
# Objs are all the sources, with .cpp replaced by .o
OBJS := $(SOURCES:.cpp=.o)

all: raytracer

# Compile the binary by calling the compiler with cflags, lflags, and any libs (if defined) and the list of objects.
raytracer: $(OBJS)
	$(CXX) $(CFLAGS) -o raytracer $(OBJS) $(LFLAGS) $(LIBS)

# Get a .o from a .cpp by calling compiler with cflags and includes (if defined)
.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<