all: terrain
	
raytracer: raytracer.o
	$(CXX) $(CXXFLAGS) -o shaderbox raytracer.o $(LFLAGS) $(LIBS)
raytracer.o : app.cpp app_raytracer.h main.h def.h util.h util_optics.h material.h light.h intersect.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DAPP_RAYTRACER -c -o raytracer.o $<

egg: egg.o
	$(CXX) $(CXXFLAGS) -o shaderbox egg.o $(LFLAGS) $(LIBS)
egg.o : app.cpp app_egg.h main.h def.h util.h IK.h sdf.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DAPP_EGG -c -o egg.o $<

sdf_ao: sdf_ao.o
	$(CXX) $(CXXFLAGS) -o shaderbox sdf_ao.o $(LFLAGS) $(LIBS)
sdf_ao.o : app.cpp app_sdf_ao.h main.h def.h util.h sdf.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DAPP_SDF_AO -c -o sdf_ao.o $<

terrain: terrain.o
	$(CXX) $(CXXFLAGS) -o shaderbox terrain.o $(LFLAGS) $(LIBS)
terrain.o : app.cpp app_terrain.h def.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DAPP_TERRAIN -c -o terrain.o $<

2d: 2d.o
	$(CXX) $(CXXFLAGS) -o shaderbox 2d.o $(LFLAGS) $(LIBS)
2d.o : app.cpp app_2d.h def.h noise.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DAPP_2D -c -o 2d.o $<

clean:
	rm *.o
	rm shaderbox