all: sdf_ao
	
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

clean:
	rm *.o
	rm shaderbox