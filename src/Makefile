raytracer: raytracer.o
	$(CXX) $(CXXFLAGS) -o raytracer raytracer.o $(LFLAGS) $(LIBS)

raytracer.o : app.cpp app_raytracer.h main.h def.h util.h util_optics.h material.h light.h intersect.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DAPP_RAYTRACER -c -o raytracer.o $<

egg: egg.o
	$(CXX) $(CXXFLAGS) -o egg egg.o $(LFLAGS) $(LIBS)

egg.o : app.cpp app_egg.h main.h def.h util.h IK.h sdf.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DAPP_EGG -c -o egg.o $<