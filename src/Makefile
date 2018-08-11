# mini help
# $@ is the name of the file to be made.
# $? is the names of the changed dependents.
# $< the name of the related file that caused the action.
# $* the prefix shared by target and dependent files.

all: planet
CXXFLAGS += -Ofast -march=native -funroll-loops

OUT = shaderbox
SIZE = -DSCREEN_WIDTH=100 -DSCREEN_HEIGHT=100
	
raytracer: raytracer.o
	$(CXX) $(CXXFLAGS) -o $(OUT) raytracer.o $(LFLAGS) $(LIBS)
raytracer.o : app.cpp app_raytracer.h main.h def.h util.h util_optics.h material.h light.h intersect.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DAPP_RAYTRACER $(SIZE) -c -o raytracer.o $<

egg: egg.o
	$(CXX) $(CXXFLAGS) -o $(OUT) egg.o $(LFLAGS) $(LIBS)
egg.o : app.cpp app_egg.h main.h def.h util.h IK.h sdf.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DAPP_EGG $(SIZE) -c -o egg.o $<

sdf_ao: sdf_ao.o
	$(CXX) $(CXXFLAGS) -o $(OUT) sdf_ao.o $(LFLAGS) $(LIBS)
sdf_ao.o : app.cpp app_sdf_ao.h main.h def.h util.h sdf.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DAPP_SDF_AO $(SIZE) -c -o sdf_ao.o $<

clouds: clouds.o
	$(CXX) $(CXXFLAGS) -o $(OUT) clouds.o $(LFLAGS) $(LIBS)
clouds.o : app.cpp app_clouds.h def.h fbm.h util.h volumetric.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DAPP_CLOUDS $(SIZE) -c -o clouds.o $<

atmosphere: atmosphere.o
	$(CXX) $(CXXFLAGS) -o $(OUT) atmosphere.o $(LFLAGS) $(LIBS)
atmosphere.o : app.cpp app_atmosphere.h def.h util.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DAPP_ATMOSPHERE $(SIZE) -c -o atmosphere.o $<

2d: 2d.o
	$(CXX) $(CXXFLAGS) -o $(OUT) 2d.o $(LFLAGS) $(LIBS)
2d.o : app.cpp app_2d.h def.h noise.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DAPP_2D $(SIZE) -c -o 2d.o $<

planet: planet.o
	$(CXX) $(CXXFLAGS) -o $(OUT) planet.o $(LFLAGS) $(LIBS)
planet.o : app.cpp app_planet.h main.h def.h util.h intersect.h noise_iq.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DAPP_PLANET $(SIZE) -c -o planet.o $<

func: func.o
	$(CXX) $(CXXFLAGS) -o $(OUT) func.o $(LFLAGS) $(LIBS)
func.o : app.cpp app_func.h def.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DAPP_FUNC $(SIZE) -c -o func.o $<

vinyl: vinyl.o
	$(CXX) $(CXXFLAGS) -o $(OUT) vinyl.o $(LFLAGS) $(LIBS)
vinyl.o : app.cpp app_vinyl.h main.h def.h util.h sdf.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DAPP_VINYL $(SIZE) -c -o vinyl.o $<


clean:
	rm *.o
	rm $(OUT)