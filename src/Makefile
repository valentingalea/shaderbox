# mini help
# $@ is the name of the file to be made.
# $? is the names of the changed dependents.
# $< the name of the related file that caused the action.
# $* the prefix shared by target and dependent files.

all: app

APP = -DAPP_PLANET
SIZE = -DSCR_W8=240 -DSCR_H8=240

CXXFLAGS += -Wall -std=c++17 -DC4DROID -fsingle-precision-constant
CXXFLAGS += -Ofast -march=native -funroll-loops

INCLUDES += -I../src -I../lib/vml/test/SDL_app/SDL-1.2.15/include
LIBS += -lSDL
OUT = shaderbox
	
app: app.o
	$(CXX) $(CXXFLAGS) -o $(OUT) app.o $(LIBS)
app.o: ../lib/vml/test/SDL_app/SDL_app.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(APP) $(SIZE) -c -o app.o $<

clean:
	rm *.o
	rm $(OUT)