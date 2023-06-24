%.obj : %.cpp
	cl /c /EHsc $<

SRCS = grid.cpp grid2.cpp grid3.cpp grid4.cpp

OBJS = $(subst cpp,obj, $(SRCS))

grid.lib : $(OBJS)
	lib /OUT:grid.lib $(OBJS)

getPhi.exe : getPhi.cpp grid.lib
	cl /c /openmp /I.. -EHsc getPhi.cpp
	link getPhi.obj grid.lib agama.lib /LIBPATH:..\x64\release /NOIMPLIB

relax.exe : relax.cpp grid.lib
	cl /c /openmp /I.. -EHsc relax.cpp
	link relax.obj grid.lib agama.lib /LIBPATH:..\x64\release /NOIMPLIB

