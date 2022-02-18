PROGRAM.O = domain.o radiation.o particle.o detector.o pixel.o namelist.o IO.o SIRC.o

M_MAKE = $(MAKE) -j4

FC = mpic++ 
LD = mpic++ 

Optimized = -O3

FCFLAGS = -std=c++11 $(Optimized) $(DEFINES) 
LDFLAGS = 

INC = -I /usr/local/include
LIB = -L /usr/local/lib -lhdf5 -lz

EXE = SIRC

all : program

program: $(EXE)

%.o : %.cpp
		$(FC) $(FCFLAGS) -c -o $@ $< $(INC)

$(EXE): $(PROGRAM.O)
		$(LD) $(LDFLAGS) -o $@ $^ $(LIB)

debug: 

	$(M_MAKE) DEFINES="${DEFINE4ALL} -D_DEBUG"


clean:
		rm -r $(PROGRAM.O) $(EXE)
