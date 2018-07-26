
CXX = g++
CFLAGS = -O3 -Wall

EXE = iprimal
EXPT = expt

OBJS = Assignment.o iPrimal.o Utils.o
EXE_OBJS = $(OBJS) main.o
EXPT_OBJS = $(OBJS) main_expt.o

all: $(EXE)
$(EXE): $(EXE_OBJS)
	$(CXX) $(CFLAGS) -o $@ $(EXE_OBJS)
	@echo Done.

$(EXPT): $(EXPT_OBJS)
	$(CXX) $(CFLAGS) -o $@ $(EXPT_OBJS)
	@echo Done.

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $<

clean:
	rm -f $(EXE) $(EXPT) *.o *~ *.swp


