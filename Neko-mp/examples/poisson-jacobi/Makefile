FFLAGS  = `pkg-config --cflags neko`
LIBS    = `pkg-config --libs neko`
FC      = `pkg-config --variable=compiler neko`

DEST    = poisson_mp
SRC	= ax_poisson_mp.f90 setup_mp.f90 driver_mp.f90
OBJ	= ${SRC:.f90=.o}

all: $(DEST)

install:

clean:
	-rm -f *.o core *.core $(OBJ) $(DEST) *.mod

$(DEST): ${OBJ}
	$(FC) $(FFLAGS) ${OBJ} -o $@  $(LIBS)

%.o: %.f90
	${FC} ${FFLAGS} -c $<


