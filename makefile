# Makefile for crustcool.cc 

# common directories
CDIR = c
ODIR = o
# local directories
LOCCDIR = c
LOCODIR = o

# compiler
#CC = xcrun c++ 
##CC = clang 
CC=c++
#CC=icpc
#FORTRAN=ifort
FORTRAN=gfortran -m64 -O3
#FORTRAN=gfortran -m64 -O3
CFLAGS = -O3 -pipe
#CFLAGS = -lm -parallel -fast 

# main code
OBJS = $(LOCODIR)/crustcool.o $(LOCODIR)/crust.o $(ODIR)/root.o $(ODIR)/nr.o $(ODIR)/odeint.o $(ODIR)/eos.o $(ODIR)/spline.o $(LOCODIR)/condegin13.o $(LOCODIR)/eosmag12.o $(LOCODIR)/eos12.o $(LOCODIR)/timer.o $(LOCODIR)/data.o $(LOCODIR)/ns.o
OBJS2 = $(LOCODIR)/ocean.o $(ODIR)/root.o $(ODIR)/nr.o $(ODIR)/odeint.o $(ODIR)/eos.o $(ODIR)/spline.o
OBJS3 = $(LOCODIR)/makegrid.o $(ODIR)/root.o $(ODIR)/nr.o $(ODIR)/odeint.o $(ODIR)/eos.o $(ODIR)/spline.o $(LOCODIR)/condegin13.o $(LOCODIR)/eosmag12.o $(LOCODIR)/eos12.o $(LOCODIR)/envelope.o

crustcool : $(OBJS)
	$(CC) -o crustcool $(OBJS) $(CFLAGS) -lm -lgfortran -lgsl

$(LOCODIR)/crustcool.o : $(LOCCDIR)/crustcool.cc
	$(CC) -c $(LOCCDIR)/crustcool.cc -o $(LOCODIR)/crustcool.o $(CFLAGS)

makegrid : $(OBJS3)
	$(CC) -o makegrid $(OBJS3) $(CFLAGS) -lgfortran -lgsl

$(LOCODIR)/makegrid.o : $(LOCCDIR)/makegrid.cc
	$(CC) -c $(LOCCDIR)/makegrid.cc -o $(LOCODIR)/makegrid.o $(CFLAGS) 

$(LOCODIR)/condegin13.o : $(LOCCDIR)/condegin13.f
	$(FORTRAN) -c $(LOCCDIR)/condegin13.f -o $(LOCODIR)/condegin13.o

$(LOCODIR)/eosmag12.o : $(LOCCDIR)/eosmag12.f
	$(FORTRAN) -c $(LOCCDIR)/eosmag12.f -o $(LOCODIR)/eosmag12.o

$(LOCODIR)/eos12.o : $(LOCCDIR)/eos12.f
	$(FORTRAN) -c $(LOCCDIR)/eos12.f -o $(LOCODIR)/eos12.o

# compile routines from the common directory

$(ODIR)/root.o : $(CDIR)/root.c
	$(CC) -c $(CDIR)/root.c -o $(ODIR)/root.o $(CFLAGS)

$(ODIR)/timer.o : $(CDIR)/timer.c
		$(CC) -c $(CDIR)/timer.c -o $(ODIR)/timer.o $(CFLAGS)

$(ODIR)/data.o : $(CDIR)/data.cc
		$(CC) -c $(CDIR)/data.cc -o $(ODIR)/data.o $(CFLAGS)

$(ODIR)/ns.o : $(CDIR)/ns.c
		$(CC) -c $(CDIR)/ns.c -o $(ODIR)/ns.o $(CFLAGS)

$(ODIR)/eos.o : $(CDIR)/eos.cc
	$(CC) -c $(CDIR)/eos.cc -o $(ODIR)/eos.o $(CFLAGS)

$(ODIR)/envelope.o : $(CDIR)/envelope.cc
	$(CC) -c $(CDIR)/envelope.cc -o $(ODIR)/envelope.o $(CFLAGS)

$(ODIR)/crust.o : $(CDIR)/crust.cc
	$(CC) -c $(CDIR)/crust.cc -o $(ODIR)/crust.o $(CFLAGS)

$(ODIR)/nr.o : $(CDIR)/nr.c
	$(CC) -c $(CDIR)/nr.c -o $(ODIR)/nr.o $(CFLAGS)

$(ODIR)/odeint.o : $(CDIR)/odeint.cc
	$(CC) -c $(CDIR)/odeint.cc -o $(ODIR)/odeint.o $(CFLAGS)

$(ODIR)/spline.o : $(CDIR)/spline.cc
	$(CC) -c $(CDIR)/spline.cc -o $(ODIR)/spline.o $(CFLAGS)

# clean up

clean:
	rm -f $(LOCODIR)/*.o

cleanprecalc:
	rm -f gon_out/precalc*

movie:
#	ffmpeg -qscale 1 -r 20 -b 9600 
	ffmpeg -sameq  -i png/%3d.png movie.mp4

cleanpng:
	rm -f png/*.png
