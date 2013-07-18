# Makefile for crustcool.cc 

# common directories
CDIR = c
ODIR = o
# local directories
LOCCDIR = c
LOCODIR = o

# compiler
#CC = /usr/local/bin/g++ 
#CC = clang 
CC=g++
#CC=icpc
#FORTRAN=ifort
FORTRAN=gfortran -m64
CFLAGS = -O3 -pipe
#CFLAGS = -lm -parallel -fast 

# main code
OBJS = $(LOCODIR)/crustcool.o $(ODIR)/root.o $(ODIR)/nr.o $(ODIR)/odeint.o $(ODIR)/eos.o $(ODIR)/spline.o $(ODIR)/special.o $(LOCODIR)/condegin13.o $(LOCODIR)/eosmag12.o $(LOCODIR)/eos12.o
OBJS2 = $(LOCODIR)/ocean.o $(ODIR)/root.o $(ODIR)/nr.o $(ODIR)/odeint.o $(ODIR)/eos.o $(ODIR)/spline.o $(ODIR)/special.o
OBJS3 = $(LOCODIR)/makegrid.o $(ODIR)/root.o $(ODIR)/nr.o $(ODIR)/odeint.o $(ODIR)/eos.o $(ODIR)/spline.o $(ODIR)/special.o $(LOCODIR)/condegin13.o $(LOCODIR)/eosmag12.o $(LOCODIR)/eos12.o

crustcool : $(OBJS)
	$(CC) -o crustcool $(OBJS) $(CFLAGS) -lm -lgfortran

$(LOCODIR)/crustcool.o : $(LOCCDIR)/crustcool.cc
	$(CC) -c $(LOCCDIR)/crustcool.cc -o $(LOCODIR)/crustcool.o $(CFLAGS)

makegrid : $(OBJS3)
	$(CC) -o makegrid $(OBJS3) $(CFLAGS) -lgfortran

$(LOCODIR)/makegrid.o : $(LOCCDIR)/makegrid.cc
	$(CC) -c $(LOCCDIR)/makegrid.cc -o $(LOCODIR)/makegrid.o $(CFLAGS) 

$(LOCODIR)/condegin13.o : $(LOCCDIR)/condegin13.f
	$(FORTRAN) -c $(LOCCDIR)/condegin13.f -o $(LOCODIR)/condegin13.o

$(LOCODIR)/eosmag12.o : $(LOCCDIR)/eosmag12.f
	$(FORTRAN) -c $(LOCCDIR)/eosmag12.f -o $(LOCODIR)/eosmag12.o

$(LOCODIR)/eos12.o : $(LOCCDIR)/eos12.f
	$(FORTRAN) -c $(LOCCDIR)/eos12.f -o $(LOCODIR)/eos12.o

$(LOCODIR)/condegin08.o : $(LOCCDIR)/condegin08.f
	$(FORTRAN) -c $(LOCCDIR)/condegin08.f -o $(LOCODIR)/condegin08.o

$(LOCODIR)/conduct.o : $(LOCCDIR)/conduct.f
	$(FORTRAN) -c $(LOCCDIR)/conduct.f -o $(LOCODIR)/conduct.o


# compile routines from the common directory

$(ODIR)/root.o : $(CDIR)/root.c
	$(CC) -c $(CDIR)/root.c -o $(ODIR)/root.o $(CFLAGS)

$(ODIR)/eos.o : $(CDIR)/eos.cc
	$(CC) -c $(CDIR)/eos.cc -o $(ODIR)/eos.o $(CFLAGS)

$(ODIR)/nr.o : $(CDIR)/nr.c
	$(CC) -c $(CDIR)/nr.c -o $(ODIR)/nr.o $(CFLAGS)

$(ODIR)/odeint.o : $(CDIR)/odeint.cc
	$(CC) -c $(CDIR)/odeint.cc -o $(ODIR)/odeint.o $(CFLAGS)

$(ODIR)/special.o : $(CDIR)/special.cc
	$(CC) -c $(CDIR)/special.cc -o $(ODIR)/special.o $(CFLAGS)

$(ODIR)/spline.o : $(CDIR)/spline.cc
	$(CC) -c $(CDIR)/spline.cc -o $(ODIR)/spline.o $(CFLAGS)

# clean up

cleanall: 
	rm -f $(ODIR)/*.o
	rm -f $(LOCODIR)/*.o

clean:
	rm -f $(LOCODIR)/*.o

cleanprecalc:
	rm -f gon_out/precalc*

movie:
#	ffmpeg -qscale 1 -r 20 -b 9600 
	ffmpeg -qscale 1 -r 20 -b 9600 -i png/%3d.png movie.mp4

cleanpng:
	rm -f png/*.png
