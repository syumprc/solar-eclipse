.KEEP_STATE:

FC=g77
CC=gcc
CCC=g++

CFLAGS = -O2

INCLUDES= \
safelib.h \
any_list.h \
pipeback.h

OBS= \
safelib.o \
any_list.o \
pipeback.o \
pshinit.o \
error.o

all:	$(INCLUDES) $(OBS) libsafe.a

libsafe.a:	$(INCLUDES) $(OBS)
	        ar -rc libsafe.a any_list.o safelib.o pipeback.o pshinit.o \
				 error.o


install:	libsafe.a
		cp -f libsafe.a /home/charlesp/lib
		cp -f *.h /home/charlesp/lib
		cp -f *.o /home/charlesp/lib

CHARLES:	libsafe.a
		cp -f libsafe.a /export/home/charlesp/lib
		cp -f *.h /export/home/charlesp/lib
clean:
	rm *.o
	rm libsafe.a
