# ==============================================================================
# Makefile for program homo



PROGRAM		= homo
SRCS		= homo.c utils.c
OBJS		= homo.o utils.o
HEADERS         = homo.h utils.h
MAKEFILE	= Makefile

CC		= gcc
CFLAGS		= -ansi -pedantic -Wall -Werror -O
LINKER		= gcc
LDFLAGS		= -O
LIBS		= -lm 

PTR		= 



$(PROGRAM): 	$(OBJS)
		$(LINKER) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
		@echo "program $(PROGRAM) built"

$(OBJS): 	$(HEADERS)

clean:
		rm -f $(OBJS)
		@echo "clean"

reallyclean:
		rm -f $(OBJS) $(PROGRAM)
		@echo "really clean"

print:
		enscript --landscape --columns=2 -M Letter -d $(PTR) $(MAKEFILE) $(HEADERS) $(SRCS)
		@echo "printing $(MAKEFILE) $(HEADERS) $(SRCS)"
		#lp -P$(PTR) $(MAKEFILE) $(HEADERS) $(SRCS)
