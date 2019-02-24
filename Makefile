NAME	= 	MPI

DEBUG	= 	-O3
CC		= 	mpicxx
LD 		=	mpicxx
CFLAGS	=	-c -Wall -Werror -Wshadow -pedantic-errors \
			-Wpointer-arith  \
			-Wmissing-declarations \
			-Wno-long-long -Winline -Wredundant-decls \
			-Wcast-qual -Wcast-align -D__STRICT_ANSI__
LDFLAGS	=	$(DEBUG)
LIBS	=	-lm
SOURCES	=	main.cpp create_matrix.cpp args.cpp func_eval.cpp
OBJECTS	=   $(SOURCES:.cpp=.o)

all: $(SOURCES) $(NAME)
	
$(NAME): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm -rf *.o $(NAME)
