CC=gcc
# CFLAGS=-gx -DBASE_REPRESENTATION
#CFLAGS=-gx 
CFLAGS=-Wall
LDFLAGS=-ll

YACC= bison  -y
#YFLAGS=-Nl1200 -d -v -t
YFLAGS= -g -k -d -v -t

LEX= flex -s -p

OBJS=y.tab.o lex.yy.o main.o

all: fparse

fparse: $(OBJS)
	$(CC) -o fparse $(CFLAGS) $(OBJS) $(LDFLAGS)

y.tab.o: y.tab.c y.tab.h
	$(CC) -c $(CFLAGS) y.tab.c

y.tab.c y.tab.h: C99-parser.yacc
	$(YACC) $(YFLAGS) C99-parser.yacc

lex.yy.o: lex.yy.c y.tab.h
	$(CC) -c $(CFLAGS) lex.yy.c

lex.yy.c: C99-scanner.lex
	$(LEX) C99-scanner.lex

main.o: main.c
	$(CC) -c $(CFLAGS) main.c

clean:
	rm -f $(OBJS) core y.* lex.yy.? fparse

