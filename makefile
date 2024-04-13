CC=g++
# CFLAGS=-gx -DBASE_REPRESENTATION
#CFLAGS=-gx
CFLAGS=-Wall
LDFLAGS=-ll

YACC= bison  -o y.tab.cc
#YFLAGS=-Nl1200 -d -v -t
YFLAGS= -g -k -d -v -t

LEX= flex -s -p -o lex.yy.cc

OBJS=y.tab.o lex.yy.o main.o

all: fparse

fparse: $(OBJS)
	$(CC) -o fparse $(CFLAGS) $(OBJS) $(LDFLAGS)

y.tab.o: y.tab.cc y.tab.h
	$(CC) -c $(CFLAGS) y.tab.cc

y.tab.cc y.tab.h: C99-parser.yacc
	$(YACC) $(YFLAGS) C99-parser.yacc

lex.yy.o: lex.yy.cc y.tab.h
	$(CC) -c $(CFLAGS) lex.yy.cc

lex.yy.cc: C99-scanner.lex
	$(LEX) C99-scanner.lex

main.o: main.cc
	$(CC) -c $(CFLAGS) main.cc

clean:
	rm -f $(OBJS) core y.* lex.yy.? fparse

