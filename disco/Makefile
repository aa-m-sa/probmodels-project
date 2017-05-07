SHELL   = /bin/sh
CXX     = g++
#CXX     = clang++
PACKAGE = BEANDisco
PROGRAM = beand
VERSION = $(shell grep 'define BEAND_VERSION_STRING' beand.cpp | sed 's/.*"\(.*\)"/\1/')
FILES   = $(wildcard *.hpp) $(wildcard dist/*.hpp)
CFLAGS  = -std=c++11 -Wall -Wextra -Werror -Werror=return-type -Wno-unused-parameter -Wno-unused-variable -Wno-comment -Wno-sign-compare
LDFLAGS = -lboost_program_options -lboost_iostreams -lboost_regex
#CFLAGS  = -I$$HOME/.local/include -std=c++0x $(FLAGS)
#LDFLAGS = -L$$HOME/.local/lib -lboost_program_options -lboost_iostreams -lboost_regex

release: CFLAGS += -O2 -DNDEBUG 
release: $(PROGRAM)

debug: CFLAGS += -g
debug: $(PROGRAM)

profiled2: CFLAGS += -g -pg -O2 -DNDEBUG
profiled2: $(PROGRAM)

profiled1: CFLAGS += -g -pg -O1 -DNDEBUG
profiled1: $(PROGRAM)

$(PROGRAM): $(PROGRAM).cpp $(FILES)
	$(CXX) $(CFLAGS) $(PROGRAM).cpp $(LDFLAGS) -o $(PROGRAM)

clean:
	rm -f $(PROGRAM) $(PACKAGE)-$(VERSION).tar.gz gmon.out

package: clean
	tar -cpz --exclude-backups -f $(PACKAGE)-$(VERSION).tar.gz \
			beand.cpp *.hpp dist/*.hpp runtests \
			ChangeLog COPYING Makefile README example.cfg example.dat

