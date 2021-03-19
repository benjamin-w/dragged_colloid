# makefile for repelling_line.c

CC = gcc
LDFLAGS = -lm -lgsl -lgslcblas 
OPTIM = -O3

ALL_SRCS=$(wildcard *.c)
SRCS = $(filter-out mt199337.c, $(ALL_SRCS))
PROGS = $(patsubst %.c,%,$(SRCS))

all: $(PROGS)

%: %.c

	$(CC) $(CFLAGS) $(LDFLAGS) $(OPTIM) -o $@ $<
