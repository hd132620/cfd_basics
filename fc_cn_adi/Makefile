# Makefile
CC = gcc
CFLAGS = -Wall -O2
TARGET = main

SOURCES = main.c setting.c navier_stokes.c poisson.c
OBJECTS = $(SOURCES:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJECTS) -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

setting.o: setting.h
navier_stokes.o : setting.h navier_stokes.h
poisson.o : poisson.h

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) $(OBJECTS)
