CC = gcc
CFLAGS := -O3 -c
SOURCES=main.c classicrbn.c experiments.c rbnio.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=netDynamics

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ -lm

.c.o:
	$(CC) $(CFLAGS) $< -o $@ -lm

