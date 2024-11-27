CFLAGS= -g -O2 -Wall -Wextra -fopenmp #-fsanitize=address
CPPFLAGS= $(shell pkg-config --cflags libjpeg libturbojpeg)
LDLIBS= $(shell pkg-config --libs libjpeg libturbojpeg) -lm

SRCDIR=src
OBJDIR=obj
PREFIX?=/usr/local
EXECUTABLE=program

build: $(OBJDIR) $(OBJECTS) $(EXECUTABLE)

OBJECTS= $(OBJDIR)/image.o $(OBJDIR)/processing.o $(OBJDIR)/color.o $(OBJDIR)/main.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $^ -o $@ $(LDLIBS)

$(OBJDIR):
	@mkdir $(OBJDIR)

$(OBJDIR)/main.o: $(SRCDIR)/main.c $(SRCDIR)/main.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJDIR)/image.o: $(SRCDIR)/image.c $(SRCDIR)/image.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJDIR)/color.o: $(SRCDIR)/color.c $(SRCDIR)/color_private.h $(SRCDIR)/color.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJDIR)/processing.o: $(SRCDIR)/processing.c $(SRCDIR)/processing_private.h $(SRCDIR)/processing.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

.PHONEY: clean

clean:
	@-rm -r $(OBJDIR)
	@-rm $(EXECUTABLE)

