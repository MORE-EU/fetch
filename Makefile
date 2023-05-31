# Name of the module
MODULE_NAME = subann



CC = g++
CFLAGS =  -std=c++14  -fopenmp -O3 -g

SRCDIR = src
OBJDIR = obj
BINDIR = bin

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

TARGET = $(BINDIR)/motifs

MODULE_TARGET = $(BINDIR)/$(MODULE_NAME).so

all: $(TARGET) $(MODULE_TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ -lpython3.8
	
$(MODULE_TARGET): $(OBJECTS)
	$(CC) -shared $(CFLAGS) -o $@ $^ -lpython3.8	

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) -fPIC $(CFLAGS) -c  -o $@ $<



clean:
	rm -f $(OBJDIR)/*.o $(TARGET)

