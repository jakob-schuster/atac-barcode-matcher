CC = gcc
CFLAGS = -Wall -std=c++17
LDFLAGS = -lstdc++ -lz -lm

CPP_SOURCES = src/ssw/ssw_cpp.cpp src/edlib/src/edlib.cpp src/base.cpp src/seq.cpp src/barcode.cpp src/atac_barcode_matcher.cpp
C_SOURCES = src/ssw/ssw.c
OBJECTS = $(CPP_SOURCES:.cpp=.o) $(C_SOURCES:.c=.o)
TARGET = atac-barcode-matcher

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean:
	@rm -f $(TARGET) $(OBJECTS)