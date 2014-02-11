
CXX := g++
TARGET := generate_input.ex

all: $(TARGET)

$(TARGET): generate_input.cc
	g++ $^ -o $@

.PHONEY: clean
clean:
	rm -f $(TARGET)
