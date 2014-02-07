# The main targets are dachshund.ex and test
#

# The platform file must define the following:
# INC_PATHS : the include paths
# LIB_PATHS : the library paths
# CXX : the c++ compiler
# CXXFLAGS :
include platform.make

SRCS = $(wildcard src/*.cc)
OBJS = $(patsubst src/%.cc, build/%.o, $(SRCS))
TARGET = dachshund.ex

TEST_SRCS = $(wildcard test/*.cc)
TEST_OBJS = $(patsubst test/%.cc, build/test_%.o, $(TEST_SRCS))
TEST_TARGET = test_dachshund.ex

# Targets
all: $(TARGET)

$(TARGET): build $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

$(TEST_TARGET): $(OBJS) $(TEST_OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) $(TEST_OBJS) $(LIBS) -o $@

test: $(TEST_TARGET)
	./$(TEST_TARGET)

build/%.o: src/%.cc
	$(CXX) $(DEFINES) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

.PHONEY: build
build:
	@mkdir -p build

.PHONEY: clean
clean:
	-rm -rf build $(TARGET) $(TEST_TARGET)
