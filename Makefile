# The platform file must define the following:
# CXX : the C++ compiler
# CPPFLAGS : the preprocessor flags
# CXXFLAGS : the c++ compile flags
# LDFLAGS : the link flags
include platform.make

SOURCES=$(wildcard src/*.cc)
OBJECTS=$(patsubst %.cc,%.o,$(SOURCES))
TARGET=dachshund.ex

TEST_SRC=$(wildcard tests/*_tests.cc)
TESTS=$(patsubst %.cc,%,$(TEST_SRC))

# Targets
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@

.PHONY: tests
tests: $(OBJECTS) $(TESTS)
	sh ./tests/runtests.sh

.PHONY: clean
clean:
	rm -rf $(OBJECTS) $(TESTS)
	rm -f tests/tests.log
