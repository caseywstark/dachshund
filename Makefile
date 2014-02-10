# The main targets are dachshund.ex and test
#

# The platform file must define the following:
# CXX : the C++ compiler
# CPPFLAGS : the preprocessor flags
# CXXFLAGS : the c++ compile flags
# LDFLAGS : the link flags
include platform.make

LIB_SRCS = $(wildcard lib/*.cc)
LIB_OBJS = $(patsubst lib/%.cc, lib/%.o, $(LIB_SRCS))
MAIN_OBJS = dachshund/main.o
TARGET = dachshund.ex

TEST_SRCS = $(wildcard test/*.cc)
TEST_OBJS = $(patsubst test/%.cc, test/%.o, $(TEST_SRCS))
TEST_TARGET = test_dachshund.ex

# Targets
all: $(TARGET)

$(TARGET): $(LIB_OBJS) $(MAIN_OBJS)
	$(CXX) $(LDFLAGS) $^ $(LIBS) -o $@

$(TEST_TARGET): $(LIB_OBJS) $(TEST_OBJS)
	$(CXX) $(LDFLAGS) $^ $(LIBS) -o $@

%.o: %.cc
	$(CXX) $(DEFINES) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

test: $(TEST_TARGET)
	./$(TEST_TARGET)

.PHONEY: clean
clean:
	-rm -f $(TARGET) $(TEST_TARGET) $(TEST_OBJS) $(MAIN_OBJS) $(LIB_OBJS)
