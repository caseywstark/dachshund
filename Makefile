# The platform file must define the following:
# CXX : the C++ compiler
# CPPFLAGS : the preprocessor flags
# CXXFLAGS : the c++ compile flags
# LDFLAGS : the link flags
include platform.make

CPPFLAGS += -I./lib -I./ext -I./ext/gtest -I./ext/gtest/include

LIB_SRCS = $(wildcard lib/*.cc)
LIB_OBJS = $(patsubst %.cc,%.o,$(LIB_SRCS))
LIB_TARGET = libdachshund.a

APP_SRCS = $(wildcard app/*.cc)
APP_OBJS = $(patsubst %.cc,%.o,$(APP_SRCS))
APP_TARGET = dachshund.exe

TEST_SRCS = $(wildcard tests/*.cc)
TEST_OBJS = $(patsubst %.cc,%.o,$(TEST_SRCS))
TEST_TARGET = tests/run_tests.exe

# Targets
all: $(LIB_TARGET) $(APP_TARGET) tests

lib/%.o: lib/%.cc lib/%.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(LIB_TARGET): $(LIB_OBJS)
	@echo [MAKE] Archiving static library.
	ar cr $(LIB_TARGET) $(LIB_OBJS)
	ranlib $(LIB_TARGET)

$(APP_TARGET): $(LIB_TARGET) $(APP_OBJS)
	@echo [MAKE] Linking app.
	$(CXX) $(LDFLAGS) $(APP_OBJS) $(LIB_TARGET) $(LIBS) -o $@

$(TEST_TARGET): $(LIB_TARGET) $(TEST_OBJS)
	@echo [MAKE] Linking test runner.
	$(CXX) $(LDFLAGS) $(TEST_OBJS) $(LIB_TARGET) $(LIBS) -o $@

.PHONY: tests
tests: $(TEST_TARGET)
	@echo [MAKE] Running tests.
	@./$(TEST_TARGET)

.PHONY: clean
clean:
	rm -rf $(LIB_OBJS) $(LIB_TARGET)
	rm -rf $(APP_OBJS) $(APP_TARGET)
	rm -rf $(TEST_OBJS) $(TEST_TARGET)
