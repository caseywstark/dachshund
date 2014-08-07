# The platform file must define the following:
# CXX : the C++ compiler
# CPPFLAGS : the preprocessor flags
# CXXFLAGS : the c++ compile flags
# LDFLAGS : the link flags
include platform.make

CPPFLAGS += -Ilib -Ieigen

LIB_SRCS = $(wildcard lib/*.cc)
LIB_OBJS = $(patsubst %.cc,%.o,$(LIB_SRCS))
LIB_TARGET = libdachshund.a

APP_SRCS = $(wildcard app/*.cc)
APP_OBJS = $(patsubst %.cc,%.o,$(APP_SRCS))
APP_TARGET = dachshund.ex

TEST_SRCS = $(wildcard tests/*.cc)
TEST_OBJS = $(patsubst %.cc,%.o,$(TEST_SRCS))
TEST_TARGET = tests/run_tests.ex

# Targets
all: $(LIB_TARGET) $(APP_TARGET) tests

lib/%.o: lib/%.cc lib/%.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(LIB_TARGET): $(LIB_OBJS)
	@echo [MAKE] Archiving static library.
	ar cr $(LIB_TARGET) $(LIB_OBJS)
	ranlib $(LIB_TARGET)

app/%.o: app/%.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(APP_TARGET): $(LIB_TARGET) $(APP_OBJS)
	@echo [MAKE] Linking app.
	$(CXX) $(LDFLAGS) $(APP_OBJS) $(LIB_TARGET) $(LIBS) -o $@

tests/%.o: tests/%.cc
	$(CXX) $(CPPFLAGS) -Itests/UnitTest++/src $(CXXFLAGS) -c $< -o $@

$(TEST_TARGET): $(LIB_TARGET) $(TEST_OBJS)
	@echo [MAKE] Linking test runner.
	$(CXX) $(LDFLAGS) -Ltests/UnitTest++ $(TEST_OBJS) $(LIB_TARGET) -lUnitTest++ $(LIBS) -o $@

.PHONY: tests
tests: $(TEST_TARGET)
	@echo [MAKE] Running tests.
	@./$(TEST_TARGET)

.PHONY: test
test: tests

.PHONY: clean
clean:
	rm -rf $(LIB_OBJS) $(LIB_TARGET) $(APP_OBJS) $(APP_TARGET) $(TEST_OBJS) $(TEST_TARGET)
