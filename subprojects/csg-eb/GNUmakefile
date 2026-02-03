###
### Makefile for building libcsg-eb.a with GNU Make.
### Only works with GCC compiler
### Expects to find paths for dependencies defined in: BOOST_HOME, CATCH2_HOME, CGAL_HOME, PEGTL_HOME
###



TARGET ?= libcsg-eb.a
BUILD_DIR ?= ./build
CONFIG_FILE ?= CsgEbConfig.mk

SRCFILENAMES := csg.cpp \
	levelset_2d.cpp \
	levelset_3d.cpp \
	matrix_functions.cpp \
	parser.cpp \
	exception.cpp \
	$(if $(ENABLE_CGAL),cgal_helper_polygon.cpp,) \
	$(if $(ENABLE_CGAL),cgal_helper_polyhedron.cpp,)

SRCS := $(SRCFILENAMES:%=./src/csg/%)
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS := src include
INC_DIRS += $(BOOST_HOME)/include
INC_DIRS += $(CATCH2_HOME)/include
INC_DIRS += $(CGAL_HOME)/include
INC_DIRS += $(PEGTL_HOME)/include

INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP

ifdef $(ENABLE_CGAL)
	CPPFLAGS += -DUSE_CGAL
endif

$(TARGET): $(OBJS)
	@ar rcs $(TARGET) $(OBJS)

$(BUILD_DIR)/$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 -c $< -o $@

$(CONFIG_FILE): $(TARGET)
	@echo "CSGEB_GIT_HASH := $(shell git describe --abbrev=12 --dirty --always --tags)" > $@

.PHONY: install
install: $(TARGET) $(CONFIG_FILE)
	@$(MKDIR_P) $(DESTDIR)/lib $(DESTDIR)/include
	install -m 0755 $< $(DESTDIR)/lib
	install -m 0755 include/*.hpp $(DESTDIR)/include
	install -m 0755 $(CONFIG_FILE) $(DESTDIR)

.PHONY: clean
clean:
	$(RM) -r $(BUILD_DIR)


-include $(DEPS)

MKDIR_P ?= mkdir -p
