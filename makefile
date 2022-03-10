INC_DIR=include
APPS_DIR=app
DEPS_DIR=include
SRCS_DIR=src
OBJS_DIR=obj

CXX=g++
FLAGS= -lm -lgsl -g -lgslcblas -Wall -I $(INC_DIR)

APPS = $(wildcard $(APPS_DIR)/*.cpp)
SRCS = $(wildcard $(SRCS_DIR)/*.cpp)
OBJS = $(patsubst $(SRCS_DIR)/%.cpp,$(OBJS_DIR)/%.o,$(SRCS))

TARGETS = $(if $(ONLY),$(ONLY),$(patsubst $(APPS_DIR)/%.cpp,%,$(APPS)))

all: $(TARGETS)

help:
	@echo "To compile all scripts in the directory $(APPS_DIR), enter simply with make"
	@echo "Or to compile only one ($(TARGETS)), enter with ONLY=name"

$(TARGETS): %: $(OBJS) $(APPS_DIR)/%.cpp 
	$(CXX) -o $@ $(OBJS) $(APPS_DIR)/$@.cpp $(FLAGS)

$(OBJS_DIR)/%.o: $(SRCS_DIR)/%.cpp $(DEPS_DIR)/%.h | $(OBJS_DIR)
	$(CXX) -c -o $@ $< $(FLAGS)

$(OBJS_DIR):
	mkdir $(OBJS_DIR)

clear:
	rm -r $(OBJS_DIR)

.PHONY: all clean help
