# lifted from https://gist.github.com/D3r3k23/b2174dbdc8c256958bf480abc8117ab2

SRC_DIR   := src
BUILD_DIR := build
EXE 	  := $(BUILD_DIR)/rocketsolve
GCC_FLAGS := -g -O2 -std=gnu11

SRCS := $(shell find $(SRC_DIR) -name '*.c')
OBJS := $(subst $(SRC_DIR), $(BUILD_DIR), $(SRCS:.c=.o))

all : $(OBJS) $(EXE)

$(EXE) : $(OBJS) | $(BUILD_DIR)
	@echo "------ Make $(EXE) ------"
	rm -f $(EXE)
	gcc $(GCC_FLAGS) -o $(EXE) $(OBJS)

$(BUILD_DIR)/%.o : $(SRC_DIR)/%.c | $(BUILD_DIR)
	@echo "------ Make $(@) ------"
	rm -f $@
	gcc $(GCC_FLAGS) -c -o $@ $<

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

-include $(BUILD_DIR)/*.d

clean:
	rm -rf $(BUILD_DIR)/*
