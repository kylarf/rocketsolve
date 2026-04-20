# lifted from https://gist.github.com/D3r3k23/b2174dbdc8c256958bf480abc8117ab2

SRC_DIR   := src
BUILD_DIR := build
EXE 	  := $(BUILD_DIR)/rocketsolve
CFLAGS    := -O2 -std=gnu11
LDFLAGS   :=
LDLIBS    := -lm

SRCS := $(shell find $(SRC_DIR) -name '*.c')
OBJS := $(subst $(SRC_DIR), $(BUILD_DIR), $(SRCS:.c=.o))

all : $(OBJS) $(EXE)

$(EXE) : $(OBJS) | $(BUILD_DIR)
	@echo "------ Make $(EXE) ------"
	rm -f $(EXE)
	gcc $(CFLAGS) $(LDFLAGS) $(LDLIBS) -o $(EXE) $(OBJS)

$(BUILD_DIR)/%.o : $(SRC_DIR)/%.c | $(BUILD_DIR)
	@echo "------ Make $(@) ------"
	rm -f $@
	gcc $(CFLAGS) -c -o $@ $<

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

-include $(BUILD_DIR)/*.d

clean:
	rm -rf $(BUILD_DIR)/*
