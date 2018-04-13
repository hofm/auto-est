.PHONY: all

SRC_DIR		=	src
BUILD_DIR	=	build
ASM_DIR		=   assembly
TARGET		=   25pt-stencil
LIKWID_PATH	= 	/mnt/opt/likwid-4.3.1/

Q	   	   ?=   @
CC			=	icc
CFLAGS  	=	-qopenmp -O3 -xHost -fno-alias -DLIKWID_PERFMON -I/$(LIKWID_PATH)/include
LFLAGS  	=	-mkl -lrt -pthread -llikwid -L$(LIKWID_PATH)/lib
ASFLAGS		=	-S -masm=intel
OBJ			= $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.c))
ASM     	= $(patsubst $(SRC_DIR)/%.c, $(ASM_DIR)/%.s,$(wildcard $(SRC_DIR)/*.c))


default:	$(TARGET)

all:		default

asm: $(ASM_DIR) $(ASM)

$(BUILD_DIR):
	@mkdir $(BUILD_DIR)

$(ASM_DIR):
	@mkdir $(ASM_DIR)

$(TARGET):	$(BUILD_DIR) $(OBJ)
	@echo "===>  LINK $@"
	$(Q)$(CC) $(CFLAGS) -o $@-$(HOSTNAME) $(OBJ) $(LFLAGS)

$(BUILD_DIR)/%.o:  $(SRC_DIR)/%.c
	@echo "===>  COMPILE  $@"
	$(Q)$(CC) -c $(CFLAGS) $< -o $@

$(BUILD_DIR)/%.o:  $(SRC_DIR)/%.s
	@echo "===>  ASSEMBLE  $@"
	$(Q)$(AS) -c $(CFLAGS) $< -o $@

$(ASM_DIR)/%.s: $(SRC_DIR)/%.c
	@echo "===>  CREATING ASSEMBLY  $@"
	$(Q)$(CC) -c $(CFLAGS) $(ASFLAGS) $< -o $@

clean:
	@echo "===>  CLEAN"
	@rm -rf $(BUILD_DIR) $(ASM_DIR)

distclean: clean
	@echo "===>  DIST CLEAN"
	@rm -f $(TARGET)
