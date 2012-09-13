LIBS_PATH = ../libs
STRING_BUF_PATH = $(LIBS_PATH)/string_buffer
UTILITY_LIB_PATH = $(LIBS_PATH)/utility_lib

ifndef CC
	CC = gcc
endif

ifdef DEBUG
	CFLAGS := -DDEBUG=1 --debug
else
	CFLAGS := -O3
endif

CFLAGS := -Wall -Wextra

INCS = -I $(STRING_BUF_PATH) $(STRING_BUF_PATH)/string_buffer.c \
       -I $(UTILITY_LIB_PATH) $(UTILITY_LIB_PATH)/utility_lib.c

LIBS = -lz

all:
	$(CC) -o density_around $(CFLAGS) $(INCS) density_around.c $(LIBS)

clean:
	if test -e density_around; then rm density_around; fi
	if test -e density_around.dSYM; then rm -r density_around.dSYM; fi
	if test -e density_around.greg; then rm density_around.greg; fi
