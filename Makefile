ifdef DEBUG
	FLAGS := -DDEBUG=1 --debug
else
	FLAGS := -O3
endif

STRING_BUF_PATH=../libs/string_buffer
UTILITY_LIB_PATH=../libs/utility_lib

LIB_FILES := -I $(STRING_BUF_PATH) $(STRING_BUF_PATH)/string_buffer.c \
             -I $(UTILITY_LIB_PATH) $(UTILITY_LIB_PATH)/utility_lib.c

all:
	gcc -o density_around $(FLAGS) -Wall -lz $(LIB_FILES) density_around.c

clean:
	if test -e density_around; then rm density_around; fi
	if test -e density_around.dSYM; then rm -r density_around.dSYM; fi
	if test -e density_around.greg; then rm density_around.greg; fi
