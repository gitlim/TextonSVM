CC = g++
CFLAGS = -Wall -g -Wl,--rpath
LDFLAGS = 

CFLAGS += `pkg-config --cflags gtk+-2.0`
CFLAGS += `pkg-config --cflags gtkimageview`
CFLAGS += -I.

LDFLAGS += `pkg-config --libs gtk+-2.0`
LDFLAGS += `pkg-config --libs gtkimageview`
LDFLAGS += -L.

all: imresize.o TEXTONSVM_ds.o TEXTONSVM_mrc.o selector.o painter.o $(TEXTONSVM_BIN_DIR)/input_ui

imresize.o: imresize.c imresize.h
	g++ $(CFLAGS) $(LDFLAGS) imresize.c -c

TEXTONSVM_ds.o: TEXTONSVM_ds.cc TEXTONSVM.h
	g++ $(CFLAGS) $(LDFLAGS) TEXTONSVM_ds.cc -c

TEXTONSVM_mrc.o: TEXTONSVM_mrc.c TEXTONSVM.h TEXTONSVM_ds.o
	g++ $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) TEXTONSVM_mrc.c -c

selector.o: selector.c selector.h TEXTONSVM_mrc.o
	g++ $(CFLAGS) $(LDFLAGS) selector.c -c

painter.o: painter.c painter.h selector.o
	g++ $(CFLAGS) $(LDFLAGS) painter.c -c

$(TEXTONSVM_BIN_DIR)/input_ui: imresize.o TEXTONSVM_ds.o TEXTONSVM_mrc.o selector.o painter.o interactive.c
ifdef TEXTONSVM_BIN_DIR 
	g++ -o $(TEXTONSVM_BIN_DIR)/input_ui interactive.c imresize.o TEXTONSVM_ds.o TEXTONSVM_mrc.o selector.o painter.o $(CFLAGS) $(LDFLAGS) -ltiff
else
	@echo ERROR: TEXTONSVM_BIN_DIR is not initialized. Check [TEXTONSVM_DIR/textonsvm_setup];
	@exit 2;
endif


.PHONY: clean
clean:
	rm -f *.o
	rm -f *~
	rm -f $(TEXTONSVM_BIN_DIR)/input_ui



