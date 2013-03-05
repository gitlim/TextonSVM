CC = g++
CFLAGS = -Wall -g 
#-Wl,--rpath
LDFLAGS = 

CFLAGS += `pkg-config --cflags gtk+-2.0`
CFLAGS += `pkg-config --cflags gtkimageview`
LDFLAGS += `pkg-config --libs gtk+-2.0`
LDFLAGS += `pkg-config --libs gtkimageview`

CFLAGS += -I.
LDFLAGS += -L.

all: TEXTONSVM_ds.o TEXTONSVM_mrc.o $(TEXTONSVM_BIN_DIR)/output_ui

TEXTONSVM_ds.o:  TEXTONSVM_ds.c TEXTONSVM.h
	g++ $(CFLAGS) $(LDFLAGS) TEXTONSVM_ds.c -c

TEXTONSVM_mrc.o: TEXTONSVM_mrc.c TEXTONSVM.h
	g++ $(CFLAGS) $(LDFLAGS) TEXTONSVM_mrc.c -c

$(TEXTONSVM_BIN_DIR)/output_ui: TEXTONSVM_ds.o TEXTONSVM_mrc.o interactive.c
ifdef TEXTONSVM_BIN_DIR
	g++ -o $(TEXTONSVM_BIN_DIR)/output_ui interactive.c $(CFLAGS) $(LDFLAGS) TEXTONSVM_ds.o TEXTONSVM_mrc.o $(CFLAGS) $(LDFLAGS)
	cp ./color.map $(TEXTONSVM_BIN_DIR)/.
else
	@echo ERROR: TEXTONSVM_BIN_DIR is not initialized. Check [TEXTONSVM_DIR/textonsvm_setup];
	@exit 2;
endif

clean:
	rm -f *.o
	rm -f *~
	rm -f $(TEXTONSVM_BIN_DIR)/output_ui



