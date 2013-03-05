.PHONY: all
all:
ifdef TEXTONSVM_BIN_DIR
	make -C input
	make -C output
	make -C backend
	make -C "backend/textonsvm"
else
	@echo ERROR: TEXTONSVM_BIN_DIR is not initialized. Check [TEXTONSVM_DIR/textonsvm_setup];
	@exit 2;
endif

.PHONY: clean
clean:
	make -r -C input clean
	make -r -C output clean
	make -r -C backend clean
	make -r -C "backend/textonsvm" clean
