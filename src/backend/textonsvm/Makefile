$(TEXTONSVM_BIN_DIR)/main_final: main_final.o ../bin/pcap/Image.o
ifdef TEXTONSVM_BIN_DIR
	g++ -ansi -Wall -Wno-long-long -Wundef -Wpointer-arith -Wold-style-cast -Woverloaded-virtual -Wsign-promo -pthread -fexceptions -fPIC -O3 -rdynamic -o $(TEXTONSVM_BIN_DIR)/main_final main_final.o ../bin/pcap/Image.o  ../bin/pcap/libsvm/svm.o -iquote "../include" ../bin/concurrent/threads/child_thread.o  ../bin/concurrent/threads/runnable.o  ../bin/concurrent/threads/synchronization/synchronizables/synchronizable.o  ../bin/concurrent/threads/synchronization/synchronizables/unsynchronized.o  ../bin/concurrent/threads/thread.o   ../bin/io/serialization/serial_input_stream.o  ../bin/io/serialization/serial_output_stream.o  ../bin/lang/array.o  ../bin/lang/exceptions/ex_bad_cast.o  ../bin/lang/exceptions/ex_file_error.o  ../bin/lang/exceptions/ex_file_open_error.o  ../bin/lang/exceptions/ex_file_read_error.o  ../bin/lang/exceptions/ex_file_write_error.o  ../bin/lang/exceptions/ex_index_out_of_bounds.o  ../bin/lang/exceptions/ex_invalid_argument.o  ../bin/lang/exceptions/ex_not_found.o  ../bin/lang/exceptions/ex_null_pointer_dereference.o  ../bin/lang/exceptions/exception.o  ../bin/lang/exceptions/throwable.o  ../bin/math/math.o  ../bin/math/matrices/exceptions/ex_matrix_dimension_mismatch.o  ../bin/math/matrices/exceptions/ex_matrix_singular.o  ../bin/math/matrices/matrix.o  ../bin/math/random/generators/rand_gen_uniform.o  ../bin/math/random/sources/mersenne_twister_64.o  ../bin/math/random/sources/rand_source.o  ../bin/math/random/sources/rand_source_64.o  ../bin/math/random/sources/rand_source_default.o  ../bin/math/random/sources/system_entropy.o  ../bin/math/random/util/randperm.o  ../bin/mlearning/clustering/clusterers/abstract/clusterer.o  ../bin/mlearning/clustering/clusterers/abstract/weighted_clusterer.o  ../bin/mlearning/clustering/clusterers/kmeans/basic_clusterer.o  
else
	@echo ERROR: TEXTONSVM_BIN_DIR is not initialized. Check [TEXTONSVM_DIR/textonsvm_setup];
	@exit 2;
endif

$(TEXTONSVM_BIN_DIR)/main_boost: main_boost.o
ifdef TEXTONSVM_BIN_DIR
	g++ -ansi -Wall -Wno-long-long -Wundef -Wpointer-arith -Wold-style-cast -Woverloaded-virtual -Wsign-promo -pthread -fexceptions -fPIC -O3 -rdynamic -ljpeg -lpng -iquote "../include" -o $(TEXTONSVM_BIN_DIR)/main_boost main_boost.o  ../bin/concurrent/threads/child_thread.o  ../bin/concurrent/threads/runnable.o  ../bin/concurrent/threads/synchronization/synchronizables/synchronizable.o  ../bin/concurrent/threads/synchronization/synchronizables/unsynchronized.o  ../bin/concurrent/threads/thread.o  ../bin/io/formats/image/jpeg.o  ../bin/io/formats/image/png.o  ../bin/io/serialization/serial_input_stream.o  ../bin/io/serialization/serial_output_stream.o  ../bin/lang/array.o  ../bin/lang/exceptions/ex_bad_cast.o  ../bin/lang/exceptions/ex_file_error.o  ../bin/lang/exceptions/ex_file_open_error.o  ../bin/lang/exceptions/ex_file_read_error.o  ../bin/lang/exceptions/ex_file_write_error.o  ../bin/lang/exceptions/ex_index_out_of_bounds.o  ../bin/lang/exceptions/ex_invalid_argument.o  ../bin/lang/exceptions/ex_not_found.o  ../bin/lang/exceptions/ex_null_pointer_dereference.o  ../bin/lang/exceptions/exception.o  ../bin/lang/exceptions/throwable.o  ../bin/math/math.o  ../bin/math/matrices/exceptions/ex_matrix_dimension_mismatch.o  ../bin/math/matrices/exceptions/ex_matrix_singular.o  ../bin/math/matrices/matrix.o  ../bin/math/random/generators/rand_gen_uniform.o  ../bin/math/random/sources/mersenne_twister_64.o  ../bin/math/random/sources/rand_source.o  ../bin/math/random/sources/rand_source_64.o  ../bin/math/random/sources/rand_source_default.o  ../bin/math/random/sources/system_entropy.o  ../bin/math/random/util/randperm.o  ../bin/mlearning/clustering/clusterers/abstract/clusterer.o  ../bin/mlearning/clustering/clusterers/abstract/weighted_clusterer.o  ../bin/mlearning/clustering/clusterers/kmeans/basic_clusterer.o  ../bin/pcap/Image.o  ../bin/pcap/libsvm/svm.o
else
	@echo ERROR: TEXTONSVM_BIN_DIR is not initialized. Check [TEXTONSVM_DIR/textonsvm_setup];
	@exit 2;
endif

$(TEXTONSVM_BIN_DIR)/main: main.o
ifdef TEXTONSVM_BIN_DIR
	g++ -ansi -Wall -Wno-long-long -Wundef -Wpointer-arith -Wold-style-cast -Woverloaded-virtual -Wsign-promo -pthread -fexceptions -fPIC -O3 -rdynamic -ljpeg -lpng -iquote "../include" -o $(TEXTONSVM_BIN_DIR)/main main.o  ../bin/concurrent/threads/child_thread.o  ../bin/concurrent/threads/runnable.o  ../bin/concurrent/threads/synchronization/synchronizables/synchronizable.o  ../bin/concurrent/threads/synchronization/synchronizables/unsynchronized.o  ../bin/concurrent/threads/thread.o  ../bin/io/formats/image/jpeg.o  ../bin/io/formats/image/png.o  ../bin/io/serialization/serial_input_stream.o  ../bin/io/serialization/serial_output_stream.o  ../bin/lang/array.o  ../bin/lang/exceptions/ex_bad_cast.o  ../bin/lang/exceptions/ex_file_error.o  ../bin/lang/exceptions/ex_file_open_error.o  ../bin/lang/exceptions/ex_file_read_error.o  ../bin/lang/exceptions/ex_file_write_error.o  ../bin/lang/exceptions/ex_index_out_of_bounds.o  ../bin/lang/exceptions/ex_invalid_argument.o  ../bin/lang/exceptions/ex_not_found.o  ../bin/lang/exceptions/ex_null_pointer_dereference.o  ../bin/lang/exceptions/exception.o  ../bin/lang/exceptions/throwable.o  ../bin/math/math.o  ../bin/math/matrices/exceptions/ex_matrix_dimension_mismatch.o  ../bin/math/matrices/exceptions/ex_matrix_singular.o  ../bin/math/matrices/matrix.o  ../bin/math/random/generators/rand_gen_uniform.o  ../bin/math/random/sources/mersenne_twister_64.o  ../bin/math/random/sources/rand_source.o  ../bin/math/random/sources/rand_source_64.o  ../bin/math/random/sources/rand_source_default.o  ../bin/math/random/sources/system_entropy.o  ../bin/math/random/util/randperm.o  ../bin/mlearning/clustering/clusterers/abstract/clusterer.o  ../bin/mlearning/clustering/clusterers/abstract/weighted_clusterer.o  ../bin/mlearning/clustering/clusterers/kmeans/basic_clusterer.o  ../bin/pcap/Image.o  ../bin/pcap/libsvm/svm.o
else
	@echo ERROR: TEXTONSVM_BIN_DIR is not initialized. Check [TEXTONSVM_DIR/textonsvm_setup];
	@exit 2;
endif

main_final.o: main_final.cc
	g++ -ansi -Wall -Wno-long-long -Wundef -Wpointer-arith -Wold-style-cast -Woverloaded-virtual -Wsign-promo -pthread -fexceptions -fPIC -O3 -rdynamic -iquote "../include" -c -o main_final.o main_final.cc 

main_boost.o: main_boost.cc
	g++ -ansi -Wall -Wno-long-long -Wundef -Wpointer-arith -Wold-style-cast -Woverloaded-virtual -Wsign-promo -pthread -fexceptions -fPIC -O3 -rdynamic -iquote "../include" -c -o main_boost.o main_boost.cc

main.o: main.cc
	g++ -ansi -Wall -Wno-long-long -Wundef -Wpointer-arith -Wold-style-cast -Woverloaded-virtual -Wsign-promo -pthread -fexceptions -fPIC -O3 -rdynamic -iquote "../include" -c -o main.o main.cc



.PHONY: clean
clean:
	rm -f *.o
	rm -f *~
	rm -f $(TEXTONSVM_BIN_DIR)/main