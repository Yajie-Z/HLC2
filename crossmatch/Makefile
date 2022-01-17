# CC          := clang
CXX         := gcc
AR          := ar -rc
CUDA_PATH   ?=/usr/local/cuda-10.2
NVCC        :=$(CUDA_PATH)/bin/nvcc
INCLUDE     := -I/home/zyj/.local/include\
			        -I./
			  
LIBRARIES   :=-L/usr/local/cuda-10.2/lib64 -lcudart -lcufft\
			        -L/home/zyj/.local/lib\

NVCC_FLAGS := -Wno-deprecated-gpu-targets -dc
CXX_FLAGS := -std=c++11 -lstdc++ 
Crossmatch: manager.o helper.o system.o newstruc.o 
	@ echo ./$@ $+
	$(NVCC) -O3 -ccbin $(CXX) -Wno-deprecated-gpu-targets $(INCLUDE) -o $@ $+ $(CUDA_ARCH) $(LIBRARIES) $(CXX_FLAGS)
 
manager.o: manager.cu
	$(NVCC) -O3 -ccbin $(CXX) $(NVCC_FLAGS) $(INCLUDE) -o $@ -c $<
 
helper.o: helper.cu
	$(NVCC) -O3 -ccbin $(CXX) $(NVCC_FLAGS) $(INCLUDE) -o $@ -c $<

system.o: system.cu
	$(NVCC) -O3 -ccbin $(CXX) $(NVCC_FLAGS) $(INCLUDE) -o $@ -c $<
 
newstruc.o: match.cu
	$(NVCC) -O3 -ccbin $(CXX) $(NVCC_FLAGS) $(CXX_FLAGS) $(INCLUDE) -o $@ -c $<


.PHONY : clean
clean:
	rm -rf *.o Crossmatch


