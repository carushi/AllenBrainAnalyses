#Makefile
CC = g++
CFLAGS = -O2 
LDFLAGS = -lpthread
#INCLUDES = -I/opt/local/include
TARGET = mousebrain

OBJS = main.o filename.o fileformat.o planefileformat.o genesequence.o getrawsequence.o miRNA.o annotationtree.o structurenode.o getstring.o suffixarray.o burrowswheelertransform.o psoutput.o csvoutput.o
.PHONY: all
all: $(TARGET)
test: $(TARGET)
		./$(TARGET) -1 1 0
$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)
.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
main.o: main.cpp filename.o
filename.o: filename.cpp filename.h planefileformat.o fileformat.o genesequence.o miRNA.o thread.h uniq.h
fileformat.o: fileformat.cpp fileformat.h getstring.o
getstring.o: getstring.cpp getstring.h
planefileformat.o: planefileformat.cpp planefileformat.h annotationtree.o genesequence.o mymutex.h
annotationtree.o: annotationtree.cpp annotationtree.h getstring.o structurenode.o
genesequence.o: genesequence.cpp genesequence.h suffixarray.o getrawsequence.o baseString.h
miRNA.o: miRNA.cpp miRNA.h getstring.o planefileformat.o genesequence.o mymutex.h pvalue.h
suffixarray.o: suffixarray.cpp suffixarray.h burrowswheelertransform.o
burrowswheelertransform.o: burrowswheelertransform.cpp burrowswheelertransform.h getstring.o
getrawsequence.o: getrawsequence.cpp getrawsequence.h getstring.o

.PHONY: clean
clean:
	rm -f *.o mousebrain *.gch *~
