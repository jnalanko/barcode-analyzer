.PHONY: barcode_demultiplexer
all: barcode_demultiplexer

barcode_demultiplexer:
	g++ barcode_demultiplexer.cpp SeqIO.cpp -o barcode_demultiplexer -O3 -Wall -std=c++17
