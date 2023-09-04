.PHONY: barcode_demultiplexer
all: barcode_demultiplexer

barcode_demultiplexer:
	g++ barcode_analyzer.cpp SeqIO.cpp -o barcode_analyzer -O3 -Wall -std=c++17
