srcL:=$(wildcard ../type4/*.root)
dstL:=$(srcL:%.root=%.h5)

.PHONY: all
all: $(dstL)

%.h5: %.root
	python3 ConvertTruth.py $^ $@
