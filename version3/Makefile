.PHONY: all

E:= $(shell seq -f '%+.2f' -0.64 0.01 0.64)
files:=$(E:%=file%.h5)

all:$(files)

file%.h5:
	mkdir -p calib
	python3 main_calib.py $* calib/$@

