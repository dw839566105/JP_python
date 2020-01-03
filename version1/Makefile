srcL:=$(wildcard input/data/a/*.root)
dstL:=$(srcL:input/data/%.root=output/%.h5)

.PHONY: all
all: $(dstL)

bs=bsub -J $@ -oo $@.log -eo $@.err

output/%.h5: input/data/%.root
	mkdir -p output/$*
	$(bs) python3 main_recon.py $@ $^

# Delete partial files when the processes are killed.
.DELETE_ON_ERROR:
# Keep intermediate files around
.SECONDARY:
