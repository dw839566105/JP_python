srcL:=$(wildcard type1/ele*.root)
dstL:=$(srcL:type1/%.root=output/%.h5)

.PHONY: all
all: $(dstL)

bs=bsub -J $@ -oo $@.log -eo $@.err

output/%.h5: type1/%.root
	mkdir -p output
	$(bs) python3 main_recon.py $@ $^

# Delete partial files when the processes are killed.
.DELETE_ON_ERROR:
# Keep intermediate files around
.SECONDARY:
