srcL:=$(wildcard type1/ele*.root)
dstL:=$(srcL:type1/%.root=8/%.h5)

.PHONY: all
all: $(dstL)

iterl:=$(shell seq 7)

bs=bsub -J $@ -oo $@.log -eo $@.err

1/%.h5: type1/%.root
	mkdir -p $(dir $@)
	echo $(bs) python3 Recon.py 0 $@ $^

define tpl_iter
$(guile (+ $(1) 1))/%.h5: data/type2/%.root $(1)/%.h5
	mkdir -p $$(dir $$@)
	$$(bs) -w 'done($$(word 2,$$^))' python3 Recon.py $(1) $$@ $$^
endef

$(foreach i,$(iterl),$(eval $(call tpl_iter,$(i))))

# Delete partial files when the processes are killed.
.DELETE_ON_ERROR:
# Keep intermediate files around
.SECONDARY:
