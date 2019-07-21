iterl:=$(shell seq 7)

1/%.h5: data/type1/%.root
	mkdir -p $(dir $@)
	python3 Recon.py 0 $@ $^ > $@.log

define tpl_iter
$(guile (+ $(1) 1))/%.h5: data/type1/%.root $(1)/%.h5
	mkdir -p $$(dir $$@)
	python3 Recon.py $(1) $$@ $$^ > $$@.log
endef

$(foreach i,$(iterl),$(eval $(call tpl_iter,$(i))))

# Delete partial files when the processes are killed.
.DELETE_ON_ERROR:
# Keep intermediate files around
.SECONDARY:
