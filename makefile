UNAME_S := $(shell uname -s)
USERNAME := $(shell whoami)

ifeq ($(UNAME_S), Linux)
	PLATFORM := linux
	ifeq ($(USERNAME), stereonote)
		PYTHON := python
	else
		PYTHON := /hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/tools/anaconda/envs/stereopy-rapids/bin/python
else ifeq ($(findstring MINGW, $(UNAME_S)), MINGW)
	PLATFORM := windows
	PYTHON := F:/Program\ Data/condaEnvs/stereo/python
else
	PLATFORM := unknown
	PYTHON := python
endif

# init dirs
init:
	mkdir -p out/{annotation,hvg,leiden,marker,preprocess,temp,umap}

# ouput platform info
info:
	@echo "Platform: $(PLATFORM)"
	@echo "Python Interpreter: $(PYTHON)"

install:
	$(PYTHON) -m pip install -r requirements.txt

# run tissueAna
tissueAna:
	bash ./misc/bash/run_and_log.sh $(PYTHON) tissueAna

# reload log monitor
reload_log:
	bash ./misc/bash/continue_watch.sh

# use qsub upload job
up_ti_gpu:
	qsub -cwd -l num_proc=32,vf=128G -P P21Z10200N0096 -q st_gpu.q /hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer/misc/bash/tissueAna_job.sh
up_ti:
	qsub -cwd -l num_proc=32,vf=129G -P P21Z10200N0096 -q st.q /hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer/misc/bash/tissueAna_job.sh -o /hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer/log/tissueAna_o.log -e /hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer/log/tissueAna_e.log
up_anno_r:
	qsub -cwd -l num_proc=32,vf=129G -P P21Z10200N0096 -q st.q /hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer/misc/bash/tissueAnno_job.sh -o /hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer/log/tissueAnno_o.log -e /hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer/log/tissueAnno_e.log
