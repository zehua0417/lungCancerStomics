UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S), Linux)
	PLATFORM := linux
	PYTHON := /hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/tools/anaconda/envs/stereo/bin/python
else ifeq ($(findstring MINGW, $(UNAME_S)), MINGW)
	PLATFORM := windows
	PYTHON := F:/Program\ Data/condaEnvs/stereo/python
else
	PLATFORM := unknown
	PYTHON := python
endif

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
